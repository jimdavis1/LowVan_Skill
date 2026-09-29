#!/usr/bin/env python3
"""Shared binning engine for build_collections.py.

Extracted from the per-taxon scripts after the fourth taxon reproduced the
same 180 lines with a different RULES table. A taxon script now declares only
what is specific to it -- the rules, the windows, the chemistry -- and calls
run(). Everything the earlier scripts did is preserved, including every
tracking file, so their output is byte-comparable.

The one capability the earlier template lacked is WINDOW: an EXPLICIT
(min, max) length window that replaces the median-relative RUNT/GIANT test
for a named feature.

Why that had to exist: RUNT/GIANT are computed from the feature's own median,
which assumes the feature is unimodal. Capillovirus MP is not. It has a 320 aa
class and a 463 aa class, 42% identical over 212 residues, both genuinely ORF2.
The pooled median is 320, so 0.7-1.3x median = 224-416 and the ENTIRE 463 aa
class is discarded as "overlong" -- 119 sequences, silently, with every
number self-consistent. This is the Citrivirus CP failure from
annotation-triage.md arriving through the length filter instead of the
string rules.

  Rule: if a feature's length histogram is bimodal and both modes are real,
  it needs a WINDOW. Never let a pooled median arbitrate between two modes.
"""
import collections
import os
import re
import subprocess
import sys
import tempfile

STD = set("ACDEFGHIKLMNPQRSTVWY")
FATAL = set("*")          # a premature stop means the frame is wrong
RUNT, GIANT = 0.70, 1.30
X_SOLE_MAX = 0.01


def compiled(rules):
    return [(k, re.compile(p, re.I)) for k, p in rules]


def classify(anno, rules, notmod):
    a = " ".join(anno.split())
    for rx, why in notmod:
        if re.search(rx, a, re.I):
            return None, why
    for k, rx in rules:
        if rx.search(a):
            return k, None
    return None, None


def accessory_clusters(pool, min_members, min_seq_id, cov, min_med_len=45):
    """Split a pool of (id, anno, seq) into homology groups by clustering.

    For taxa whose accessory proteins are named by MASS. A mass is not a
    homology group -- measured on Crinivirus, "p22" is FOUR non-homologous
    groups at 188/191/192/193 aa, while "p59" and "p60" are ONE group of 78.
    Binning those on the string gives features that are neither, so the
    string is used only to NAME a cluster, never to form one.

    Naming: the cluster's most common label, reduced to P<mass> when that
    label is a mass and to an uppercase slug otherwise. Two clusters that
    want the same name get A/B/C suffixes in descending member count, which
    is the Ampelovirus P20A/P20B precedent made general.
    """
    if not pool:
        return {}
    t = tempfile.mkdtemp()
    with open(t + "/q.faa", "w") as fh:
        for i, row in enumerate(pool):
            fh.write(">%d\n%s\n" % (i, row[2]))
    subprocess.run(["mmseqs", "easy-cluster", t + "/q.faa", t + "/c", t + "/tmp",
                    "--min-seq-id", str(min_seq_id), "-c", str(cov),
                    "--cov-mode", "1", "--threads", "4"],
                   capture_output=True)
    groups = collections.defaultdict(list)
    for l in open(t + "/c_cluster.tsv"):
        r, m = l.rstrip("\n").split("\t")
        groups[r].append(int(m))
    #  A cluster must be big enough to build a profile AND long enough for
    #  that profile to be specific. Below ~45 aa a PSSM cannot be given a
    #  bit_cutoff that both fires on its own members and stays silent across
    #  a 9 kb genome, so a 33 aa "p4" is dropped rather than shipped as a
    #  feature that can only be wrong in one of two directions.
    def medlen(g):
        L = sorted(len(pool[i][2]) for i in g)
        return L[len(L) // 2]
    big = [g for g in groups.values()
           if len(g) >= min_members and medlen(g) >= min_med_len]
    big.sort(key=len, reverse=True)

    def base(idxs):
        lab = collections.Counter(pool[i][1] for i in idxs).most_common(1)[0][0]
        m = re.match(r"^p(\d{1,3})\b", lab.strip(), re.I)
        if m:
            return "P" + m.group(1)
        w = re.findall(r"[A-Za-z]{3,}", lab)
        return (w[0].upper()[:12] if w else
                re.sub(r"[^A-Z0-9]", "", lab.upper())[:12]) or "UNC"

    want = collections.defaultdict(list)
    for g in big:
        want[base(g)].append(g)
    out = {}
    for name, gs in want.items():
        if len(gs) == 1:
            out[name] = gs[0]
        else:
            for n, g in enumerate(gs):
                out["%s%s" % (name, chr(ord("A") + n))] = g
    return {k: [pool[i] for i in v] for k, v in out.items()}


def run(here, T, MODULE, RULES, NOT_MODELLED=(), LENGTH_ROUTE=None,
        PRE_FLOOR=None, WINDOW=None, JUNCTION=None, AMB_MAX=2, ACCESSORY=None):
    LENGTH_ROUTE = LENGTH_ROUTE or {}
    PRE_FLOOR = PRE_FLOOR or {}
    WINDOW = WINDOW or {}
    JUNCTION = JUNCTION or {}
    rules = compiled(RULES)
    notmod = list(NOT_MODELLED)

    md5s = [l.strip() for l in open(os.path.join(here, "%s.uniq.md5" % T))]
    ann = [l.rstrip("\n").split("\t") for l in open(os.path.join(here, "%s.uniq.id_ann" % T))]
    seq = [l.rstrip("\n").split("\t") for l in open(os.path.join(here, "%s.uniq.seq" % T))]
    n = len(md5s)
    if not (len(ann) == len(seq) == n):
        sys.exit("uniq.* files are not the same length -- run check_dump.py")

    per_md5 = collections.Counter()
    for l in open(os.path.join(here, "%s.id_md5" % T), errors="replace"):
        p = l.rstrip("\n").split("\t")
        if len(p) >= 2 and p[1].strip():
            per_md5[p[1].strip()] += 1

    out = collections.defaultdict(list)
    unassigned = collections.Counter()
    notmodelled = collections.Counter()
    feat_counts = collections.Counter()
    acc_pool = []
    for i in range(n):
        a = ann[i][1] if len(ann[i]) > 1 else ""
        s = seq[i][1] if len(seq[i]) > 1 else ""
        if not s:
            continue
        key, why = classify(a, rules, notmod)
        nfeat = per_md5.get(md5s[i], 1)
        if key:
            out[key].append((ann[i][0] if ann[i] else md5s[i], a, s))
            feat_counts[key] += nfeat
        elif why:
            notmodelled[(a, why)] += nfeat
        else:
            unassigned[a] += nfeat
            if ACCESSORY and re.search(ACCESSORY["pool"], " ".join(a.split()), re.I):
                acc_pool.append((ann[i][0] if ann[i] else md5s[i], " ".join(a.split()), s, nfeat))

    #  ---- accessory discovery: cluster the mass-named pool into real groups
    acc_named = {}
    if ACCESSORY and acc_pool:
        acc_named = accessory_clusters(acc_pool, ACCESSORY.get("min_members", 12),
                                       ACCESSORY.get("min_seq_id", 0.25),
                                       ACCESSORY.get("cov", 0.5),
                                       ACCESSORY.get("min_med_len", 45))
        for k, rows in acc_named.items():
            out[k] = [(fid, a, sq) for fid, a, sq, _nf in rows]
            feat_counts[k] = sum(nf for _f, _a, _s, nf in rows)
            for _fid, a, _s, nf in rows:
                #  decrement by the FEATURE count, not by one: unassigned was
                #  incremented by nfeat, so decrementing by one leaves a
                #  phantom residue and understates the binned fraction.
                if unassigned.get(a):
                    unassigned[a] -= nf
                    if unassigned[a] <= 0:
                        del unassigned[a]

    #  ---- length routing: strings that name two different proteins
    routed = collections.Counter()
    for k, (lo, hi, why) in LENGTH_ROUTE.items():
        if k not in out:
            continue
        keep = [r for r in out[k] if lo <= len(r[2]) <= hi]
        for r in out[k]:
            if not (lo <= len(r[2]) <= hi):
                notmodelled[(r[1], why)] += 1
        routed[k] = len(out[k]) - len(keep)
        out[k] = keep

    #  ---- ambiguity: any non-standard residue, sole-representative exception
    amb = lambda q: sum(1 for c in q if c not in STD)
    xdropped = collections.Counter(); xfatal = collections.Counter(); xkept = []
    for k in list(out):
        rows = out[k]
        fatal = [r for r in rows if set(r[2]) & FATAL]
        xfatal[k] = len(fatal)
        rows = [r for r in rows if not (set(r[2]) & FATAL)]
        keep = [r for r in rows if amb(r[2]) <= AMB_MAX]
        kept_tax = {r[0].split("|")[1].split(".")[0] for r in keep if "|" in r[0]}
        bytax = collections.defaultdict(list)
        for r in rows:
            t = r[0].split("|")[1].split(".")[0] if "|" in r[0] else "?"
            if t not in kept_tax:
                bytax[t].append(r)
        for t, rws in bytax.items():
            best = min(rws, key=lambda r: amb(r[2]) / len(r[2]))
            if amb(best[2]) / len(best[2]) <= X_SOLE_MAX:
                keep.append(best)
                xkept.append((k, t, best[0], amb(best[2]), len(best[2])))
        xdropped[k] = len(out[k]) - len(keep)
        out[k] = keep

    #  ---- length: explicit WINDOW if given, else median-relative
    dropped = collections.Counter(); overlong = collections.Counter()
    prefloored = collections.Counter(); windowed = {}
    for k in list(out):
        if k in PRE_FLOOR:
            before = len(out[k])
            out[k] = [r for r in out[k] if len(r[2]) >= PRE_FLOOR[k]]
            prefloored[k] = before - len(out[k])
        lens = sorted(len(s) for _, _, s in out[k])
        if not lens:
            continue
        if k in WINDOW:
            lo, hi = WINDOW[k]
            windowed[k] = (lo, hi)
        else:
            med = lens[len(lens) // 2]
            lo, hi = RUNT * med, GIANT * med
        dropped[k] = sum(1 for r in out[k] if len(r[2]) < lo)
        overlong[k] = sum(1 for r in out[k] if len(r[2]) > hi)
        out[k] = [r for r in out[k] if lo <= len(r[2]) <= hi]

    #  ---- cleavage chemistry / terminal residue
    jdropped = collections.Counter(); jdetail = []
    for k in list(out):
        if k not in JUNCTION:
            continue
        fs, ls = JUNCTION[k]
        keep = []
        for r in out[k]:
            s = r[2]
            badn = fs is not None and s[0] not in fs
            badc = ls is not None and s[-1] not in ls
            if badn or badc:
                jdetail.append((k, r[0], len(s), s[0], s[-1],
                                "N" if badn and not badc else "C" if badc and not badn else "NC"))
            else:
                keep.append(r)
        jdropped[k] = len(out[k]) - len(keep)
        out[k] = keep

    cdir = os.path.join(here, "collections", MODULE)
    os.makedirs(cdir, exist_ok=True)
    for k, rows in sorted(out.items()):
        with open(os.path.join(cdir, k + ".fasta"), "w") as fh:
            for fid, _a, s in rows:
                fh.write(">%s\n%s\n" % (fid, s))

    #  Remove collections this run did not produce. Re-running with changed
    #  rules leaves the PREVIOUS rule set's fastas on disk, and every
    #  downstream step globs the directory -- so a feature that the rules no
    #  longer create still gets a PSSM built, still gets declared, and still
    #  fires at annotation time. Renaming SUPPRESSOROF to SUPPRESSOR left
    #  both; so did every accessory cluster that changed suffix.
    stale = sorted(f for f in os.listdir(cdir)
                   if f.endswith(".fasta") and f[:-6] not in out)
    for f in stale:
        os.remove(os.path.join(cdir, f))
        for ext in (".phr", ".pin", ".psq", ".pdb", ".pot", ".ptf", ".pto"):
            side = os.path.join(cdir, f + ext)
            if os.path.exists(side):
                os.remove(side)
    if stale:
        print("  removed %d stale collection(s) from a previous rule set: %s"
              % (len(stale), ", ".join(x[:-6] for x in stale)))

    td = os.path.join(here, "collections")
    def dump(name, header, rowsrc):
        with open(os.path.join(td, name), "w") as fh:
            fh.write(header + "\n")
            for r in rowsrc:
                fh.write("\t".join(str(x) for x in r) + "\n")
    dump("LENGTH_ROUTED.tsv", "feature\tremoved\twindow\treason",
         ((k, routed[k], "%d-%d" % LENGTH_ROUTE[k][:2], LENGTH_ROUTE[k][2])
          for k in sorted(routed)))
    dump("PREFLOOR_DROPPED.tsv", "feature\tdropped\tfloor",
         ((k, prefloored[k], PRE_FLOOR[k]) for k in sorted(prefloored)))
    dump("EXPLICIT_WINDOWS.tsv", "feature\tmin\tmax",
         ((k, windowed[k][0], windowed[k][1]) for k in sorted(windowed)))
    dump("RUNTS_DROPPED.tsv", "feature\tdropped\tthreshold",
         ((k, dropped[k], "explicit" if k in WINDOW else "%d%% of median" % (RUNT * 100))
          for k in sorted(dropped)))
    dump("OVERLONG_DROPPED.tsv", "feature\tdropped\tthreshold",
         ((k, overlong[k], "explicit" if k in WINDOW else "%d%% of median" % (GIANT * 100))
          for k in sorted(overlong)))
    dump("AMBIGUOUS_DROPPED.tsv", "feature\tdropped\tof_which_premature_stop",
         ((k, xdropped[k], xfatal[k]) for k in sorted(xdropped)))
    dump("AMBIGUOUS_KEPT.tsv", "feature\ttaxon\tfeature_id\tn_amb\tlength", sorted(xkept))
    dump("JUNCTION_DROPPED.tsv", "feature\tfeature_id\tlength\tfirst\tlast\tbad_end", sorted(jdetail))
    dump("JUNCTION_SUMMARY.tsv", "feature\tdropped", ((k, jdropped[k]) for k in sorted(jdropped)))
    dump("ACCESSORY_CLUSTERS.tsv", "feature\tmembers\tmin_len\tmed_len\tmax_len\tlabels",
         ((k, len(v), min(len(r[2]) for r in v),
           sorted(len(r[2]) for r in v)[len(v) // 2],
           max(len(r[2]) for r in v),
           "; ".join("%s x%d" % (a, c) for a, c in
                     collections.Counter(r[1] for r in v).most_common(4)))
          for k, v in sorted(acc_named.items())))
    dump("UNASSIGNED_TRACKING.tsv", "features\tannotation",
         ((c, a) for a, c in unassigned.most_common()))
    dump("NOT_MODELLED.tsv", "features\tannotation\treason",
         ((c, a, w) for (a, w), c in notmodelled.most_common()))

    print("  %-9s %8s %9s %8s %8s %8s %8s %9s" %
          ("feature", "uniq", "features", "prefloor", "runt", "overlong", "amb", "junction"))
    tot = 0
    for k in sorted(out):
        print("  %-9s %8d %9d %8d %8d %8d %8d %9d" %
              (k, len(out[k]), feat_counts[k], prefloored[k], dropped[k], overlong[k],
               xdropped[k], jdropped[k]))
        tot += feat_counts[k]
    allf = tot + sum(notmodelled.values()) + sum(unassigned.values())
    print("  %-9s %8s %9d" % ("NOTMOD", "", sum(notmodelled.values())))
    print("  %-9s %8d %9d" % ("UNASSIG", len(unassigned), sum(unassigned.values())))
    print("\n  binned %d of %d features = %.1f%%" % (tot, allf, 100.0 * tot / allf))
    return out
