#!/usr/bin/env python3
"""Find features that are called on the same locus as another feature.

Two profiles that both fire on one ORF produce two features on identical or
nested coordinates. The annotation looks rich and is wrong, and nothing in
the build reports it: each profile recovers its own collection perfectly, the
collections may share no sequences at all, and the call rate of the offending
feature looks excellent -- which is the trap, because a feature that fires
inside a longer one is called on nearly every genome and so EARNS copy_num.

Found on Crinivirus, where it invalidated five of eleven copy_num
assignments:

    RNASE3 inside SUPPRESSOR  97% of calls   99% identical, 100% coverage
    P9     inside HSP90       98%            one protein, three mass names
    P26    inside P28         83%
    P8A    inside HSP90       95%

Run it after run_gto_eval.py, on the scored genomes:

    python3 qc_duplicate_calls.py --eval <M>/eval --json <M>/<M>_Viral_PSSM.json \
            --module <M>

Reports three things, because they need different answers:

  EXACT      two features on identical coordinates. Almost always one
             protein under two names -- merge them, or drop one.
  NESTED     a shorter feature inside a longer one. If the two collections
             are homologous, they are one family; if they are not, the
             shorter profile is firing promiscuously and needs a higher
             bit_cutoff or to be dropped.
  PARALOGUE  both features long, both called, collections disjoint but
             similar. Genuine paralogues that the profiles cannot separate
             -- document it, do not "fix" it.

Distinguishing NESTED from PARALOGUE needs the collections, so pass
--collections to have it blastp each offending pair and print the identity.
"""
import argparse, collections, glob, json, os, subprocess, sys, tempfile


def span(f):
    l = (f.get("location") or [[None, 0, "+", 0]])[0]
    try:
        c, s, st, d = l[0], int(l[1]), l[2], int(l[3])
    except Exception:
        return None
    return (c, s - d + 1, s) if st == "-" else (c, s, s + d - 1)


def relate(a, b):
    """median identity and coverage of collection a against collection b"""
    t = tempfile.mkdtemp()
    subprocess.run(["makeblastdb", "-in", b, "-dbtype", "prot", "-out", t + "/db"],
                   capture_output=True)
    r = subprocess.run(["blastp", "-query", a, "-db", t + "/db",
                        "-outfmt", "6 qseqid pident qcovs bitscore",
                        "-max_target_seqs", "1", "-evalue", "10", "-num_threads", "2"],
                       capture_output=True, text=True)
    best = {}
    for line in r.stdout.strip().split("\n"):
        if not line:
            continue
        p = line.split("\t")
        best.setdefault(p[0], (float(p[1]), float(p[2]), float(p[3])))
    if not best:
        return None
    import statistics as st
    return (st.median([v[0] for v in best.values()]),
            st.median([v[1] for v in best.values()]),
            st.median([v[2] for v in best.values()]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eval", required=True)
    ap.add_argument("--json", required=True)
    ap.add_argument("--module", required=True)
    ap.add_argument("--collections", default=None,
                    help="collections/<Module> dir; enables the blastp column")
    ap.add_argument("--nested-frac", type=float, default=0.35)
    ap.add_argument("--expected", default="",
                    help="comma-separated CHILD:PARENT overlaps that are real "
                         "genome architecture, e.g. CP:POLY,MP:REP. These are "
                         "reported as EXPECTED rather than as findings.")
    a = ap.parse_args()

    expected = set()
    for pair in a.expected.split(","):
        if ":" in pair:
            c, par = pair.split(":", 1)
            expected.add((c.strip(), par.strip()))

    blk = json.load(open(a.json))[a.module]
    a2k = {v["anno"]: k for k, v in blk["features"].items()}

    exact = collections.Counter()
    nested = collections.Counter()
    who = collections.defaultdict(collections.Counter)
    total = collections.Counter()
    n = 0
    for q in sorted(glob.glob(os.path.join(a.eval, "*.qual.gto"))):
        try:
            d = json.load(open(q))
        except Exception:
            continue
        fs = [(a2k.get(f.get("function") or ""), span(f), (f.get("location") or [[None]])[0])
              for f in (d.get("features") or [])]
        fs = [x for x in fs if x[0] and x[1]]
        if not fs:
            continue
        n += 1
        byloc = collections.defaultdict(list)
        for k, s, raw in fs:
            total[k] += 1
            byloc[tuple(raw)].append(k)
        for _loc, ks in byloc.items():
            if len(ks) > 1:
                #  the same key twice on one locus is a different problem from
                #  two different keys: it means one feature was emitted twice,
                #  not that two profiles collided.
                u = sorted(set(ks))
                exact[tuple(u) if len(u) > 1 else (u[0] + " (x%d, same feature)" % len(ks),)] += 1
        for i, (k, s, _r) in enumerate(fs):
            for j, (k2, s2, _r2) in enumerate(fs):
                if i == j or k2 == k:
                    continue
                if s[0] == s2[0] and s[1] <= s2[2] and s2[1] <= s[2] \
                   and (s2[2] - s2[1]) > (s[2] - s[1]):
                    nested[k] += 1
                    who[k][k2] += 1
                    break

    print("  %s: %d scored genome(s) carrying features\n" % (a.module, n))
    if exact:
        print("  EXACT -- two features on identical coordinates:")
        for ks, c in exact.most_common():
            print("    %-30s %3d genome(s)  %.0f%%" % (" + ".join(ks), c, 100.0 * c / n))
    else:
        print("  EXACT -- none")

    hits = [(k, nested[k] / total[k]) for k in sorted(total)
            if total[k] and nested[k] / total[k] >= a.nested_frac]
    exp = [(k, f) for k, f in hits if (k, who[k].most_common(1)[0][0]) in expected]
    hits = [(k, f) for k, f in hits if (k, who[k].most_common(1)[0][0]) not in expected]
    if exp:
        print("\n  EXPECTED -- declared genome architecture, not a finding:")
        for k, frac in sorted(exp, key=lambda x: -x[1]):
            print("    %-11s %5.0f%% of %d calls inside %s"
                  % (k, 100 * frac, total[k], who[k].most_common(1)[0][0]))
    print("\n  NESTED -- called inside a longer feature at or above %.0f%%:" % (100 * a.nested_frac))
    if not hits:
        print("    none")
    for k, frac in sorted(hits, key=lambda x: -x[1]):
        par = who[k].most_common(1)[0][0]
        extra = ""
        if a.collections:
            fa = os.path.join(a.collections, k + ".fasta")
            fb = os.path.join(a.collections, par + ".fasta")
            if os.path.exists(fa) and os.path.exists(fb):
                r = relate(fa, fb)
                if r:
                    extra = "   blastp vs %s: id=%.0f%% cov=%.0f%% bits=%.0f" % (par, r[0], r[1], r[2])
                    if r[0] >= 90 and r[1] >= 90:
                        extra += "   -> SAME PROTEIN, merge"
                    elif r[1] >= 80:
                        extra += "   -> one family, merge"
                    else:
                        extra += "   -> promiscuous profile, raise bit_cutoff or drop"
        print("    %-11s %5.0f%% of %d calls inside %-11s%s"
              % (k, 100 * frac, total[k], par, extra))
    print("\n  A feature that fires inside a longer one is called on nearly every")
    print("  genome and so EARNS copy_num. Resolve these BEFORE apply_copy_num.py.")


main()
