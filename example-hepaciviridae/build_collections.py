#!/usr/bin/env python3
"""Bin the Hepacivirus dump into feature collections.

Module name is **Hepaciviridae**, not Flaviviridae. ICTV ratified the split of
the old family in 2025 (2025.006S.Amarillovirales_3reorgfam) into Flaviviridae,
Hepaciviridae and Pestiviridae. BV-BRC has not caught up and still files every
Hepacivirus record under `family = Flaviviridae`; the module name is written
into every output GTO as viral_family, so it follows ICTV, not the dump.

Layout: one polyprotein, C-E1-E2-p7-NS2-NS3-NS4A-NS4B-NS5A-NS5B.

The structure of this file follows example-orthoflavivirus: rules, then a
JUNCTION table of cleavage chemistry, then ambiguity / re-bin / length /
junction filters in that order.
"""
import collections
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
T = "Hepaci"
MODULE = "Hepaciviridae"

#  Strings that name no single protein. These are excluded BEFORE the feature
#  rules, because several of them contain a feature's own name as a substring
#  and would otherwise be binned as that feature.
NOT_MODELLED = [
    #  Precursors spanning two products. Binning one as its first component
    #  puts a sequence twice the right length into that collection and drags
    #  its length window with it.
    (r"NS3[-/ ]?4[Aa]", "precursor NS3-4A, spans two products"),
    (r"NS5AB", "precursor NS5A-NS5B, spans two products"),
    (r"^NS4$|^NS4 protein$|^non-?structural protein 4$",
     "NS4 is NS4A+NS4B, spans two products"),
    (r"E1[-/]E2|envelope protein E1/E2|core-E2",
     "precursor spanning two products"),
    #  'NS5' is genuinely ambiguous in this taxon: HCV has NS5A and NS5B and
    #  no NS5. The records carrying it split into a ~364 aa group and a
    #  ~1000 aa group -- a truncated something and the NS5A-NS5B precursor.
    #  Neither can be assigned from the string.
    (r"^(put\.? )?NS5( protein)?$|^non[- ]?structural( protein)? 5$",
     "NS5 is ambiguous: HCV has NS5A and NS5B, not NS5"),
    #  No protein named.
    (r"^hypothetical|^unnamed|^source:|^nonstructural gene$|^protein$",
     "names no protein"),
    (r"Fusion protein, Feo", "not a hepacivirus protein"),
]

#  key -> pattern. Order matters: the first match wins, so the longer and more
#  specific names come first.
#
#  The NS1 trap. In the early HCV literature E2 was called NS1, so `E2/NS1`,
#  `NS1/E2 protein`, `NS1 protein` and `boundary between NS1 and E2` all mean
#  E2 -- median 346-352 aa against E2's 363. This is a genus-dependent synonym
#  of exactly the kind references/annotation-triage.md warns about, and it is
#  worse than most because `NS1` names a real and completely different protein
#  in Orthoflavivirus, a module in the same repository. Bin on the string, but
#  the junction and length filters below are what actually keep it honest.
RULES = [
    ("POLY",  r"^(putative |elongated |precursor |flavivirus |[\w/.-]+ )?poly ?protein"
              r"( precursor| \(EC [\d.]+\))?$|^ORF1$"),
    ("NS5B",  r"NS ?5 ?B|RNA[- ]dependa?e?nt RNA (polymerase|ploymerase)|^RdRp$|"
              r"^polymerase( NS5B)?$|truncated RNA-dependent"),
    ("NS5A",  r"NS ?5 ?A|non[- ]?structural(( protein)? 5 ?A| 5 A)|NS5A phosphoprotein"),
    ("NS4B",  r"NS ?4 ?[Bb]|non[- ]?structural( protein)? 4 ?[Bb]"),
    ("NS4A",  r"NS ?4 ?[Aa]|non[- ]?structural( protein)? 4 ?[Aa]|NS4A peptide|"
              r"NS3/4A protease cofactor|NS3/4A proteinase cofactor"),
    ("NS3",   r"NS ?3|protease[- /]helicase|helicase|protease$|^protease|proteinase"),
    ("NS2",   r"NS ?2|non[- ]?structural( protein)? 2|non structural protein 2"),
    ("P7",    r"^p7|^P7|p7 protein|p7 peptide|protein p7|canal protein|protein p6|"
              r"processing product p7|p7_protein"),
    #  E2 first, because several E2 synonyms contain 'E1' (E1/E2 is already
    #  excluded above) and because the NS1 synonyms must not fall through to
    #  anything else.
    ("E2",    r"E ?2\b|E2_protein|e2 protein|gp70|envelop(e)? ?(protein )?2|"
              r"envelope glycoprotein E2|NS1|short E2"),
    ("E1",    r"E ?1\b|E1_protein|gp35|envelop(e)? ?(protein )?1|"
              r"envelope glycoprotein E1|^envelope"),
    ("C",     r"^core|core protein|core_protein|^C$|^C protein$|capsid|p22|"
              r"core nucleocapsid"),
    #  F / ARFP: a ribosomal frameshift product of the core coding region.
    #  Binned so it can be counted and examined at step 6b; whether it ships
    #  as a transcript_edit feature is decided there, not here.
    ("F",     r"protein F$|^F protein$|ARFP|alternate reading frame"),
]

#  Cleavage chemistry. HCV is cut by three enzymes with published specificities:
#
#    host signal peptidase   C/E1, E1/E2, E2/p7, p7/NS2   small residue at P1
#    NS2-3 autoprotease      NS2/NS3                      NS3 begins with A
#    NS3-4A serine protease  NS3/4A, 4A/4B, 4B/5A, 5A/5B  Cys or Thr at P1,
#                                                         Ser or Ala at P1'
#
#  Note this is the OPPOSITE expectation to Orthoflavivirus, where NS2B-NS3
#  cuts after a dibasic. Carrying over the R/K intuition would reject the whole
#  taxon and look like a finding.
SMALL = set("GSATC")     # signalase P1
CYS   = set("CT")        # NS3-4A P1
SERALA= set("SA")        # NS3-4A P1'
MET   = set("M")
JUNCTION = {
    "POLY": (MET,    None),      # initiator Met / polyprotein C-terminus
    "C":    (MET,    SMALL),     # initiator Met / C-E1 signalase
    "E1":   (None,   SMALL),     # signalase both sides
    "E2":   (None,   SMALL),
    "P7":   (None,   SMALL),
    "NS2":  (None,   None),      # signalase P1' / NS2-3 autoprotease, no motif
    "NS3":  (set("A"), CYS),     # NS2-3 autoprotease gives A / NS3-4A gives C
    "NS4A": (SERALA, CYS),
    "NS4B": (SERALA, CYS),
    "NS5A": (SERALA, CYS),
    "NS5B": (SERALA, None),      # NS3-4A P1' / polyprotein C-terminus
    "F":    (MET,    None),      # shares core's initiator Met
}

STD = set("ACDEFGHIKLMNPQRSTVWY")
RUNT, GIANT = 0.70, 1.30

#  A floor applied BEFORE the median is computed, for features whose record
#  set is mostly partial. POLY needs one and needs it badly: 73.4% of the
#  12,263 records matching a polyprotein string are under 2,000 aa against a
#  true HCV polyprotein of ~3,011, so the median lands at 1,400 -- inside a
#  mode made of half-polyproteins -- and the 0.7-1.3x window it implies
#  (979-1,820) keeps the halves and excludes every real polyprotein.
#
#  The distribution is tri-modal with a clean empty band: 4,900 fragments at
#  250-750, 3,570 halves at 1,250-1,750, NOTHING between 2,000 and 2,750, then
#  3,212 full length at 2,750-3,250. Anywhere in that gap separates them, so
#  the cut is 2,400 rather than a tuned number.
PRE_FLOOR = {"POLY": 2400}

#  Ambiguity policy for this taxon: keep a sequence carrying at most AMB_MAX
#  ambiguous residues, drop the rest. Curator's call, 22 September 2026.
#
#  Orthoflavivirus uses zero tolerance and that is right there -- its worst
#  offenders are genuinely junk, Orthoflavivirus tembusu VN0613 carrying 1,261
#  X in one polyprotein, 37% of the sequence. Hepacivirus is different in a way
#  that is a property of how the data is produced, not of the module: HCV is
#  sequenced from clinical samples, each of which is a quasispecies, so Sanger
#  reads give overlapping peaks at polymorphic sites and the submitter records
#  an IUPAC code that translates to X. Measured on NS3: 97.4% of ambiguous
#  positions are ISOLATED single residues and 94% of affected sequences have no
#  run at all, which is the signature of population polymorphism rather than an
#  assembly gap. Only 2% carry a run of >=5, which is what a real gap looks like.
#
#  The cost of zero tolerance here is not in sequence count but in diversity:
#  it discards 286 of NS3's 384 distinct 90%-identity clusters, 74.5%, because
#  whole HCV genotypes exist in BV-BRC only as population sequences.
AMB_MAX = 2

#  A premature stop is not ambiguity. It means the reading frame or the
#  annotation is wrong and everything downstream of it is suspect, so it
#  disqualifies regardless of count -- 1,108 of them in NS3 alone.
FATAL = set("*")

X_SOLE_MAX, X_SOLE_ABS = 0.01, 2


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


def main():
    rules = compiled(RULES)
    notmod = [(p, w) for p, w in NOT_MODELLED]
    md5s = [l.strip() for l in open(os.path.join(HERE, "%s.uniq.md5" % T))]
    ann = [l.rstrip("\n").split("\t") for l in open(os.path.join(HERE, "%s.uniq.id_ann" % T))]
    seq = [l.rstrip("\n").split("\t") for l in open(os.path.join(HERE, "%s.uniq.seq" % T))]
    n = len(md5s)
    if not (len(ann) == len(seq) == n):
        sys.exit("uniq.* files are not the same length -- run check_dump.py")

    per_md5 = collections.Counter()
    for l in open(os.path.join(HERE, "%s.id_md5" % T), errors="replace"):
        p = l.rstrip("\n").split("\t")
        if len(p) >= 2 and p[1].strip():
            per_md5[p[1].strip()] += 1

    out = collections.defaultdict(list)
    unassigned = collections.Counter()
    notmodelled = collections.Counter()
    feat_counts = collections.Counter()
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

    #  ---- ambiguity: any non-standard residue, sole-representative exception
    amb = lambda q: sum(1 for c in q if c not in STD)
    xdropped = collections.Counter(); xfatal = collections.Counter(); xkept = []
    for k in list(out):
        rows = out[k]
        fatal = [r for r in rows if set(r[2]) & FATAL]
        xfatal[k] = len(fatal)
        rows = [r for r in rows if not (set(r[2]) & FATAL)]
        keep = [r for r in rows if amb(r[2]) <= AMB_MAX]
        #  the compelling-reason clause: a taxon whose every sequence for this
        #  feature exceeds the limit keeps its cleanest one rather than vanishing
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

    #  ---- length: floor and ceiling, both relative to the feature's median
    dropped = collections.Counter(); overlong = collections.Counter()
    prefloored = collections.Counter()
    for k in list(out):
        if k in PRE_FLOOR:
            before = len(out[k])
            out[k] = [r for r in out[k] if len(r[2]) >= PRE_FLOOR[k]]
            prefloored[k] = before - len(out[k])
        lens = sorted(len(s) for _, _, s in out[k])
        if not lens:
            continue
        med = lens[len(lens) // 2]
        keep = [r for r in out[k] if RUNT * med <= len(r[2]) <= GIANT * med]
        dropped[k] = sum(1 for r in out[k] if len(r[2]) < RUNT * med)
        overlong[k] = sum(1 for r in out[k] if len(r[2]) > GIANT * med)
        out[k] = keep

    #  ---- cleavage chemistry
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

    cdir = os.path.join(HERE, "collections", MODULE)
    os.makedirs(cdir, exist_ok=True)
    for k, rows in sorted(out.items()):
        with open(os.path.join(cdir, k + ".fasta"), "w") as fh:
            for fid, _a, s in rows:
                fh.write(">%s\n%s\n" % (fid, s))

    td = os.path.join(HERE, "collections")
    def dump(name, header, rowsrc):
        with open(os.path.join(td, name), "w") as fh:
            fh.write(header + "\n")
            for r in rowsrc:
                fh.write("\t".join(str(x) for x in r) + "\n")
    dump("PREFLOOR_DROPPED.tsv", "feature\tdropped\tfloor",
         ((k, prefloored[k], PRE_FLOOR[k]) for k in sorted(prefloored)))
    dump("RUNTS_DROPPED.tsv", "feature\tdropped\tthreshold",
         ((k, dropped[k], "%d%% of median" % (RUNT * 100)) for k in sorted(dropped)))
    dump("OVERLONG_DROPPED.tsv", "feature\tdropped\tthreshold",
         ((k, overlong[k], "%d%% of median" % (GIANT * 100)) for k in sorted(overlong)))
    dump("AMBIGUOUS_DROPPED.tsv", "feature\tdropped\tof_which_premature_stop",
         ((k, xdropped[k], xfatal[k]) for k in sorted(xdropped)))
    dump("AMBIGUOUS_KEPT.tsv", "feature\ttaxon\tfeature_id\tn_amb\tlength", sorted(xkept))
    dump("JUNCTION_DROPPED.tsv", "feature\tfeature_id\tlength\tfirst\tlast\tbad_end", sorted(jdetail))
    dump("JUNCTION_SUMMARY.tsv", "feature\tdropped", ((k, jdropped[k]) for k in sorted(jdropped)))
    dump("UNASSIGNED_TRACKING.tsv", "features\tannotation",
         ((c, a) for a, c in unassigned.most_common()))
    dump("NOT_MODELLED.tsv", "features\tannotation\treason",
         ((c, a, w) for (a, w), c in notmodelled.most_common()))

    print("  %-7s %8s %9s %8s %8s %8s %8s %9s" %
          ("feature", "uniq", "features", "prefloor", "runt", "overlong", "amb", "junction"))
    tot = 0
    for k in sorted(out):
        print("  %-7s %8d %9d %8d %8d %8d %8d %9d" %
              (k, len(out[k]), feat_counts[k], prefloored[k], dropped[k], overlong[k],
               xdropped[k], jdropped[k]))
        tot += feat_counts[k]
    allf = tot + sum(notmodelled.values()) + sum(unassigned.values())
    print("  %-7s %8s %9d" % ("NOTMOD", "", sum(notmodelled.values())))
    print("  %-7s %8d %9d" % ("UNASSIG", len(unassigned), sum(unassigned.values())))
    print("\n  binned %d of %d features = %.1f%%" % (tot, allf, 100.0 * tot / allf))


main()
