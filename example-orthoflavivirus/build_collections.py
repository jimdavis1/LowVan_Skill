#!/usr/bin/env python3
"""Bin Orthoflavivirus protein annotations into feature collections.

The genome is one long ORF cleaved into ten mature products, so the module
follows the Togaviridae pattern: a polyprotein CDS plus mat_peptides for the
products, rather than pretending each product is its own gene.

    C - prM/M - E - NS1 - NS2A - NS2B - NS3 - NS4A - 2K - NS4B - NS5

Two traps in this taxon's vocabulary, both of which cost real sequences if
the rules are written carelessly:

**prM, pr and M are three different things.** prM is the precursor; furin
cleaves it into pr (released) and M (retained in the virion). BV-BRC uses
"membrane glycoprotein precursor" and "membrane glycoprotein precursor prM"
for the precursor, "protein pr" for the released peptide, and "membrane
glycoprotein" and "matrix protein M" for the mature chain. A rule matching
plain /membrane glycoprotein/ first would swallow all three.

**"anchored capsid" is not "capsid".** The anchored form retains the
C-terminal signal that anchors it to the ER membrane before NS2B-3 cleaves
it; the mature capsid does not. They differ by ~20 residues, so binning them
together produces a cluster with a jagged C-terminus.

Ordering therefore matters: every rule is anchored and the most specific
pattern is tested first. RULES is an ordered list, not a dict.
"""
import collections
import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
T = "Ortho"
MODULE = "Orthoflavivirus"

#  Ordered. First match wins, so specific before general.
RULES = [
    # --- the whole ORF -------------------------------------------------
    ("POLY",  r"^(putative )?(genome )?polyprotein"),
    ("POLY",  r"^polyprotein\b"),

    # --- capsid: anchored form BEFORE mature ---------------------------
    ("ANCHC", r"^anchored (capsid|core) protein( C)?$"),
    ("ANCHC", r"^anchored capsid protein ancC$"),
    ("C",     r"^(capsid|core) protein( C)?$"),
    ("C",     r"^capsid$"),

    # --- prM / pr / M : precursor, released peptide, mature chain ------
    ("PRM",   r"^(membrane glycoprotein )?precursor( prM)?$"),
    ("PRM",   r"^membrane glycoprotein precursor( prM)?$"),
    ("PRM",   r"^prM\b"),
    ("PRM",   r"^protein prM$"),
    ("PR",    r"^(protein )?pr$"),
    ("PR",    r"^peptide pr$"),
    ("M",     r"^(matrix |membrane )?(glyco)?protein M$"),
    ("M",     r"^membrane glycoprotein( M)?$"),
    ("M",     r"^matrix protein( M)?$"),

    # --- envelope ------------------------------------------------------
    ("E",     r"^envelope (glyco)?protein( E)?$"),
    ("E",     r"^protein E$"),
    ("E",     r"^E protein$"),

    # --- non-structural. NS2A/NS2B and NS4A/NS4B before bare NS2/NS4 ---
    ("NS1",   r"^(non-?structural )?protein ?NS1$"),
    ("NS1",   r"^nonstructural protein 1$"),
    ("NS2A",  r"^(non-?structural )?protein ?NS2[aA]$"),
    ("NS2A",  r"^nonstructural protein 2A$"),
    ("NS2B",  r"^(non-?structural )?protein ?NS2[bB]$"),
    ("NS2B",  r"^nonstructural protein 2B$"),
    ("NS3",   r"^(non-?structural )?protein ?NS3$"),
    ("NS3",   r"^nonstructural protein 3$"),
    ("NS3",   r"^(NS3 )?(protease|helicase)(/helicase)?( NS3)?$"),
    ("NS4A",  r"^(non-?structural )?protein ?NS4[aA]$"),
    ("NS4A",  r"^nonstructural protein 4A$"),
    ("NS4B",  r"^(non-?structural )?protein ?NS4[bB]$"),
    ("NS4B",  r"^nonstructural protein 4B$"),
    ("NS5",   r"^(non-?structural )?protein ?NS5$"),
    ("NS5",   r"^nonstructural protein 5$"),
    ("NS5",   r"^RNA-dependent RNA polymerase( NS5)?$"),
    ("NS5",   r"^(methyltransferase|RdRp)( NS5)?$"),

    # --- forms the first pass missed ----------------------------------
    #  "membrane glycoprotein precursor M" and "PreM" are prM, not M:
    #  median length 166-168 aa against 75 for the mature chain. Binning
    #  them by wording rather than by length would have put the precursor
    #  in the M collection.
    ("PRM",   r"^membrane glycoprotein precursor M$"),
    ("PRM",   r"^PreM( protein)?$"),
    ("PRM",   r"^pre-?membrane( protein)?( prM)?$"),
    ("M",     r"^small envelope protein M$"),
    #  bare and suffixed non-structural names
    ("NS1",   r"^NS1( protein)?$"),
    ("NS2A",  r"^NS2[aA]( protein)?$"),
    ("NS2B",  r"^NS2[bB]( protein)?$"),
    ("NS3",   r"^NS3( protein)?$"),
    ("NS4A",  r"^NS4[aA]( protein)?$"),
    ("NS4B",  r"^NS4[bB]( protein)?$"),
    ("NS5",   r"^NS5( protein)?$"),
    ("E",     r"^E$"),
    ("C",     r"^C( protein)?$"),

    # --- NS1', the -1 frameshift product ------------------------------
    #  Real, and not collinear with the polyprotein: JEV-serogroup viruses
    #  (JEV, West Nile, Usutu) produce NS1' by a -1 ribosomal frameshift
    #  near the start of NS2A, giving NS1 plus a ~52-residue extension. It
    #  gets its own collection rather than being folded into NS1, because a
    #  profile trained on both would have a ragged C-terminus, and because
    #  step 6b of the workflow has to decide whether it ships as a
    #  transcript_edit feature.
    ("NS1P",  r"^(non-?structural protein )?NS1'?( prime)?( protein)?$".replace("NS1'?( prime)?", "NS1('| prime)")),
    ("NS1P",  r"^NS1'( protein)?$"),
    ("NS1P",  r"^nonstructural protein NS1 prime$"),

    # --- remaining wording variants, resolved against the same products
    ("PRM",   r"^M protein precursor$"),
    ("PRM",   r"^membrane precursor \(prM\) protein$"),
    ("PRM",   r"^prM protein$"),
    ("E",     r"^envelope( \(E\))?( protein)?$"),
    ("NS5",   r"^RNA-directed RNA polymerase( NS5)?$"),
    ("NS3",   r"^serine protease( NS3)?$"),
    ("NS3",   r"^NS3 (serine )?protease$"),
    ("C",     r"^virion capsid( \(virC\))?( protein)?$"),
    ("ANCHC", r"^anchored capsid( \(anchC\))?( protein)?$"),
    ("ANCHC", r"^anc[CD]$"),
    ("M",     r"^membrane( \(M\))?( protein)?$"),

    # --- the 2K signal peptide between NS4A and NS4B -------------------
    ("2K",    r"^(protein )?2K( peptide| protein)?$"),
    ("2K",    r"^peptide 2K$"),
]

#  Strings that are real but deliberately not modelled, so they do not sit
#  in UNASSIGNED looking like an oversight.
NOT_MODELLED = [
    (r"fragment", "explicit fragment"),
    (r",\s*[NC]-termin", "positional half, not a whole product"),
    (r"^(truncated|partial)\b", "explicit partial"),
    (r"^hypothetical protein$", "no functional claim"),
    (r"^unnamed protein product$", "no functional claim"),
    (r"^ORF\d*$", "positional label, not a homology group"),
]


def compiled(rules):
    return [(k, re.compile(p, re.I)) for k, p in rules]


def classify(anno, rules, notmod):
    a = " ".join(anno.split())
    #  Partials are excluded BEFORE the feature rules, not after. A rule
    #  anchored with ^ still matches a prefix, so "putative genome
    #  polyprotein, C-terminal" satisfies the POLY pattern and lands in the
    #  collection as though it were a whole ORF. The full Orthoflavivirus
    #  polyprotein is ~3392 aa; the C-terminal and N-terminal halves are
    #  ~1736 and ~1776, and "fragment" is ~327. Training a profile on that
    #  mixture gives a jagged model of a gene that is actually uniform.
    for rx, why in notmod:
        if re.search(rx, a, re.I):
            return None, why
    for k, rx in rules:
        if rx.match(a):
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

    #  count features, not unique sequences: a protein shared by 400 dengue
    #  genomes is one uniq row and 400 features.
    per_md5 = collections.Counter()
    for l in open(os.path.join(HERE, "%s.id_md5" % T)):
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
            fid = ann[i][0] if ann[i] else md5s[i]
            out[key].append((fid, a, s))
            feat_counts[key] += nfeat
        elif why:
            notmodelled[(a, why)] += nfeat
        else:
            unassigned[a] += nfeat

    #  Re-binning runs FIRST, before any length or composition filter.
    #  A record's length is only meaningful relative to the feature it
    #  actually belongs to: a 3,391-residue polyprotein filed as NS5 is
    #  4x overlong for NS5 and exactly right for POLY, and a 139-residue
    #  anchored capsid filed as C is overlong for C (median 100) and normal
    #  for ancC (median 114). Filtering before re-binning silently discarded
    #  36 of the 108 ancC records and both polyprotein records.
    #  Re-bin records that homology identifies as a different feature. This
    #  is a second pass by construction: REASSIGNED.tsv is produced by
    #  reassign.py from the JUNCTION_DROPPED.tsv that the first pass wrote,
    #  so the order is build -> identify -> rebuild. Absent the file this is
    #  a no-op and the records simply stay dropped.
    #
    #  Only records whose sequence covers a clean member of the target
    #  feature at >=80% both ways, at >=50% identity, are listed there. 108
    #  of the 113 are anchored capsid filed as mature C: ancC and C overlap
    #  in length (C p90 114, ancC p10 110) so no length rule separates them,
    #  and the annotation string says "capsid protein" for both. Homology
    #  plus the junction chemistry does separate them, because their
    #  C-termini are cut by different enzymes.
    #
    #  The other 1,309 junction failures are NOT re-binned: 1,094 match
    #  their own feature best -- the right protein with broken ends, which
    #  has no bin to move to and would need trimming -- 188 match nothing at
    #  >=80% coverage, and 27 fail the target feature's own junction rule.
    rmap = {}
    rfile = os.path.join(HERE, "REASSIGNED.tsv")
    if os.path.exists(rfile):
        for line in list(open(rfile))[1:]:
            f = line.rstrip("\n").split("\t")
            rmap[(f[0], f[2])] = f[1]
    reassigned = collections.Counter()
    if rmap:
        for k in list(out):
            stay = []
            for r in out[k]:
                dest = rmap.get((k, r[0]))
                if dest and dest != k:
                    out.setdefault(dest, []).append(r)
                    reassigned[(k, dest)] += 1
                else:
                    stay.append(r)
            out[k] = stay

    #  Drop runts. The partial filter catches sequences the source *labels*
    #  as fragments; it cannot catch one labelled "polyprotein" that is 4
    #  residues long. Every feature here is a cleavage product of one
    #  polyprotein and so is tightly sized -- E is 495 aa at both the median
    #  and the 5th percentile -- which makes a fraction-of-median cut safe
    #  and a profile trained on the runts badly ragged. At 70% this removes
    #  233 of 18,374 POLY sequences (1.3%) and under 1% of most features.
    #  ...and a matching bound above. The pipeline had only a floor, which is
    #  why positional halves leaking into POLY needed a hand-written fix and
    #  why M -- a 75-residue protein -- still held sequences up to 500 aa
    #  after every other filter. At 1.3x median this removes 50 sequences
    #  across all fifteen features (M 26, NS4A 13, POLY 10, NS2A 1, NS2B 1)
    #  and there is a clean gap either side of the threshold: the longest
    #  legitimate variant is Culex flavivirus ancC at 1.22x, and the junk
    #  starts at 3.6x (a POLY of 12,309 residues, roughly four concatenated
    #  polyproteins). 1.2x would discard the Culex ancC; 1.5x removes the
    #  same 50 as 1.3x, so nothing real lives between them.
    RUNT = 0.70
    GIANT = 1.30
    dropped = collections.Counter()
    overlong = collections.Counter()
    for k in list(out):
        lens = sorted(len(s) for _, _, s in out[k])
        med = lens[len(lens) // 2]
        keep = [r for r in out[k]
                if RUNT * med <= len(r[2]) <= GIANT * med]
        dropped[k] = sum(1 for r in out[k] if len(r[2]) < RUNT * med)
        overlong[k] = sum(1 for r in out[k] if len(r[2]) > GIANT * med)
        out[k] = keep

    #  Drop ambiguous sequence. A PSSM column built on X carries no
    #  information, and 17.2% of this taxon's sequences (11,705 of 68,202)
    #  contain at least one. The worst are not marginal: Orthoflavivirus
    #  tembusu VN0613 contributes a POLY of 1,261 X, an E of 392 and an NS4A
    #  that is 126 X out of ~130 residues. Curator's instruction is to blow
    #  away any sequence with an X unless there is a compelling reason to
    #  keep it.
    #
    #  The one compelling reason is coverage: a sequence that is the only one
    #  its taxon has for that feature, spoiled by an incidental residue or
    #  two, is the sole training data for a distinct species. Dropping it
    #  buys nothing and costs a species. Sokoluk virus NS4A, Dakar bat virus
    #  M and NS3, Phnom Penh bat virus NS3/E/NS5/POLY and Flavivirus sp. M138
    #  POLY (4 X in ~3,400) are kept on that basis; tembusu VN0613 is not,
    #  because its X fraction is 7-37% and the species is represented by
    #  other strains and a rep contig.
    #  Two thresholds, whichever is kinder, because neither alone is right
    #  at both ends of the length range. A fraction is the correct measure
    #  for a 3,400-residue polyprotein and wrong for a 75-residue M, where a
    #  single incidental X is already 1.33%; an absolute count is correct for
    #  M and far too loose for POLY. Dakar bat virus M -- one X, sole
    #  representative -- is the case that requires the absolute arm.
    #  X is not the only ambiguity code the source emits. B is Asx (D or N),
    #  Z is Glx (E or Q), J is Leu-or-Ile, and one record carries a literal
    #  '*' -- a stop codon in the middle of a protein. They are the same
    #  class of thing as X and are filtered by the same rule: 26 B, 2 Z, 1 J
    #  and 1 '*' survived an X-only filter.
    STD = set("ACDEFGHIKLMNPQRSTVWY")
    X_SOLE_MAX = 0.01   # sole representative kept if this clean by fraction
    X_SOLE_ABS = 2      # ...or this clean by count
    xdropped = collections.Counter()
    xkept = []
    for k in list(out):
        bytax = collections.defaultdict(list)
        for r in out[k]:
            bytax[r[0].split("|")[1].split(".")[0]].append(r)
        keep = []
        for t, rows in bytax.items():
            clean = [r for r in rows if not (set(r[2]) - STD)]
            if clean:
                keep.extend(clean)
                continue
            #  every sequence this taxon has for this feature carries an X
            amb = lambda q: sum(1 for c in q if c not in STD)
            best = min(rows, key=lambda r: amb(r[2]) / len(r[2]))
            nx = amb(best[2])
            if nx <= X_SOLE_ABS or nx / len(best[2]) <= X_SOLE_MAX:
                keep.append(best)
                xkept.append((k, t, best[0], nx, len(best[2])))
        xdropped[k] = len(out[k]) - len(keep)
        out[k] = keep

    #  Drop features whose termini contradict the published cleavage
    #  chemistry. Every feature here is a cleavage product of one
    #  polyprotein, so its two ends are not a matter of opinion: they are
    #  where a named enzyme cuts.
    #
    #    Chambers 1990 (PMID 2174669)  the polyprotein cleavage map
    #    Lin 1993     (PMID 8445732)   the NS4A/2K site is cut by NS2B-NS3
    #                                  and sits exactly 23 residues upstream
    #                                  of the 2K/NS4B signalase site, which
    #                                  is why 2K is a 23-mer
    #    Wahaab 2021  (PMID 33494395)  NS2B-NS3 requires two basic residues
    #                                  (K-R, R-R, R-K, occasionally Q-R) at
    #                                  P2-P1 and a small residue (G, S, A)
    #                                  at P1'
    #
    #  So an NS2B-NS3-generated C-terminus ends in R or K, a signalase-
    #  generated one ends in a small residue, and an NS2B-NS3-generated
    #  N-terminus begins with a small residue. A feature that breaks its own
    #  rule is mis-annotated or mis-placed, and is thrown out rather than
    #  trimmed -- curator's instruction, and the right call: a trimmed
    #  boundary is a guess, while the 56k sequences that pass are evidence.
    #
    #  This is what keeps the alignment ends stable. It catches what neither
    #  the wording rules nor the length rules can:
    #    * 18-residue chunks of the NS3 helicase filed as 2K, because the
    #      annotation string said 2K and 18 aa sits inside 2K's window
    #    * anchored capsid filed as mature C -- 187 sequences whose real
    #      C-terminus is the ancC/prM signalase site. C and ancC overlap in
    #      length (C p90 114, ancC p10 110) so no length rule separates
    #      them, but their C-termini are cut by different enzymes and so
    #      differ by residue class
    SMALL = set("GSATC")   # signalase P1, and NS2B-NS3 P1'
    BASIC = set("RK")      # NS2B-NS3 and furin P1
    MET   = set("M")       # polyprotein initiator
    JUNCTION = {
        "ANCHC": (MET,   SMALL), "C":    (MET,   BASIC),
        "PRM":   (None,  SMALL), "PR":   (None,  BASIC),
        "M":     (None,  SMALL), "E":    (None,  SMALL),
        "NS1":   (None,  SMALL), "NS2A": (None,  BASIC),
        "NS2B":  (SMALL, BASIC), "NS3":  (SMALL, BASIC),
        "NS4A":  (SMALL, BASIC), "2K":   (SMALL, SMALL),
        "NS4B":  (None,  BASIC), "NS5":  (SMALL, None),
        "POLY":  (MET,   None),
    }
    jdropped = collections.Counter()
    jdetail = []
    for k in list(out):
        if k not in JUNCTION:
            continue
        fs, ls = JUNCTION[k]
        keep = []
        for r in out[k]:
            seq = r[2]
            badn = fs is not None and seq[0] not in fs
            badc = ls is not None and seq[-1] not in ls
            if badn or badc:
                jdetail.append((k, r[0], len(seq), seq[0], seq[-1],
                                "N" if badn and not badc else
                                "C" if badc and not badn else "NC"))
            else:
                keep.append(r)
        jdropped[k] = len(out[k]) - len(keep)
        out[k] = keep

    cdir = os.path.join(HERE, "collections", MODULE)
    os.makedirs(cdir, exist_ok=True)
    for k, rows in sorted(out.items()):
        with open(os.path.join(cdir, k + ".fasta"), "w") as fh:
            for fid, a, s in rows:
                fh.write(">%s\n%s\n" % (fid, s))
    with open(os.path.join(HERE, "collections", "RUNTS_DROPPED.tsv"), "w") as fh:
        fh.write("feature\tdropped\tthreshold\n")
        for k in sorted(dropped):
            fh.write("%s\t%d\t%.0f%% of median\n" % (k, dropped[k], RUNT * 100))
    with open(os.path.join(HERE, "collections", "OVERLONG_DROPPED.tsv"), "w") as fh:
        fh.write("feature\tdropped\tthreshold\n")
        for k in sorted(overlong):
            fh.write("%s\t%d\t%.0f%% of median\n" % (k, overlong[k], GIANT * 100))
    with open(os.path.join(HERE, "collections", "AMBIGUOUS_DROPPED.tsv"), "w") as fh:
        fh.write("feature\tdropped_for_X\n")
        for k in sorted(xdropped):
            fh.write("%s\t%d\n" % (k, xdropped[k]))
    with open(os.path.join(HERE, "collections", "REASSIGNMENTS_APPLIED.tsv"), "w") as fh:
        fh.write("from_feature\tto_feature\tmoved\n")
        for (a, b), c in sorted(reassigned.items()):
            fh.write("%s\t%s\t%d\n" % (a, b, c))
    with open(os.path.join(HERE, "collections", "JUNCTION_DROPPED.tsv"), "w") as fh:
        fh.write("feature\tfeature_id\tlength\tfirst\tlast\tbad_end\n")
        for row in sorted(jdetail):
            fh.write("%s\t%s\t%d\t%s\t%s\t%s\n" % row)
    with open(os.path.join(HERE, "collections", "JUNCTION_SUMMARY.tsv"), "w") as fh:
        fh.write("feature\tdropped_for_junction\n")
        for k in sorted(jdropped):
            fh.write("%s\t%d\n" % (k, jdropped[k]))
    with open(os.path.join(HERE, "collections", "AMBIGUOUS_KEPT.tsv"), "w") as fh:
        fh.write("feature\ttaxon\tfeature_id\tn_X\tlength\treason\n")
        for k, t, fid, nx, ln in sorted(xkept):
            fh.write("%s\t%s\t%s\t%d\t%d\tsole representative of taxon "
                     "for this feature, %.2f%% X\n" % (k, t, fid, nx, ln,
                                                       100.0 * nx / ln))

    tdir = os.path.join(HERE, "collections")
    with open(os.path.join(tdir, "UNASSIGNED_TRACKING.tsv"), "w") as fh:
        fh.write("features\tannotation\n")
        for a, c in unassigned.most_common():
            fh.write("%d\t%s\n" % (c, a))
    with open(os.path.join(tdir, "NOT_MODELLED.tsv"), "w") as fh:
        fh.write("features\tannotation\treason\n")
        for (a, why), c in notmodelled.most_common():
            fh.write("%d\t%s\t%s\n" % (c, a, why))

    tot = sum(feat_counts.values()) + sum(unassigned.values()) + sum(notmodelled.values())
    print("  %-8s %8s %8s" % ("feature", "uniq", "features"))
    for k in sorted(out):
        print("  %-8s %8d %8d" % (k, len(out[k]), feat_counts[k]))
    print("  %-8s %8s %8d" % ("NOTMOD", "", sum(notmodelled.values())))
    print("  %-8s %8s %8d" % ("UNASSIG", len(unassigned), sum(unassigned.values())))
    print("\n  binned %d of %d features = %.1f%%"
          % (sum(feat_counts.values()), tot, 100.0 * sum(feat_counts.values()) / max(tot, 1)))


if __name__ == "__main__":
    main()
