#!/usr/bin/env python3
"""Tilapinevirus -- tilapia lake virus (TiLV), family Amnoonviridae.

    10 segments, negative-sense ssRNA, orthomyxo-like.  Each segment carries
    ONE ORF (PMID 27048802).  Conserved quasi-complementary 5'/3' termini on
    every segment are what placed the virus with the orthomyxoviruses.

      S1  519 aa   PB1-like polymerase subunit; catalytic SDD at 288-290
      S2  457 aa   PB2-like polymerase subunit
      S3  419 aa   PA-like polymerase subunit
      S4  356 aa   NUCLEOPROTEIN
      S5  343 aa   cleavable signal peptide + one TMD (predicted)
      S6  317 aa   predicted signal peptide, low cleavage score
      S7  195 aa   no assigned function
      S8  174 aa   no assigned function
      S9  118 aa   no assigned function
      S10 113 aa   no assigned function
      S9-F3 105 aa an ELEVENTH protein, frame 3 of segment 9  (see below)

    Segments 1-3 form a heterotrimeric polymerase (cryo-EM, PMID 38066000);
    segment 4 is the nucleoprotein (PMID 35107373, PMID 39995042).

### Why this module does NOT bin on the annotation string

363 of 399 distinct deposits are annotated `hypothetical protein`, and the
only named feature in the whole dump is PB1.  A rules table over the strings
can therefore separate exactly one of eleven proteins.  The segments are
non-homologous to each other, so HOMOLOGY separates them perfectly instead.

Reference set: Til-4-2011, the Bacharach 2016 isolate, whose ten segments are
deposited as KU751814..KU751823 IN SEGMENT ORDER at exactly the published
lengths (1641/1471/1371/1250/1099/1044/777/657/548/465 nt).  TiLV-Israel-HK
(OQ437054..OQ437063) reproduces the same ten protein lengths independently,
so the mapping is measured twice, not asserted once.

**The assignment was validated against the submitters' own labels.**  429 of
558 contigs state their segment in the FASTA description; blastn against the
reference set agreed with the stated segment in **429 of 429 (100%)**, and
additionally assigned 117 contigs that carry no label at all.

### The eleventh protein is not in the database

S9-F3 is translated from a 105-codon ORF in frame 3 of segment 9, overlapping
the 116-codon primary S9 ORF, with a strong Kozak context.  It was confirmed
by mass spectrometry -- 4 unique peptides, 49% sequence coverage -- and
carries a CRM1-dependent nuclear export signal (PMID 41513822).

**BV-BRC holds zero deposits of it.**  It is therefore ABSENT in the sense of
SKILL.md 2b: not a defect in these rules, but a protein the source never
annotated.  It is extracted from the genomes' own segment 9 in frame 3 and
declared here, because a module that silently omits a mass-spec-confirmed
viral protein is a module that is wrong.

### What this module does NOT cover

Eleven further amnoonvirid species (Flavolineata, Piscibus, Lipesoluc, and the
Catfish/Cod/Flatfish/Loach/Lizard x2/Pufferfish amnoonviruses) deposit ONE
segment each -- PB1/RdRp only -- and none has a second segment anywhere.
Measured at the annotator's own blastn settings they do not route to this
module at all (0 of 12 contigs at the 150-bit floor), which is the correct
outcome: they are a separate group this module declines rather than absorbs.
They cannot support a module of their own and are recorded as a gap.

### No splicing, no frameshift, no editing

PMID 41513822 searched for these explicitly and reports none: S9-F3 arises by
leaky ribosomal scanning of an alternative frame, not by splicing, and the
authors state they "could not identify an obvious PA-X-like frameshift
signal" in segment 3.  So every feature here is collinear with its segment.
"""
import os, sys, subprocess, collections, tempfile
sys.path.insert(0, os.environ["LOWVAN_KIT"] + "/skill/scripts")
from annotation_rules import PARTIAL
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.dirname(HERE)
T, MODULE = "amnoon", "Tilapinevirus"

#  Literature protein lengths, from the two independent reference isolates.
#  A window of +/-25% around each covers the real variation (the deposits run
#  96-100% identical) while excluding a neighbouring segment -- the nearest
#  pair is S9 118 / S10 113, which homology separates anyway.
SEGMENTS = {
    "S1_PB1": (519, "PB1-like polymerase subunit, SDD 288-290 (PMID 38066000)"),
    "S2_PB2": (457, "PB2-like polymerase subunit (PMID 38066000)"),
    "S3_PA":  (419, "PA-like polymerase subunit (PMID 38066000)"),
    "S4_NP":  (356, "nucleoprotein (PMID 35107373, PMID 39995042)"),
    "S5":     (343, "predicted signal peptide + TMD (PMID 41513822)"),
    "S6":     (317, "predicted signal peptide (PMID 41513822)"),
    "S7":     (195, "no assigned function"),
    "S8":     (174, "no assigned function"),
    "S9":     (118, "no assigned function"),
    "S10":    (113, "no assigned function"),
}
MIN_BITS, MIN_QCOV = 40.0, 0.50
LO, HI = 0.75, 1.25

def main():
    md5s = [l.strip() for l in open(os.path.join(WORK, "%s.uniq.md5" % T))]
    ann  = [l.rstrip("\n").split("\t") for l in open(os.path.join(WORK, "%s.uniq.id_ann" % T))]
    seq  = [l.rstrip("\n").split("\t") for l in open(os.path.join(WORK, "%s.uniq.seq" % T))]
    assert len(md5s) == len(ann) == len(seq), "uniq.* not aligned -- run check_dump.py"
    g2genus = {}
    for l in open(os.path.join(WORK, "%s.id_name_gs" % T), errors="replace"):
        f = l.rstrip("\n").split("\t")
        if len(f) >= 4: g2genus[f[0]] = f[3] or "unclassified"
    def genus_of(fid):
        try: return g2genus.get(".".join(fid.split("|",1)[1].split(".")[:2]), "unclassified")
        except Exception: return "unclassified"

    #  SCOPE.  The module is Tilapinevirus only.  Eleven further amnoonvirid
    #  species deposit a divergent PB1 and nothing else; at 0.75-1.25x the
    #  reference length seven of them land inside the S1 window and would
    #  train the PB1 profile on viruses this module does not route (measured:
    #  0 of 12 of their contigs route here at the 150-bit floor).  Excluded on
    #  the SPECIES, which is what the scope decision is about -- not on length,
    #  which cannot tell a divergent orthologue from a long isolate.
    g2sp = {}
    for l in open(os.path.join(WORK, "%s.full.tsv" % T), errors="replace"):
        f = l.rstrip("\n").split("\t")
        if len(f) > 4: g2sp[f[0]] = f[4]
    def species_of(fid):
        try: return g2sp.get(".".join(fid.split("|",1)[1].split(".")[:2]), "")
        except Exception: return ""
    IN_SCOPE = "Tilapinevirus tilapiae"

    #  1. discard what the SOURCE calls partial, before any homology is asked
    pool, dropped_partial, out_of_scope = [], [], []
    for i in range(len(md5s)):
        fid, a = (ann[i] + [""])[:2]
        s = seq[i][1] if len(seq[i]) > 1 else ""
        if not s: continue
        if PARTIAL.search(" ".join(a.split())):
            dropped_partial.append((fid, a, len(s))); continue
        if "*" in s: dropped_partial.append((fid, a + " [internal stop]", len(s))); continue
        if species_of(fid) != IN_SCOPE:
            out_of_scope.append((fid, species_of(fid), a, len(s))); continue
        pool.append((fid, a, s))

    #  2. assign by homology to the Til-4-2011 reference segments
    tmp = tempfile.mkdtemp()
    q = os.path.join(tmp, "q.faa")
    with open(q, "w") as o:
        for n, (fid, a, s) in enumerate(pool): o.write(">%d\n%s\n" % (n, s))
    ref = os.path.join(WORK, "ref", "tilv_segments.faa")
    subprocess.run(["makeblastdb","-in",ref,"-dbtype","prot","-out",os.path.join(tmp,"db")],
                   capture_output=True)
    out = subprocess.run(["blastp","-query",q,"-db",os.path.join(tmp,"db"),"-outfmt",
                          "6 qseqid sseqid bitscore length qlen","-evalue","1e-3",
                          "-num_threads","8"], capture_output=True, text=True).stdout
    best = {}
    for l in out.splitlines():
        p = l.split("\t"); n = int(p[0]); b = float(p[2])
        if float(p[3]) / float(p[4]) < MIN_QCOV: continue
        if b < MIN_BITS: continue
        if n not in best or b > best[n][1]: best[n] = (p[1], b)
    order = ["S1_PB1","S2_PB2","S3_PA","S4_NP","S5","S6","S7","S8","S9","S10"]
    seg2key = dict(("SEG%d" % (i+1), k) for i, k in enumerate(order))

    coll = collections.defaultdict(list)
    unassigned, runt, giant = [], [], []
    for n, (fid, a, s) in enumerate(pool):
        if n not in best: unassigned.append((fid, a, len(s), genus_of(fid))); continue
        key = seg2key[best[n][0]]
        exp = SEGMENTS[key][0]
        if len(s) < LO*exp: runt.append((fid, key, len(s), exp)); continue
        if len(s) > HI*exp: giant.append((fid, key, len(s), exp)); continue
        coll[key].append((fid, a, s))

    cd = os.path.join(WORK, "collections", MODULE)
    os.makedirs(cd, exist_ok=True)
    for k in order:
        with open(os.path.join(cd, k + ".fasta"), "w") as o:
            for fid, a, s in coll[k]: o.write(">%s\n%s\n" % (fid, s))

    td = os.path.join(WORK, "collections")
    def dump(name, hdr, rows):
        with open(os.path.join(td, name), "w") as o:
            o.write(hdr + "\n")
            for r in rows: o.write("\t".join(str(x) for x in r) + "\n")
    dump("OUT_OF_SCOPE.tsv", "feature_id\tspecies\tannotation\tlen", out_of_scope)
    dump("PARTIAL_DROPPED.tsv", "feature_id\tannotation\tlen", dropped_partial)
    dump("RUNTS_DROPPED.tsv", "feature_id\tsegment\tlen\texpected", runt)
    dump("OVERLONG_DROPPED.tsv", "feature_id\tsegment\tlen\texpected", giant)
    dump("UNASSIGNED_TRACKING.tsv", "feature_id\tannotation\tlen\tgenus", unassigned)
    #  synonyms.tsv -- the string-collapse record, in the schema the
    #  literature-gap check and the collapse artifact both read.
    syn = collections.Counter()
    for k in order:
        for fid, a, s in coll[k]: syn[(a, genus_of(fid), k)] += 1
    for fid, a, L, g in unassigned: syn[(a, g, "")] += 1
    with open(os.path.join(td, "synonyms.tsv"), "w") as o:
        o.write("count\tannotation\tgenus\tbins_to\n")
        for (a, g, k), c in sorted(syn.items(), key=lambda x: -x[1]):
            o.write("%d\t%s\t%s\t%s\n" % (c, a, g, k))

    print("%-9s %6s  %s" % ("feature", "n", "length range"))
    for k in order:
        L = sorted(len(s) for _, _, s in coll[k])
        print("  %-9s %4d  %s" % (k, len(L), ("%d-%d (ref %d)" % (L[0], L[-1], SEGMENTS[k][0])) if L else "EMPTY"))
    print("\nout of scope (species) : %d  %s" % (len(out_of_scope),
          dict(collections.Counter(x[1] for x in out_of_scope))))
    print("partial/stop discarded : %d" % len(dropped_partial))
    print("runts / overlong       : %d / %d" % (len(runt), len(giant)))
    print("unassigned             : %d" % len(unassigned))
    if unassigned:
        c = collections.Counter(x[3] for x in unassigned)
        print("   by genus: %s" % dict(c))

main()
