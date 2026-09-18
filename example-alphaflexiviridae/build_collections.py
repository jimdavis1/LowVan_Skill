#!/usr/bin/env python3
"""Triage Alphaflexiviridae into feature collections.

Potexvirus, Allexivirus and Lolavirus in one module: replicase + triple gene
block + coat protein, with Allexivirus and Lolavirus adding a 40K serine-rich
protein and a nucleic-acid binding protein. Mandarivirus belongs here too and
arrives without being asked for -- BV-BRC labels Indian citrus ringspot virus
and Citrus yellow vein clearing virus as Potexvirus, and they are recovered by
homology rather than by the genus column. See TRIAGE_NOTES.md.

Seven features. The rules bin on (string, genus, length) together, because a
string-only rule is wrong for three of them and a length-only rule is wrong
for two:

    REP          1,200-1,800   replicase; the RdRp domain
    TGB1           200-280     triple gene block protein 1
    TGB2            85-135     triple gene block protein 2
    TGB3            50-120     triple gene block protein 3, Potexvirus
    CP             180-400     coat protein, both full-length forms
    ALLEXI_40K     320-400     40K serine-rich protein, Allexivirus/Lolavirus
    NABP            90-260     nucleic-acid binding protein, Allexi/Lolavirus
                                 and Mandarivirus

Three strings that must not be trusted on their own, all of them the U-number
trap -- a positional label covering non-homologous proteins:

  "tgb3" in Allexivirus is mostly NOT a TGB3. n=233, p5=99 but p25-p95 at
  357-364 aa against a potexvirus TGB3 of 63-89. Allexiviruses carry the 40K
  where potexviruses carry TGB3 and submitters label it by position. Binning
  on the string would call "triple gene block protein 3" on the 40K locus
  across the genus, so ALLEXI_40K is tested first and scoped to the genus.

  "rna dependent rna polymerase", unhyphenated, is bimodal at 160 and 1,417
  aa. The hyphenated form is clean.

  "polymerase" has median 230 aa -- TGB1-sized -- with real polymerases only
  in the tail.

The last two need no special rule: they classify as REP and the length window
drops them into REP.outliers.fasta, which is what outliers are for.

TGB2 (~103) overlaps NABP (~104-127), and CP overlaps TGB1 across 224-265.
Length cannot separate either pair. They are different proteins, so profiles
built from correctly labelled sequences separate them by homology; the
discriminator is bit_cutoff plus qc_cross_feature.py, and the windows stay
wide enough not to misdescribe the feature.
"""
import re, os, collections, argparse

MODULE = "Alphaflexiviridae"
DUMP   = "Alphaflexiviridae"
ALLEXI = ("Allexivirus", "Lolavirus")

def rx(p): return re.compile(p, re.I)

#  Order matters. ALLEXI_40K is first and genus-scoped so that "tgb3" on an
#  Allexivirus genome never reaches the TGB3 rule; NABP precedes CP so that
#  "nucleic acid binding protein" is not caught by a loose CP pattern; and
#  REP is last because "replicase"/"polymerase" are the least specific.
RULES = [
    ("ALLEXI_40K", ALLEXI, rx(r"^tgb ?3$|^tgbp ?3$|serine.rich|40 ?kda|^p4[02]$|^orf ?4$")),
    ("NABP", ALLEXI, rx(r"nucleic.acid.bi[dn]|nucleic acid biding|\bnabp\b|\bnbp\b|rna.bind")),
    ("NABP", None,   rx(r"nucleic.acid.bi[dn]|nucleic acid biding|\bnabp\b")),
    ("TGB1", None, rx(r"triple gene block (protein )?1|^tgbp ?1|^tgb ?1|25 ?kda")),
    ("TGB2", None, rx(r"triple gene block (protein )?2|^tgbp ?2|^tgb ?2|12 ?kda")),
    ("TGB3", None, rx(r"triple gene block (protein )?3|^tgbp ?3|^tgb ?3|6\.?4 ?kda|^7 ?kda")),
    ("CP",   None, rx(r"\bcoat\b|capsid|\bcp\b")),
    ("REP",  None, rx(r"replicas|replicat|polymeras|\brdrp\b|helicase|methyltransferase")),
]
FORBIDDEN = [rx(r"^hypothetical protein"), rx(r"^unknown"), rx(r"^orf ?\d+$"),
             rx(r"^unnamed protein product$"), rx(r"^putative protein$"),
             rx(r"^p\d+$"), rx(r"^\d+\.?\d* ?k(da)? protein$")]

EXPECTED = {
    "REP":        (1200, 1800),
    "TGB1":       ( 200,  280),
    "TGB2":       (  85,  135),
    "TGB3":       (  50,  120),
    "CP":         ( 180,  400),
    "ALLEXI_40K": ( 320,  400),
    "NABP":       (  90,  260),
}
NOT_MODELLED = [
    ("TGB1/2/3, CP (Alphaflexiviridae_noTGB)", "Botrexvirus, Platypuvirus, Sclerodarnavirus", 42,
     "Separate module: fungal and non-virion, replicase and CP with no triple "
     "gene block. 42 genomes, below the rep-contig budget; may not ship."),
]

def normalise(a):
    return re.sub(r"\s+", " ", re.sub(r"[-_]", " ", a)).strip()

def _match(ann, genus, scoped):
    for cand in (ann, normalise(ann)):
        for key, g, prx in RULES:
            if scoped and not g: continue
            if not scoped and g: continue
            if g and genus not in (g if isinstance(g, tuple) else (g,)): continue
            if prx.search(cand): return key
    return None

def genus_specific(ann, genus): return _match(ann, genus, True)
def classify(ann, genus):       return _match(ann, genus, False)
def forbidden(a):
    return any(p.search(a) or p.search(normalise(a)) for p in FORBIDDEN)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workdir", default="."); a = ap.parse_args()
    W = a.workdir
    g2g = {}
    for l in open(os.path.join(W, DUMP + ".id_name_gs")):
        p = l.rstrip("\n").split("\t"); g2g[p[0]] = (p[3] if len(p) > 3 else "") or "?"
    seq = {}
    for l in open(os.path.join(W, DUMP + ".uniq.seq")):
        p = l.rstrip("\n").split("\t")
        if len(p) > 1: seq[p[0]] = p[1]
    md5 = [l.strip() for l in open(os.path.join(W, DUMP + ".uniq.md5"))]
    ann = [l.rstrip("\n").split("\t") for l in open(os.path.join(W, DUMP + ".uniq.id_ann"))]

    out = collections.defaultdict(list); outl = collections.defaultdict(list)
    unass = collections.Counter(); syn = collections.Counter(); hit = collections.Counter()
    for m, row in zip(md5, ann):
        fid, aname = row[0], (row[1] if len(row) > 1 else "")
        gid = fid.split("|", 1)[1].rsplit(".", 2)[0]
        genus = g2g.get(gid, "?"); s = seq.get(m, "")
        if not s: continue
        key = genus_specific(aname, genus) or (None if forbidden(aname) else classify(aname, genus))
        syn[(aname, genus, key or "UNASSIGNED")] += 1
        if key is None:
            unass[(aname, genus)] += 1; continue
        lo, hi = EXPECTED[key]; rec = ">%s|%s|%s\n%s\n" % (fid, genus, m, s)
        if lo <= len(s) <= hi: out[key].append(rec); hit[key] += 1
        else: outl[key].append(rec)

    cdir = os.path.join(W, "collections", MODULE); os.makedirs(cdir, exist_ok=True)
    for key in EXPECTED:
        open(os.path.join(cdir, key + ".fasta"), "w").writelines(out.get(key, []))
        if outl.get(key):
            open(os.path.join(cdir, key + ".outliers.fasta"), "w").writelines(outl[key])
    cd = os.path.join(W, "collections")
    with open(os.path.join(cd, "synonyms.tsv"), "w") as f:
        f.write("count\tannotation\tgenus\tbins_to\n")
        for (s_, g, k), c in sorted(syn.items(), key=lambda x: -x[1]):
            f.write("%d\t%s\t%s\t%s\n" % (c, s_, g, k))
    with open(os.path.join(cd, "UNASSIGNED_TRACKING.tsv"), "w") as f:
        f.write("count\tannotation\tgenus\n")
        for (s_, g), c in sorted(unass.items(), key=lambda x: -x[1]):
            f.write("%d\t%s\t%s\n" % (c, s_, g))
    with open(os.path.join(cd, "NOT_MODELLED.tsv"), "w") as f:
        f.write("feature\ttaxon\tn_sequences\treason\n")
        for r in NOT_MODELLED: f.write("%s\t%s\t%d\t%s\n" % r)

    print("%-12s %9s %10s   %s" % ("feature", "in-range", "outliers", "window"))
    for key in EXPECTED:
        print("%-12s %9d %10d   %d-%d aa" % (key, hit[key], len(outl.get(key, [])), *EXPECTED[key]))
    print("\nunassigned: %d sequences (%d strings)   total unique: %d"
          % (sum(unass.values()), len(unass), len(md5)))
main()
