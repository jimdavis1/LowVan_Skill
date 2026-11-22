#!/usr/bin/env python3
"""Simarterivirinae module -- the EXPANDED arterivirus genome layout.

The simian arteriviruses (SHFV and relatives: Delta-, Epsilon-, Eta-, Zeta-,
Theta- and Iotaarterivirus plus Southwest baboon virus) carry the canonical
minor-glycoprotein cassette TWICE. Layout, confirmed from sequence on all 171
near-complete genomes rather than from the deposited strings:

    5'-ORF1a-ORF1b-[GP2',E',GP3',GP4']-[E,GP2,GP3,GP4]-ORF5a-GP5-M-N-3'

That second cassette is the whole reason this is a separate module from
`Arteriviridae`: no other arterivirus has it.

### How the structural proteins are binned, and why not by string

They are NOT binned by string. The commonest product string in this subfamily
is `minor glycoprotein` (359 deposits) and it is used for GP2, GP3, GP4 and
all four primed duplicates; five mutually incompatible naming schemes are in
circulation, including inside RefSeq (lit/structural.md section 5.2). GP2 and
GP2' are paralogues, so homology alone cannot separate them either.

`assign_cassette.py` aligns each genome's ordered protein list against the
canonical ordered feature list, scoring each pairing by blastp bitscore
against a labelled seed set -- position curates the orthogroups, the sequence
score calls them. Checked against the 1,177 deposits whose own string IS
unambiguous: **1,152 agree, 97.9%**. The single disagreeing sequence is a
214 aa protein deposited as `E protein` in 25 genomes, where every genuine E
in the subfamily is 75-81 aa; that is the ORF2a/ORF2b naming inversion
(PMID 10400725 vs 11504553) and the positional call is the correct one.

Only ORF1a is binned from strings here, and only because its strings are
unambiguous once the length window is applied.

### Mature peptides: none, and that is a literature ceiling

There is **no cleavage-site mapping of any kind for any simarterivirus**
beyond the disputed nsp1 region, and no simarterivirus protease has been
expressed or assayed. Projecting the EAV/PRRSV nsp map across 35-40% identity
would be an inference, not a measurement, so this module stops at the
polyprotein. nsp1 itself is disputed between two primary papers that both
measured it (PMID 24899184 three products, PMID 24928046 two) and is not
declared either.
"""
import os, shutil, subprocess, sys
sys.path.insert(0, os.environ["LOWVAN_KIT"] + "/skill/scripts")
from collection_engine import run
HERE = os.path.dirname(os.path.abspath(__file__))

RULES = [
    ("ORF1A", r"^\s*(polyprotein|1a\s+replicase(\s+protein)?|orf\s?1a(\s+polyprotein)?|"
              r"nonstructural\s+polyprotein\s+1a|replicase\s+protein\s+1a|"
              r"polyprotein\s+1a|pp1a)\s*$"),
]

LENGTH_ROUTE = {
    "ORF1A": (1900, 2200, "a `polyprotein`-named deposit outside the simarterivirus pp1a "
                          "length (2,048-2,105 aa across the subfamily): the ORF1b or pp1ab "
                          "product, both of which are transcript_edit features with no PSSM"),
}
WINDOW = {k: (v[0], v[1]) for k, v in LENGTH_ROUTE.items()}
JUNCTION = {"ORF1A": ("M", None)}

NOT_MODELLED = [
    (r"^\s*(orf\s?1ab(\s+protein)?|pp1ab|nonstructural\s+polyprotein\s+(1ab|pp1ab)|"
     r"replicase\s+1ab|1ab\s+replicase(\s+protein)?)\s*$",
     "pp1ab: -1 ribosomal frameshift product at the U_UUA_AAC heptamer; called from "
     "nucleotide references as special:transcript_edit, no PSSM"),
    (r"^\s*(orf\s?1b(\s+polyprotein)?|putative\s+(protein\s+)?1b|putative\s+1b\s+protein|"
     r"pp1b|1b)\s*$",
     "ORF1b has no initiation codon of its own; its product exists only inside pp1ab"),
    (r"^\s*(nsp2tf(\s+protein)?|tf\s+protein|orf1atf\s+polyprotein|"
     r"(putative\s+)?(protein\s+)?nsp2f(\s+protein)?)\s*$",
     "nsp2TF: -2 ribosomal frameshift product of the nsp2 region (PMID 23043113); "
     "a transcript_edit feature, no PSSM"),
]

if __name__ == "__main__":
    run(here=HERE, T="Sim", MODULE="Simarterivirinae", RULES=RULES,
        NOT_MODELLED=NOT_MODELLED, LENGTH_ROUTE=LENGTH_ROUTE, WINDOW=WINDOW,
        JUNCTION=JUNCTION)
    #  The cassette is assigned positionally, not by string, so it is produced
    #  here rather than by a rule. run() has just rewritten the collection
    #  directory, so the merge has to follow it in the same invocation --
    #  otherwise a later rebuild silently drops twelve features.
    subprocess.run([sys.executable, os.path.join(HERE, "assign_cassette.py")], check=True)
    cdir = os.path.join(HERE, "collections", "Simarterivirinae")
    for fn in sorted(os.listdir(os.path.join(HERE, "Extracted"))):
        if not fn.endswith(".faa"):
            continue
        key = fn.split(".")[1]
        shutil.copyfile(os.path.join(HERE, "Extracted", fn),
                        os.path.join(cdir, key + ".fasta"))
    #  run() wrote UNASSIGNED_TRACKING.tsv before the cassette existed, so it
    #  lists 1,943 features as uncalled that the positional assignment has
    #  just bound. That file is a GATE -- install_module.py refuses to install
    #  while a row of it reaches the member floor -- so leaving it stale would
    #  either block the module or, worse, be waved through by hand. Recompute
    #  it against what the collections actually hold.
    import collections as _c
    called = set()
    for fn in os.listdir(cdir):
        if fn.endswith(".fasta"):
            for line in open(os.path.join(cdir, fn)):
                if line.startswith(">"):
                    called.add(line[1:].split()[0].strip())
    ann = {}
    for line in open(os.path.join(HERE, "Sim.uniq.id_ann")):
        f = line.rstrip("\n").split("\t")
        ann[f[0]] = f[1] if len(f) > 1 else ""
    per = _c.Counter()
    for line in open(os.path.join(HERE, "Sim.id_md5")):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[1].strip():
            per[f[1].strip()] += 1
    md5s = [l.strip() for l in open(os.path.join(HERE, "Sim.uniq.md5"))]
    seqs = [l.rstrip("\n").split("\t") for l in open(os.path.join(HERE, "Sim.uniq.seq"))]
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(HERE, "Sim.uniq.id_ann"))]
    nm = set()
    nmp = os.path.join(HERE, "collections", "NOT_MODELLED.tsv")
    if os.path.exists(nmp):
        for line in open(nmp):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] != "features":
                nm.add(" ".join(f[1].split()))
    feat, uniq = _c.Counter(), _c.Counter()
    for i, m in enumerate(md5s):
        fid = rows[i][0] if rows[i] else ""
        a = " ".join((rows[i][1] if len(rows[i]) > 1 else "").split())
        if not fid or fid in called or a in nm:
            continue
        feat[a] += per.get(m, 1); uniq[a] += 1
    with open(os.path.join(HERE, "collections", "UNASSIGNED_TRACKING.tsv"), "w") as fh:
        fh.write("features\tuniq_seqs\tnot_already_called\tannotation\n")
        for a, c in feat.most_common():
            fh.write("%d\t%d\t%d\t%s\n" % (c, uniq[a], uniq[a], a))
    print("\n  cassette merged into collections/Simarterivirinae/")
    print("  unassigned recomputed against the merged collections: "
          "%d strings, %d features" % (len(feat), sum(feat.values())))
