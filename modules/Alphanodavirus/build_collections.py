#!/usr/bin/env python3
"""Alphanodavirus -- the insect nodaviruses (FHV, BBV, Boolarra, Nodamura,
Pariacoto, Newington).

Bipartite (+)ssRNA. The layout that defines this module, and separates it from
every other nodavirus module:

    RNA1   5'-[------------ protein A (RdRp) ------------]-3'
                                        [--- B1 ---]            same frame as A
                                        [--- B2 ---]            +1 frame vs A
           B1 and B2 are translated from subgenomic RNA3 (~370-414 nt), made
           from the 3' end of RNA1 during replication.

    RNA2   5'-[--------- alpha (capsid precursor) ---------]-3'
              [-------- beta --------][-- gamma --]
           alpha cleaves ITSELF after assembly: Asn | Ala, 44 residues from the
           C-terminus (Asn | Ser, 40 residues, in Pariacoto).

### Why this is its own module, and not merged with Betanodavirus

Both carry A + B1 + B2 + a capsid on two segments, so on a gene-count reading
they look identical. They are not: **the alphanodavirus capsid is cleaved and
the betanodavirus capsid is not** (PMID 9672600 -- the betanodavirus RNA2 ORF
"carries the catalytic D-75 residue but lacks the site for autocatalytic
proteolysis"). A module declares a protein SET, so a module holding both would
have to declare beta and gamma for genomes that do not make them, and would
flag every betanodavirus genome for two missing mature peptides. The budget
gate agrees that splitting costs nothing: measured 21 Nov 2026, Alphanodavirus
alone is 31 genomes in 10 clusters routing 100% at the 25-reference budget, and
Betanodavirus 342 genomes in 6 clusters, also 100%.

### The cleavage site is MEASURED here, not assumed at residue 363

The literature states the FHV site as Asn363 | Ala364 (PMID 1404613, PMID
9765417), and that is exactly right for FHV -- but FHV's alpha is 407 aa, and
quoting the residue NUMBER rather than the chemistry would misplace the cut in
every other species. Measured over the distinct alpha deposits in this dump,
22 Nov 2026, scanning for Asn followed by Ala/Ser with 38-52 residues left:

    Alphanodavirus flockense    (FHV)  alpha 407   N363|A   gamma 44
    Alphanodavirus heteronychi  (BBV)  alpha 407   N363|A   gamma 44
    American nodavirus pp9ad           alpha 407   N363|A   gamma 44
    Alphanodavirus HB-2007/CHN         alpha 407   N363|A   gamma 44
    Alphanodavirus boolarraense        alpha 403   N359|A   gamma 44
    Alphanodavirus nodamuraense (NoV)  alpha 399   N355|A   gamma 44
    Newington virus                    alpha 395   N351|A   gamma 44
    Alphanodavirus pariacotoense (PaV) alpha 401   N361|S   gamma 40

**Exactly one candidate site in every species, and no species has two.** Gamma
is 44 aa in seven of eight; Pariacoto is 40 aa at an Asn|Ser site, which is the
documented exception (PMID 10799587: "the amino acid sequence at the cleavage
site, which is Asn/Ala in all other alphanodaviruses, is Asn/Ser in PaV"). So
BETA and GAMMA are cut by chemistry on each genome's own alpha, and the
derivation reproduces the published FHV numbers without being given them.

### B1 and B2 are translation products, not cleavage products

PMID 15163762: "Protein B1 is in the same reading frame as protein A and thus
represents its C-terminal fragment", "whereas protein B2 is in the +1 reading
frame relative to protein A". Neither is cut out of A -- both initiate on
sgRNA3 -- so both are `CDS`, not `mat_peptide`.

BV-BRC has almost nothing for either: **3 B1 deposits and 9 B2 deposits in the
whole genus**, which is far below what a profile needs. Both collections are
therefore seeded from the deposits and expanded by extraction from the genomes
(`extract_by_ref.py`), with a structural gate on B1 that the deposits alone
could never provide: **a B1 is accepted only if it is an exact C-terminal
substring of protein A from its own genome.** 13 of 13 extractions passed.

### What is deliberately NOT modelled

- **The Pariacoto minor capsid protein** (PMID 13679619), a real virion protein
  initiated at the SECOND AUG of the capsid ORF and 24 residues shorter than
  alpha. It is alpha minus an N-terminal run, so it is 100% identical to alpha
  over its entire length and no PSSM can separate the two -- the same nested
  alternative-start trap the skill documents for VP1/VP2. Recorded in
  NOT_CALLABLE.md rather than declared and never called.
- **The BBV RNA2 second ORF** ("URF", 72 aa, base 1110 of NC_002037). PMID
  6548308 calls it "a second open reading frame for a putative protein 72 amino
  acid residues long" -- a 1984 computational prediction with no product
  reported in the 40 years since. One deposit exists. Excluded by name below.
"""
import os, sys
sys.path.insert(0, os.environ["LOWVAN_KIT"] + "/skill/scripts")
from collection_engine import run
HERE = os.path.dirname(os.path.abspath(__file__))

RULES = [
    #  --- RNA1 -------------------------------------------------------------
    ("A",     r"^\s*(protein\s+a|protein\s+a\s*\(rdrp\)|"
              r"put\.\s*protein\s+a(\s*\(m\.w\.\s*112\s*kd\))?|"
              r"rna[\s\-]?dep[ea]nd[ae]nt\s+rna[\s\-]?polymerase(\s+protein\s+a)?|"
              r"rna[\s\-]?dependant\s+rna[\s\-]?polymerase|"
              r"rna[\s\-]?directed\s+rna\s+polymerase|rdrp|rdrp\s+polymerase|"
              r"rna\s+polymerase|polymerase|putative\s+polymerase|"
              r"putative\s+rna[\s\-]?dep[ea]nd[ae]nt\s+rna[\s\-]?polymerase|"
              r"putative\s+rdrp|replicase|putative\s+replicase)\s*$"),
    ("B1",    r"^\s*(put\.\s*protein\s+b1|protein\s+b1|b1\s+protein|b1|"
              r"non-?structural\s+protein\s+b1)\s*$"),
    ("B2",    r"^\s*(protein\s+b2|b2\s+protein|b2|putative\s+protein\s+b2|"
              r"b2-?like\s+protein|protein\s+b|non-?structural\s+protein\s+b2)\s*$"),

    #  --- RNA2 -------------------------------------------------------------
    ("ALPHA", r"^\s*((capsid|coat)\s+pr[eo]cursor\s+protein(\s+alpha)?|"
              r"(capsid|coat)\s+protein\s+precursor(\s+alpha)?|"
              r"putative\s+(capsid|coat)\s+protein\s+precursor(\s+alpha)?|"
              r"putative\s+(capsid|coat)\s+precursor\s+protein(\s+alpha)?|"
              r"protein\s+alpha|capsid\s+protein\s+alpha|alpha\s+protein|"
              r"(capsid|coat)\s+protein|capsid|putative\s+capsid(\s+protein)?)\s*$"),
    ("BETA",  r"^\s*(mature\s+)?(capsid|coat)\s+protein\s+beta\s*$"),
    ("GAMMA", r"^\s*(mature\s+)?(capsid|coat)\s+protein\s+gamma\s*$"),
]

#  Every window is the literature's, with its PMID in length_refs.json. None is
#  a deposited median: this genus has 31 genomes and a median drawn from them
#  would be set by whichever isolate happened to be sequenced most.
LENGTH_ROUTE = {
    "A":     (930, 1060,
              "protein-A-named deposit outside the RdRp range (PaV 973 "
              "PMID 10799587; FHV/BBV/Boolarra 998; NoV 1043). A deposit far "
              "below this under a polymerase string is a partial RNA1"),
    "B1":    (95, 150,
              "B1-named deposit outside B1 length (FHV/BBV 102 deposited; NoV "
              "B1 begins at protein A residue 908, PMID 15163762). B1 is the "
              "C-terminal fragment of A, so anything near A's own length is A"),
    "B2":    (80, 150,
              "B2-named deposit outside B2 length (FHV 106 and NoV 134/137 "
              "stated in PMID 15163762; PaV 90 in PMID 10799587)"),
    "ALPHA": (380, 420,
              "capsid-named deposit outside alpha length (BBV/FHV 407 PMID "
              "6548308; PaV 401 and Boolarra 403; NoV 399; Newington 395). "
              "A capsid deposit BELOW this window is a partial RNA2, and one "
              "at ~363 would be mature beta -- a different feature"),
    "BETA":  (340, 370,
              "beta-named deposit outside mature beta length (FHV beta = alpha "
              "1-363 = 363 aa, PMID 9765417; the other species' alphas are "
              "4-12 residues shorter, so beta tracks them down to ~351)"),
    "GAMMA": (36, 52,
              "gamma-named deposit outside gamma length (FHV gamma = alpha "
              "364-407 = 44 aa, PMID 9765417; PaV gamma is 40 aa at the Asn/Ser "
              "site, PMID 10799587)"),
}

NOT_MODELLED = [
    (r"^\s*urf\s*$",
     "the BBV RNA2 second open reading frame, 72 aa from base 1110 of "
     "NC_002037 (PMID 6548308). A 1984 computational ORF call with no product "
     "reported since, and one deposit family-wide. Not a protein until a paper "
     "says so"),
    (r"^\s*(hypothetical\s+protein(\s+[12])?|unknown(\s+protein)?|"
     r"unnamed\s+protein\s+product|viral\s+protein|orf\s*[0-9]*)\s*$",
     "no product claim at all. Several of these are full-length alpha and are "
     "recovered by the homology rescue, which is the only treatment that "
     "cannot invent a homology group out of a submitter's silence"),
]

#  Chemistry at each terminus, measured on this dump's own alphas (see the
#  module docstring). alpha, A, B1 and B2 all initiate at their own Met.
#  beta starts at alpha's Met and ends at the Asn of the Asn|[AS] site;
#  gamma starts at the Ala (or Ser in PaV) immediately after it and runs to
#  alpha's own stop.
JUNCTION = {
    "A":     ("M",  None),
    "B1":    ("M",  None),
    "B2":    ("M",  None),
    "ALPHA": ("M",  None),
    "BETA":  ("M",  "N"),
    "GAMMA": ("AS", None),
}

WINDOW = {k: (v[0], v[1]) for k, v in LENGTH_ROUTE.items()}

if __name__ == "__main__":
    run(here=HERE, T="N", MODULE="Alphanodavirus", RULES=RULES,
        NOT_MODELLED=NOT_MODELLED, LENGTH_ROUTE=LENGTH_ROUTE, WINDOW=WINDOW,
        JUNCTION=JUNCTION)
