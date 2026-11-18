#!/usr/bin/env python3
"""Torovirinae (genus Torovirus) collections.

Genome organisation (+ssRNA, ~28 kb): 5'-ORF1a-ORF1b-S-M-HE-N-3'. ORF1b is
expressed only by a -1 ribosomal frameshift at the ORF1a/1b overlap, so the
pp1ab product is not collinear with the genome and is declared
`special: transcript_edit` (built from genome nucleotide references, no PSSM).
The literature for every protein is in ../notes/LITERATURE.md.

Windows (from the full-length mode, per string x genus histogram):
    ORF1A  4300-4700   pp1a: 4417-4418 (BToV/PToV), 4569 (BEV)
    S      1500-1650   1576-1587 (incl. signal peptide)
    M       220-245    233
    HE      400-440    416-430; Berne virus HE is a pseudogene (truncated) --
                       "truncated non-functional hemagglutinin esterase" is
                       discarded at triage like every truncated deposit
    N       155-175    160-167
"""
import os, sys
sys.path.insert(0, os.environ["LOWVAN_KIT"] + "/skill/scripts")
from collection_engine import run
HERE = os.path.dirname(os.path.abspath(__file__))

RULES = [
    ("ORF1A", r"^\s*(polyprotein\s*(pp)?1a|orf\s?1a|pp1a|replicase(\s+polyprotein)?(\s+1a)?|"
              r"nonstructural\s+protein)\s*$"),
    ("S",     r"^\s*(spike(\s+glyco)?(\s*protein)?|s(\s+protein)?|peplomer\s+protein.*|"
              r"spike\s+glycoprotein)\s*$"),
    ("M",     r"^\s*(m(\s+protein)?|membrane(\s+glyco)?\s*protein|membrance\s+protein|"
              r"envelope\s+protein)\s*$"),
    ("HE",    r"^\s*(he(\s+protein)?|hemagglutinin[\s\-]esterase(\s+protein)?)\s*$"),
    ("N",     r"^\s*(n(\s+protein)?|nucleocapsid(\s+phospho)?\s*protein|nucleocapsid\s+phosphoprotei|"
              r"nucleoprotein)\s*$"),
]

NOT_MODELLED = [
    (r"^\s*(polyprotein\s*1ab|1ab\s+protein|polyprotein|orf1ab|pp1ab)\s*$",
     "pp1ab: -1 ribosomal frameshift product, called by transcript_edit references, no PSSM"),
    (r"^\s*(1b\s+protein|orf1b|pp1b)\s*$",
     "ORF1b has no initiation codon; it is expressed only inside pp1ab (transcript_edit)"),
    (r"precursor\s+polypeptide", "one S record carrying its signal-peptide numbering; S is modelled"),
]

WINDOW = {
    "ORF1A": (4300, 4700),
    "S":     (1500, 1650),
    "M":     (220, 245),
    "HE":    (400, 440),
    "N":     (155, 215),   # widened: Bangali torovirus N is 211 aa and spans a 163 aa N end to end (not a fragment)
}

if __name__ == "__main__":
    run(here=HERE, T="Toro", MODULE="Torovirinae", RULES=RULES,
        NOT_MODELLED=NOT_MODELLED, WINDOW=WINDOW)
