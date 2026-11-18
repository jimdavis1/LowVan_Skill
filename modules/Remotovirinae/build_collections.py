#!/usr/bin/env python3
"""Remotovirinae (genus Bostovirus; bovine nidovirus 1) collections.

One genome paper: Tokarz 2015 (PMID 25918239), strain TCH5, KM589359,
20,261 nt -- seven ORFs: ORF1a, ORF1b, S (562 aa), M (231 aa), N (178 aa),
ORF5 "G2" (455 aa putative accessory glycoprotein) and an unnamed 88 aa ORF
near the 3' end. No HE. ORF1b is reached by a -1 frameshift at UUUAAAC
(KM589359 CDS join(487..9825,9825..15818)), so pp1ab is `transcript_edit`.

Windows from the string x genus histogram; this genus is tight (18 genomes,
one 70%-identity cluster) so every window is narrow by measurement, not by
assumption:
    ORF1A 3000-3200   3111-3114
    S       540-580   561-562
    M       225-240   231
    N       170-185   177-178
    G2      420-470   436-455
    ORF6     80-95    88
"""
import os, sys
sys.path.insert(0, os.environ["LOWVAN_KIT"] + "/skill/scripts")
from collection_engine import run
HERE = os.path.dirname(os.path.abspath(__file__))

RULES = [
    ("ORF1A", r"^\s*(polyprotein\s+pp1a|pp1a)\s*$"),
    ("S",     r"^\s*(spike\s+glycoprotein|s)\s*$"),
    ("M",     r"^\s*(membrane\s+protein\s*1|m1)\s*$"),
    ("N",     r"^\s*(nucleocapsid|n)\s*$"),
    ("G2",    r"^\s*glycoprotein\s+g2\s*$"),
    ("ORF6",  r"^\s*hypothetical\s+protein\s*$"),
]

NOT_MODELLED = [
    (r"^\s*(polyprotein\s+pp1b|pp1b)\s*$",
     "ORF1b is expressed only inside pp1ab (transcript_edit); no start codon of its own"),
    (r"^\s*pp1a/b\s*$",
     "pp1ab: -1 ribosomal frameshift product, called from transcript_edit references, no PSSM"),
]

WINDOW = {"ORF1A": (3000, 3200), "S": (540, 580), "M": (225, 240),
          "N": (170, 185), "G2": (420, 470), "ORF6": (80, 95)}

if __name__ == "__main__":
    run(here=HERE, T="Remo", MODULE="Remotovirinae", RULES=RULES,
        NOT_MODELLED=NOT_MODELLED, WINDOW=WINDOW)
