#!/usr/bin/env python3
"""Nimanivirus -- Macrobrachium rosenbergii Golda virus (MrGV), the ronivirus
of the giant freshwater prawn.  One species, Nimanivirus lahi; two variants,
LH1-2018 (Bangladesh, PMID 33023199) and Mr-18 (China, PMID 35011135), 97.8%
identical.

    5'-[--- ORF1a ---][--- ORF1b ---]--[ORF2]--[--- ORF3 ---]-3'
                    \\__ -1 PRF __/      125 aa     1503 aa

      ORF1a    4933 aa
      ORF1ab   7797 aa      -1 frameshift, slippage 5'-GGGUUUU-3' (PMID 35011135)
      ORF2      125 aa      nucleocapsid
      ORF3     1503 aa      glycoprotein precursor

    Genome 28,897-29,110 nt -- the largest in the family.

### Why this is NOT part of the Okavirus module

The ICTV profile records that MrGV lacks the triple-membrane-spanning p25
protein at the N-terminus of the okavirus pp3, and lacks ORFX.  So the
protein SET differs, which is what partitions a module.

That is a statement about alignments, and it was re-measured here.  MrGV ORF3
(1503 aa) against YHV pp3 (1666 aa) returns **no meaningful alignment at all**
-- the best HSP is 20.4 bits, 23.7% identity over 38 residues, which is noise.
The two precursors are not alignable, so no okavirus cleavage site can be
transferred to this one, and no PSSM could serve both.

Routing confirms the split independently: measured at the annotator's own
blastn settings against both modules' references at once, all 44 Nimanivirus
contigs route here and none to Okavirus.

### ORF3 ships WHOLE, and that is a documented limit

ICTV states MrGV ORF3 encodes homologues of gp116 and gp64 behind an
N-terminal signal peptide, with three predicted transmembrane domains.  But
no cleavage site is published for MrGV, and the okavirus sites cannot be
transferred (above).  So no mature peptide is derived: the precursor is
declared whole and the gp116/gp64 products are recorded as a gap rather than
cut at an invented coordinate.
"""
import os, sys
sys.path.insert(0, os.environ["LOWVAN_KIT"] + "/skill/scripts")
from collection_engine import run
HERE = os.path.dirname(os.path.abspath(__file__))

RULES = [
    ("PP1AB", r"^\s*(orf1ab(\s+protein)?|replicase\s+polyprotein\s*1ab|pp1ab)\s*$"),
    ("PP1A",  r"^\s*(orf1a(\s+protein)?|replicase\s+polyprotein\s*1a|pp1a)\s*$"),
    ("N",     r"^\s*(orf2(\s+protein)?|nucleocapsid(\s+protein)?|nucleoprotein|p20)\s*$"),
    ("ORF3",  r"^\s*(orf3(\s+protein)?|glycoprotein|structural\s+glycoprotein|pp3)\s*$"),
]

NOT_MODELLED = [
    (r"^\s*(orf1b(\s+protein|\s+polyprotein)?|replicase\s+polyprotein\s*1b)\s*$",
     "ORF1b is never translated on its own -- it is reached only by the -1 "
     "frameshift, as part of pp1ab. These 11 deposits are a submitter "
     "convention for the ORF, not a protein the virus makes"),
]

#  Every full-length deposit in this genus is at ONE length -- 4933 / 7797 /
#  125 / 1503 -- because the two variants are 97.8% identical. The windows are
#  generous around those so a future divergent isolate is not excluded.
WINDOW = {
    "PP1A":  (4400, 5400),
    "PP1AB": (7000, 8500),
    "N":     (110, 145),
    "ORF3":  (1350, 1650),
}

run(HERE, "N", "Nimanivirus", RULES, NOT_MODELLED=NOT_MODELLED, WINDOW=WINDOW)
