#!/usr/bin/env python3
"""Okavirus -- gill-associated virus (GAV), yellow head virus (YHV-1) and
YHV-8, the shrimp roniviruses.  Order Nidovirales, suborder Ronidovirineae.

    5'-[--- ORF1a (pp1a) ---][--- ORF1b ---]--[ORF2]--[--- ORF3 (pp3) ---][ORF4]-3'
                           \\__ -1 PRF __/      p20            |
                                                              +-- p25 + gp116 + gp64
                                        ORFX overlaps ORF3, +3 nt from its AUG

      pp1a    4060 aa (GAV) / 4070-4076 aa (YHV)
      pp1ab   6673 aa (GAV) / 6686-6692 aa (YHV)   -1 frameshift, 23-25%
      p20      144 aa (GAV) /  146 aa (YHV)        ORF2 nucleocapsid
      pp3     1640 aa (GAV) / 1656 (YHV-8) / 1666 aa (YHV-1)
      ORF4      83 aa  GAV ONLY
      px        86 aa  ORFX, all okaviruses

    The nucleoprotein gene sits UPSTREAM of the glycoprotein gene, which is
    unique among nidoviruses (PMID 15280504), and there is no discrete M gene.

### The ORF3 cleavage sites are published AND were re-measured here

PMID 12655087 placed the two cuts between Ala228/Thr229 and Ala1127/Leu1128
of the 1666 aa YHV pp3.  Checked against every full-length pp3 in this dump,
both sites verify exactly in the five 1666 aa isolates:

      site 1   ...LCFFITPAFA | TILSGI...      A-F-A, cut after A228
      site 2   ...SIISPASA   | LAPRQ...       A-S-A, cut after A1127

Both are the [Ala-X-Ala] signal-peptidase-type-1 motif the ICTV profile
describes, with the cut after the second Ala.

**The residue numbers do not travel and the motifs do.**  A 162 nt indel in
ORF3 separates YHV subtypes, so pp3 runs 1640-1666 aa and site 2 sits at
1101 in GAV, 1116 in YHV-8 and 1127 in YHV-1.  Measured across all ten
okavirus pp3 deposits, `ITP.FA` finds site 1 in 8 of 10 (at 228 in six of
them) and the `A[ST]A|LA` motif finds site 2 in 8 of 10.  Cutting on a
column number would mis-cut every genotype but one.

### What is NOT declared, and why

The ICTV profile states that the number and composition of the nsps released
from pp1ab "is not yet known".  BV-BRC holds deposits annotated `3C-like
protease` (335-343 aa), `helicase` (660-846 aa) and `RNA polymerase`
(92-286 aa), but these are DOMAIN fragments with no published junctions --
only the 3CLpro's own cleavage specificity VxHE|(L/V) is characterised
(PMID 12502857), and a specificity is not a coordinate.  So no mature
peptide is declared for pp1a/pp1ab.  Declaring one would mean inventing the
boundaries, which is the one thing the source annotation ceiling forbids.
"""
import os, sys
sys.path.insert(0, os.environ["LOWVAN_KIT"] + "/skill/scripts")
from collection_engine import run
HERE = os.path.dirname(os.path.abspath(__file__))

RULES = [
    ("PP1AB", r"^\s*(orf1ab\b.*|replicase\s+polyprotein\s*1ab|polyprotein\s*pp?1ab|pp1ab|1ab)\s*$"),
    ("PP1A",  r"^\s*(orf1a(\s+polyprotein)?(\s+pp1a)?|replicase\s+polyprotein\s*1a|pp1a)\s*$"),
    ("N",     r"^\s*(p20(\s+nucleocapsid(\s+protein)?)?|nucleocapsid(\s+protein)?|"
              r"n\s+protein|nucleoprotein)\s*$"),
    ("PP3",   r"^\s*((structural|envelope\s+structural|spike|mutant\s+envelope)\s+glycoprotein|"
              r"polyprotein\s*(pp3|precursor)|pp3|orf3(\s+protein)?)\s*$"),
    ("ORF4",  r"^\s*(orf4(\s+protein)?|unknown)\s*$"),
]

#  Deposits that are real proteins but cannot be modelled as declared features.
NOT_MODELLED = [
    (r"^\s*(3c-?like\s+protease|helicase|rna\s+polymerase|"
     r"rna-dependent\s+rna\s+polymerase)\s*$",
     "a DOMAIN of pp1a/pp1ab, not a declared mature peptide. The ICTV profile "
     "states the number and composition of the nsps released from pp1ab is not "
     "yet known; only the 3CLpro cleavage SPECIFICITY VxHE|(L/V) is published "
     "(PMID 12502857), and a specificity is not a coordinate. Declaring a "
     "mature peptide here would mean inventing its boundaries"),
    (r"^\s*(orf1b(\s+polyprotein|-like\s+protein)?|replicase\s+polyprotein\s*1b)\s*$",
     "ORF1b is never translated on its own -- it is reached only by the -1 "
     "frameshift, as part of pp1ab (PMID 10811931). These deposits are a "
     "submitter convention for the ORF, not a protein the virus makes"),
    (r"^\s*(glycoprotein\s*(116|64)|gp\s*(116|64))\s*$",
     "fragments of pp3 (267-338 aa against a true 899 aa gp116 and 539 aa "
     "gp64). The real products are DERIVED by cutting pp3 at its two measured "
     "[Ala-X-Ala] sites, not binned from these partial deposits"),
    (r"^\s*non-structural\s+protein\s*$",
     "names no homology group in this genus"),
]

#  Literature windows.  pp1a/pp1ab/pp3/p20 ranges are the measured full-length
#  deposits of GAV, YHV-1 and YHV-8; ORF4 is GAV's 83 aa (PMID 10811931 /
#  Walker & Sittidilokratna, Encyclopedia of Virology 2008).
WINDOW = {
    "PP1A":  (3800, 4300),
    "PP1AB": (6300, 7000),
    "N":     (130, 160),
    "PP3":   (1500, 1750),
    "ORF4":  (70, 110),
}

run(HERE, "O", "Okavirus", RULES, NOT_MODELLED=NOT_MODELLED, WINDOW=WINDOW)
