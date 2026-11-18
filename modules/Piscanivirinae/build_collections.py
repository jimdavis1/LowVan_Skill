#!/usr/bin/env python3
"""Piscanivirinae (Bafinivirus + Oncotshavirus) collections.

Genome organisation (+ssRNA, 26-27 kb), established for white bream virus by
Schuetze 2006 (PMID 16987966) and for fathead minnow nidovirus by Batts 2012
(PMID 22422065): 5'-ORF1a-ORF1b-S(ORF2)-M(ORF3)-N(ORF4)-3'. No HE, no E, no
accessory ORFs. ORF1b is expressed only through the -1 ribosomal frameshift
at the UUUAAAC slippery heptamer (NC_008516 misc_feature 14549..14555;
CDS join(906..14555,14555..21523)), so pp1ab is declared
`special: transcript_edit` and ORF1b is never a feature of its own.

One module for both genera: measured protein identity between them is
96.7-98.8% for the replicase and the layout is identical, so they share a
gene layout in the sense module-partitioning.md means.

Windows, from the string x genus histogram (modes, not medians):
    ORF1A  4400-5000   WBV 4555, FHMNV 4860, CSBV/OnV 4504-4675
    S      1150-1250   WBV 1220, FHMNV 1190, CSBV 1205
    M       215-250    225-245
    N       155-200    161-178

LENGTH_ROUTE on ORF1A is load-bearing: three deposits name an ORF1b product
"replicase polyprotein 1a" (e.g. KJ681496 15037..21879, 2280 aa) and two name
ORF1a "replicase polyprotein 1ab". The string names the wrong protein; the
length separates them, which is the documented remedy for a bimodal bin.
"""
import os, sys
sys.path.insert(0, os.environ["LOWVAN_KIT"] + "/skill/scripts")
from collection_engine import run
HERE = os.path.dirname(os.path.abspath(__file__))

RULES = [
    ("ORF1A", r"^\s*(replicase\s+polyprotein\s+(pp)?1a|polyprotein\s+1a|orf\s?1a|pp1a)\s*$"),
    ("S",     r"^\s*(spike(\s+glyco)?(\s*protein)?|s(\s+protein)?)\s*$"),
    ("M",     r"^\s*(membrane\s+protein|m(\s+protein)?)\s*$"),
    ("N",     r"^\s*(nucleocapsid\s+protein|nucleoprotein|n(\s+protein)?)\s*$"),
]

NOT_MODELLED = [
    (r"^\s*(replicase\s+polyprotein\s+(pp)?1ab|polyprotein\s+1ab|orf\s?1ab|1ab(\s+protein)?|pp1ab)\s*$",
     "pp1ab: -1 ribosomal frameshift product, called from transcript_edit references, no PSSM"),
    (r"^\s*(orf\s?1b|polyprotein\s+1b|rdrp|replicase\s+1b)\s*$",
     "ORF1b is expressed only inside pp1ab (transcript_edit); it has no start codon of its own"),
    (r"^\s*helicase\s*$",
     "294-327 aa ORF1b sub-fragments deposited as free-standing proteins; part of pp1ab"),
    (r"^\s*putative\s+polyprotein\s+1c\s*$",
     "one 742 aa hagfish bafinivirus deposit overlapping the ORF1b region; no such protein is described in the literature"),
    (r"^\s*putative\s+envelope\s+protein\s*$",
     "one 110 aa deposit from one isolate (KY130432); no envelope protein is described for bafiniviruses (PMID 16987966)"),
]

WINDOW = {"ORF1A": (4400, 5000), "S": (1150, 1250), "M": (215, 250), "N": (155, 200)}

if __name__ == "__main__":
    run(here=HERE, T="Pisc", MODULE="Piscanivirinae", RULES=RULES,
        NOT_MODELLED=NOT_MODELLED, WINDOW=WINDOW)
