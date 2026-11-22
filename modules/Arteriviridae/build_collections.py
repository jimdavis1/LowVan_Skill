#!/usr/bin/env python3
"""Arteriviridae module -- the SIMPLE (canonical) arterivirus genome layout.

Covers every in-scope Arteriviridae genome EXCEPT Simarterivirinae, whose
duplicated minor-glycoprotein cassette is a different arrangement and a
different module. Layout, confirmed from sequence on all 2,401 near-complete
genomes (layout/call_layout.py, never from the deposited product strings):

    5'-ORF1a-ORF1b-[E,GP2,GP3,GP4]-ORF5a-GP5-M-N-3'

ORF1b has no start codon of its own: pp1ab is made by a -1 ribosomal
frameshift at a U_UUA_AAC (EAV: G_UUA_AAC) heptamer immediately upstream of
the ORF1a stop, so it is `special: transcript_edit` and gets no PSSM. The same
is true of nsp2TF and nsp2N, the -2/-1 frameshift products of the nsp2 region.

### The naming trap this file exists to survive

ORF labels in this family are NOT portable and are never used here as the sole
evidence (structural literature review, lit/structural.md):

  * EAV   calls the small protein ORF2a and the glycoprotein ORF2b (PMID 10400725)
  * PRRSV calls the glycoprotein ORF2a and the small protein ORF2b (PMID 11504553)
    -- opposite decisions on the same gene pair, and the convention tracks the
    depositing lab rather than the phylogeny.
  * "E" meant GP5 in early PRRSV-1 papers (PMID 7831770): an `E` string on a
    ~200 aa product is GP5, on a ~70 aa product it is the small protein.
  * `putative glycoprotein 1..4` (Chen 2023 deposits) is off by one throughout.

So every rule here is paired with a LENGTH window taken from the literature,
and a deposit whose string says one protein and whose length says another is
dropped to NOT_MODELLED with the reason recorded, never silently re-binned.
`putative ...` strings are sent to UNASSIGNED by the engine and arbitrated by
the homology rescue, which is the correct treatment for the off-by-one set.

Length windows (aa), from RefSeq CDS recomputation across 24 family reference
genomes plus the measured sizes in lit/structural.md section 2:

    ORF1A  1700-2700   EAV 1727 .. PRRSV-2 2503
    E        45-95     WPDV 45 .. SHFV 81 (PRRSV 73, EAV 67)
    GP2     150-320    WPDV 207, EAV 227, PRRSV 249-256, APRAV 304
    GP3     130-280    KRTGV 136 .. PRRSV-1 265
    GP4      80-210    WPDV 83 .. PRRSV-1 183
    ORF5A    40-80     PRRSV-1 43 .. WPDV 78
    GP5     160-285    WPDV 178 .. EAV 255
    M       155-190    EAV 162 .. WPDV 180
    N       100-140    EAV 110 .. APRAV 134

nsp windows are the EAV (NC_002532, measured map, PMID 9371590/9971783/17098961)
and PRRSV (PMID 20696193 Table 1, Edman/proteomics-anchored) subunit lengths.
"""
import os, sys
sys.path.insert(0, os.environ["LOWVAN_KIT"] + "/skill/scripts")
from collection_engine import run
HERE = os.path.dirname(os.path.abspath(__file__))

RULES = [
    #  --- replicase, polyprotein level -------------------------------------
    ("ORF1A", r"^\s*(polyprotein\s*(pp)?\s*1a|orf\s?1a(\s+polyprotein)?|pp1a|"
              r"replicase(\s+(polyprotein|protein))?(\s+1a)?|polyprotein\s+pp1a|"
              r"1a(\s+replicase(\s+protein)?)?|nonstructural\s+polyprotein\s+1a|"
              r"nonstructural\s+proteins?\s+orf1a|non-?structural\s+protein(\s+1a)?|"
              r"polyprotein\s+orf1a|nonstructural\s+polyprotein\s+orf1a|"
              r"(non-?structural|viral)\s+replicase|viral\s+replicase\s+polyprotein|"
              r"orf1a\s+protein|"
              r"replicase\s+orf1a\s+polyprotein|proteinase(\s+1a)?|polyprotein|"
              r"papain-like\s+cysteine\s+protease|cysteine/serine\s+protease|"
              r"replicase\s+protein\s+1a)\s*$"),
    #  --- nonstructural subunits (mature peptides) --------------------------
    ("NSP1",  r"^\s*(nsp1(\s*\(pcp1b\))?|nonstructural\s+protein\s+1)\s*$"),
    ("NSP1A", r"^\s*nsp1[\s\-]?alpha\s*$"),
    ("NSP1B", r"^\s*nsp1[\s\-]?(beta|b)\s*$"),
    ("NSP2",  r"^\s*(nsp2|non-?structural\s+protein\s+2|nonstructural\s+protein\s+nsp2|"
              r"cysteine\s+proteinase\s+2|cysteine\s+protease\s+nsp2)\s*$"),
    ("NSP3",  r"^\s*nsp3\s*$"),
    ("NSP4",  r"^\s*(nsp4|3c-like\s+serine\s+proteinase)\s*$"),
    ("NSP5",  r"^\s*nsp5\s*$"),
    ("NSP6",  r"^\s*nsp6\s*$"),
    ("NSP7A", r"^\s*nsp7[\s\-]?alpha\s*$"),
    ("NSP7B", r"^\s*nsp7[\s\-]?beta\s*$"),
    ("NSP8",  r"^\s*nsp8\s*$"),
    ("NSP9",  r"^\s*(nsp9|rna[\s\-]dependent\s+rna\s+polymerase)\s*$"),
    ("NSP10", r"^\s*(nsp10|metal-binding\s+and\s+ntpase/helicase\s+domains)\s*$"),
    ("NSP11", r"^\s*nsp11\s*$"),
    ("NSP12", r"^\s*nsp12\s*$"),
    #  --- structural --------------------------------------------------------
    ("E",     r"^\s*(e|2b|orf2b|orf4a|(small\s+)?(un)?glycosylated\s+envelope\s+protein(\s+e)?|"
              r"e\s+(envelope\s+)?protein|envelope\s+protein(\s+e)?|protein\s+e|"
              r"small\s+envelope\s+protein(\s+e)?|envelope\s+glycoprotein\s+e|"
              r"orf2b\s+protein|small\s+unglycosylated\s+(envelope\s+)?protein)\s*$"),
    ("GP2",   r"^\s*(gp2[ab]?(\s+(envelope\s+)?protein)?|gp2[ab]?\s+glyco(sylated)?(\s+envelope)?\s*protein|"
              r"glyco(sylated|protein)?\s*(membrane\s+)?protein\s*2[ab]?|glycoprotein\s+2\s+gp2|"
              r"orf2a?(\s+protein)?|orf2a?\s+envelope\s+protein|orf4b|envelope\s+glycoprotein\s+2|"
              r"envelope\s+protein\s+gp2|glycoprotein\s+gp2|glycosylated\s+envelope\s+protein\s+2|"
              r"minor\s+structural\s+glycoprotein\s+gp2|glycoprotein\s+2b|gs\s+protein|"
              r"glycosylated\s+membrane\s+protein\s+gp2|envelope\s+glycoprotein\s+gp2a?|"
              r"glycoprotein\s+gp2a|envelope\s+protein\s+gp2a|glycosylated\s+protein\s+gp2a?)\s*$"),
    ("GP3",   r"^\s*(gp3(\s+(envelope\s+)?(glyco)?protein)?|glyco(sylated|protein)?\s*(membrane\s+)?protein\s*3|"
              r"glycoprotein\s+3\s+gp3|orf3(\s+protein)?|orf3\s+envelope\s+protein|"
              r"envelope\s+glycoprotein\s+3|envelope\s+protein\s+gp3|glycoprotein\s+gp3|"
              r"glycosylated\s+membrane\s+protein\s+gp3|envelope\s+glycoprotein\s+gp3|"
              r"glycosylated\s+protein\s+gp3)\s*$"),
    ("GP4",   r"^\s*(gp4(\s+(envelope\s+)?(glyco)?protein)?|glyco(sylated|protein)?\s*(membrane\s+)?protein\s*4|"
              r"glycoprotein\s+4\s+gp4|glycoprotein\s+4\s+\(gp4\)|orf4(\s+protein)?|orf4\s+envelope\s+protein|"
              r"envelope\s+glycoprotein\s+4|envelope\s+protein\s+gp4|glycoprotein\s+gp4|"
              r"glycosylated\s+membrane\s+protein\s+gp4|envelope\s+glycoprotein\s+gp4|"
              r"glycosylated\s+protein\s+gp4)\s*$"),
    ("ORF5A", r"^\s*(orf\s?5a(\s+protein)?|gp5a(\s+protein)?|5a|"
              r"unglycosylated\s+protein\s+5a)\s*$"),
    ("GP5",   r"^\s*(gp5(\s+(envelope|surface|glyco)?\s*protein)?|gp5\s+glycosylated\s+envelope\s+protein|"
              r"glyco(sylated|protein)?\s*(membrane\s+)?protein\s*5|glycoprotein\s+5\s+gp5|"
              r"orf5(\s+protein)?|orf5\s+envelope\s+protein|orf7|envelope\s+glycoprotein\s+5|"
              r"envelope\s+protein\s+gp5|glycoprotein\s+gp5|major\s+envelope\s+glycoprotein\s+5|"
              r"glycosylated\s+envelope\s+protein\s+gp5|glycosylated\s+membrane\s+protein\s+gp5|"
              r"gl\s+protein|major\s+structural\s+glycoprotein\s+gp5|large\s+glycoprotein|"
              r"envelope\s+glycoprotein\s+gp5|major\s+envelope\s+glycoprotein(\s+gp5)?|"
              r"glycosylated\s+protein\s+gp5|major\s+envelope\s+(membrane\s+)?protein|"
              r"primary\s+envelope\s+glycoprotein)\s*$"),
    #  `glycosylated membrane protein` WITH NO NUMBER IS NOT M, and this rule
    #  used to take it via an optional `(un)?`. In arterivirus nomenclature the
    #  glycosylated membrane proteins are GP2-GP5 and M is the UNglycosylated
    #  one, so the prefix is the whole meaning of the string, not a modifier.
    #
    #  Measured 19 Nov 2026. Every NUMBERED form routes to a glycoprotein --
    #  `glycosylated membrane protein 2` x30 to GP2, `3` x76 to GP3, `4` x66 to
    #  GP4, `5` x80 to GP5 -- while `unglycosylated membrane protein` x78 and
    #  `... protein m` x13 are real M. The 36 deposits of the BARE form landed
    #  in M, and 35 of them are 178 aa and >=90% identical to GP4.
    #
    #  178 aa sits inside M's literature window (155-190) AND inside GP4's
    #  (80-210), so no length test could separate them, and at 2.1% of a 1,686
    #  member bin no fraction-based audit looked at them either. They clustered
    #  on their own into M.4 -- a profile built purely from GP4 -- which then
    #  called GP4 as M in 268 of 612 genomes while the real M went uncalled.
    #
    #  The bare string is genuinely ambiguous between GP2, GP3, GP4 and GP5, so
    #  it is NOT re-pointed at GP4 here. It falls to UNASSIGNED and the homology
    #  rescue arbitrates each deposit on sequence, which is the only evidence
    #  that can tell these four apart.
    ("M",     r"^\s*(m(\s+(membrane|matrix)\s+protein)?|m\s+protein|unglycosylated\s+membran[e]?\s+protein(\s+m)?|"
              r"membrane\s+(envelope\s+)?protein(\s+m|\s+gp6)?|matrix\s+protein(\s+m)?|"
              r"orf6(\s+(protein|membrane\s+protein))?|gp6(\s+protein)?|orf8|"
              r"non-?glycosylated\s+(envelope|membrane)\s+protein(\s+m)?|"
              r"unglycosylated\s+membrane\s+protein\s+m|envelope\s+protein\s+m)\s*$"),
    ("N",     r"^\s*(n(\s+protein)?|n\s+nucleocapsid\s+protein|nucleocapsid(\s+protein)?(\s+n)?|"
              r"nucleoprotein|np|orf7(\s+protein)?|orf7\s+nucleocapsid\s+protein|orf9|"
              r"gp7(\s+protein)?|capsid\s+protein|n\s+protein\s+\(nucleocapsid\)|"
              r"nucleocapsid\s+prtein\s+n)\s*$"),
]

#  A string that names one protein on a product whose LENGTH names another is
#  dropped here with its reason, never re-binned.  See the naming trap above.
LENGTH_ROUTE = {
    "ORF1A": (1700, 2700, "ORF1a-named deposit outside pp1a length: an ORF1b or pp1ab product "
                          "under a polyprotein string; both are transcript_edit, not PSSM features"),
    "E":     (45, 95,  "'E'/'ORF2b'-named deposit outside E length: the GP5 or GP2 product under "
                       "the other ORF2 convention (PMID 10400725 vs 11504553; PMID 7831770)"),
    "GP2":   (150, 320, "GP2-named deposit outside GP2 length: the E product under the EAV ORF2 "
                        "convention (PMID 10400725)"),
    "GP3":   (130, 280, "GP3-named deposit outside GP3 length"),
    "GP4":   (80, 210,  "GP4-named deposit outside GP4 length"),
    "ORF5A": (40, 80,   "ORF5a-named deposit outside ORF5a length"),
    "GP5":   (160, 285, "GP5/ORF5-named deposit outside GP5 length: GP3 under the renumbered "
                        "simian scheme, or an E deposit under the early PRRSV-1 'E' usage"),
    "M":     (155, 190, "M/ORF6-named deposit outside M length: GP4 under the renumbered scheme"),
    "N":     (100, 140, "N/ORF7-named deposit outside N length: GP5 under the renumbered "
                        "simian scheme, where ORF7 is GP5 and N is ORF9"),
    "NSP1":  (230, 290, "nsp1-named deposit outside the EAV nsp1 length (260 aa, PMID 7769711)"),
    "NSP1A": (160, 200, "nsp1alpha outside 180 aa (PMID 20006994)"),
    "NSP1B": (180, 230, "nsp1beta outside 203 aa (PMID 20006994)"),
    "NSP2":  (560, 1300, "nsp2-named deposit outside the arterivirus nsp2 range "
                         "(EAV 571, PRRSV-1 1078, PRRSV-2 1196): a partial nsp2 or an nsp2TF/nsp2N "
                         "frameshift product, which are transcript_edit features"),
    "NSP3":  (200, 350, "nsp3 outside 230-233 aa"),
    "NSP4":  (180, 230, "nsp4 outside 203-204 aa"),
    "NSP5":  (150, 190, "nsp5 outside 162-170 aa"),
    "NSP6":  (10, 30,   "nsp6 outside 16-22 aa"),
    "NSP7A": (130, 170, "nsp7alpha outside 123-149 aa"),
    "NSP7B": (90, 135,  "nsp7beta outside 102-120 aa"),
    "NSP8":  (35, 60,   "nsp8 outside 45-50 aa"),
    "NSP9":  (640, 720, "RdRp-named deposit outside nsp9 length: the ORF1b product, which has no "
                        "start codon of its own and is called through pp1ab (transcript_edit)"),
    "NSP10": (420, 490, "nsp10 outside 441-467 aa"),
    "NSP11": (200, 240, "nsp11 outside 219-224 aa"),
    "NSP12": (110, 170, "nsp12 outside 118-153 aa"),
}

NOT_MODELLED = [
    (r"^\s*(polyprotein\s*(pp)?\s*1ab|orf\s?1ab(\s+polyprotein|\s+protein)?|pp1ab|"
     r"1ab\s+replicase(\s+protein)?|replicase\s+polyprotein\s+1ab|replicase\s+orf1ab\s+polyprotein|"
     r"polyprotein\s+pp1ab|polyprotein\s+orf1ab)\s*$",
     "pp1ab: -1 ribosomal frameshift product at the U_UUA_AAC heptamer; called from "
     "nucleotide references as special:transcript_edit, no PSSM"),
    (r"^\s*(polyprotein\s*1b|orf\s?1b(\s+polyprotein)?|pp1b|1b|replicase\s+polyprotein\s+1b|"
     r"polyprotein\s+orf1b|nonstructural\s+protein\s+orf1b|putative\s+(protein\s+)?1b|"
     r"putative\s+1b\s+protein|replicase\s+protein\s+1b|rna\s+dependent\s+rna\s+polymerase|"
     r"rna-dependent\s+rna\s+polymerase/helicase|rna\s+polymerase|dna-dependent\s+rna\s+polymerase)\s*$",
     "ORF1b has no initiation codon of its own; its product exists only inside pp1ab, "
     "which is a transcript_edit feature"),
    (r"^\s*(polyprotein\s*1a(tf|n)|1a(tf|n)\s+polyprotein|nsp2(tf|n)(\s+protein)?|"
     r"tf\s+protein|(putative\s+)?(protein\s+)?nsp2f(\s+protein)?|"
     r"orf1a(tf|n)\s+polyprotein)\s*$",
     "nsp2TF / nsp2N: -2 and -1 ribosomal frameshift products of the nsp2 region "
     "(PMID 23043113); transcript_edit features, no PSSM"),
    (r"^\s*(orf1a/b\s+polyprotein|ppab|nonstructural\s+(poly)?protein[s]?\s+orf1b|"
     r"nonstructural\s+polyprotein\s+1b|proteinase\s+1b|zinc-finger\s+protein)\s*$",
     "an ORF1b or pp1ab product: ORF1b has no initiation codon of its own and pp1ab is a "
     "transcript_edit feature. The 'zinc-finger protein' string is used for a 1,460 aa ORF1b "
     "deposit, not for a zinc-finger domain"),
    (r"^\s*nsp7\s*$",
     "uncleaved nsp7 precursor (225-259 aa = nsp7alpha + nsp7beta): the module calls the two "
     "end products, whose internal site Glu1575-Ala1576 is conserved and lethal when abolished "
     "(PMID 17098961)"),
    (r"^\s*undefined\s+cds\s*$", "no product claim at all"),
    (r"^\s*unknown\s*$", "no product claim at all"),
    (r"^\s*hydrophobic\s+domain\s*$",
     "one EAV deposit uses this string for BOTH nsp3 (233 aa) and nsp5 (162 aa); "
     "a string that names two features names neither"),
]

#  Junction chemistry, from the MEASURED cleavage map (EAV: PMID 9371590,
#  9971783, 17098961, 7769711; PRRSV: PMID 20006994, 19646449, 20696193).
#  Every arterivirus nsp4 site is Glu-down, so each interior subunit ENDS in
#  Glu, and each one BEGINS with the P1' residue the papers report. Applied
#  here as a decontamination filter, after measuring the collections against it
#  (see BUILD_NOTES section on nsp2): 968 deposits named nsp2 hold 591 that
#  begin with the Edman-measured Ala of `AGKRARKAR` (PMID 20006994) and 279
#  that begin at an arbitrary internal residue -- the hypervariable-region
#  amplicons PRRSV epidemiology sequences, which are internal fragments of
#  nsp2 and not the deletion variants their lengths suggest.
#
#  NSP9 is NOT in this table as a PSSM feature at all. The RdRp begins at the
#  nsp7beta/nsp8 site inside ORF1a and runs through the -1 frameshift into
#  ORF1b (EAV pp1ab 1678-2370, PMID 9971783; there is no nsp8/nsp9 cleavage,
#  PMID 30353315), so it is not collinear with the genome -- the same shape as
#  coronavirus NSP12, which this repository already ships as
#  `special: transcript_edit`. Its collection is kept only to build the
#  nucleotide reference set.
JUNCTION = {
    "NSP1":  ("M",  "G"),
    "NSP1A": ("M",  "MH"),
    "NSP1B": ("AS", "G"),
    "NSP2":  ("AG", "GI"),
    "NSP3":  ("G",  "E"),
    "NSP4":  ("G",  "E"),
    "NSP5":  ("GS", "E"),
    "NSP6":  ("G",  "E"),
    "NSP7A": ("S",  "E"),
    "NSP7B": ("AN", "E"),
    "NSP8":  ("AG", None),
    "NSP9":  ("AG", "E"),
    "NSP10": ("GS", "EQ"),
    "NSP11": ("GS", "E"),
    "NSP12": ("G",  None),
    "ORF1A": ("M",  None),
}

WINDOW = {k: (v[0], v[1]) for k, v in LENGTH_ROUTE.items()}

if __name__ == "__main__":
    run(here=HERE, T="Art", MODULE="Arteriviridae", RULES=RULES,
        NOT_MODELLED=NOT_MODELLED, LENGTH_ROUTE=LENGTH_ROUTE, WINDOW=WINDOW,
        JUNCTION=JUNCTION)
