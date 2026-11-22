#!/usr/bin/env python3
"""Birnaviridae -- one module for the whole family.

Bisegmented dsRNA. Segment B carries one product, VP1 (RdRp). Segment A carries
a polyprotein and, in a separate overlapping or internal reading frame, VP5:

    segment A   5'-[VP5]--pVP2--(X)--VP4--VP3-3'
    segment B   5'-VP1-3'

### Why one module and not several

Partitioning is by genome organisation, not phylogeny, and the organisation is
the same in every genus. The only layout difference is the X insert between the
pVP2 C-terminal peptides and VP4, present in Blosnavirus and Telnavirus and
absent elsewhere -- 22 of this family's ~9,200 records. That is far too small to
carry a module of its own (the curator's floor is 10 genomes and these are 2),
so X is carried HERE as a lineage-specific feature with no copy_num.

The budget gate agrees, measured 18 Nov 2026 on all 330 merged A+B genomes:
lumped gives 22 clusters and 100% routing at the 25-reference budget; splitting
Avibirnavirus out gives 3 clusters/100% and Aquabirnavirus 1 cluster/100%.
Splitting buys nothing, so the family stays whole.

### The architecture that decides this file

**Everything on segment A except VP5 is a cleavage product of one polyprotein.**
Measured 19 Nov 2026, every VP2-class deposit in the dump aligns to polyprotein
residue 1 and differs only at its C-TERMINUS:

    Avibirnavirus   441 (n=61)   452 (n=46)   465 (n=7)   512 (n=5)
    Aquabirnavirus  442 (n=52)   504 (n=20)   512 (n=19)  538 (n=17)

and the IBDV / IPNV cleavage maps place those ends exactly:

    IBDV 1012 aa   VP2 1-441  pep46 442-487  pep7a 488-494  pep7b 495-501
                   pep11 502-512  VP4 513-755  VP3 756-1012     (PMID 11836417)
    IPNV  972 aa   VP2 1-442  p1    443-486   p2    487-495
                   p3    496-508  VP4 509-734  VP3 735-972      (PMID 10666235,
                                                                 PMID 17553791)

So `VP2` as an annotation STRING names three different proteins -- mature VP2,
the pVP2 precursor, and a hypervariable-region typing amplicon -- and the string
cannot separate them. LENGTH_ROUTE here DROPS an out-of-window deposit; it does
not move it to a neighbouring feature. Therefore the cleavage products are
SEEDED from the deposits that are unambiguous and then EXPANDED by projection
from literature-anchored references (project_matpeptides.py), rather than being
binned on a string that means whichever product the submitter had in mind.

### The amplicon hazard, which dominates this family by volume

87% of Birnaviridae records are sub-2 kb fragments and 83% of the VP2 bin is
IBDV hypervariable-region typing amplicon, most of it deposited under a bare
`VP2`. The consequence is that **no length window in this family may be
anchored on a deposited median**: the Avibirnavirus polyprotein median is
328 aa against a true 1012, and the Avibirnavirus VP1 median is 255 against a
true 878. Every window here comes from length_refs.json, with a PMID.

### Masses are not a homology group

Aquabirnavirus deposits `12.1 kDa`, `14.7 kDa`, `15.4 kDa`, `15-kDa`, `16 kDa`
and `25-kDa nonstructural protein`, and Avibirnavirus `17 kDa nonstructural
protein`, plus bare `p9`, `p13` and `p44`. Measured 19 Nov 2026 the 133 aa
cluster is VP5 and the 145 aa cluster is VP5, but they are NOT binned on the
mass -- there is deliberately no rule for any of these strings. They fall to
UNASSIGNED and the homology rescue arbitrates, which is the only treatment that
cannot invent a homology group out of a molecular weight.
"""
import os, sys
sys.path.insert(0, os.environ["LOWVAN_KIT"] + "/skill/scripts")
from collection_engine import run
HERE = os.path.dirname(os.path.abspath(__file__))

RULES = [
    #  --- segment A, polyprotein level -------------------------------------
    #  POLY is the projection parent; project_matpeptides.py reads POLY.fasta.
    ("POLY", r"^\s*(polyprotein(\s+precursor)?|structural\s+polyprotein|"
             r"segment\s+a\s+polyprotein|polyprotein\s+vp0|polyprotein\s+vp243|"
             r"vp2-?4-?3(\s+pol[hy]+protein)?|vp2-vp4-vp3|"
             r"structural\s+polyprotein\s+vp2,\s*vp4,\s*vp3|"
             r"polyprotein\s+\(vp2-vp4-vp3\)|"
             r"1[01][0-9](\.[0-9])?\s*-?\s*kda\s+polyprotein|"
             r"putative\s+polyprotein)\s*$"),

    #  --- segment A cleavage products (SEEDS; expanded by projection) ------
    ("VP2",   r"^\s*(vp2|vp2\s+protein|viral\s+protein\s+2|capsid\s+protein|"
              r"capsid\s+protein\s+vp2|major\s+capsid\s+protein(\s+vp2)?|major\s+capsid|"
              r"structural\s+protein\s+vp2|vp2\s+structural\s+protein|"
              r"structural\s+vp2\s+protein|structural\s+protein\s+\(vp2\)|"
              r"major\s+structural\s+protein(\s+vp2)?|mature\s+vp2|vp2a|"
              r"outer\s+capsid\s+protein(\s+vp2)?|outer\s+viral\s+capsid\s+protein|"
              r"virion\s+structural\s+protein\s+vp2|vp2\s+capsid\s+protein|"
              r"viral\s+capsid\s+protein\s+2|virus\s+capsid\s+protein|"
              r"structural\s+protein|structural\s+viral\s+protein|"
              r"glycoprotein\s+vp2|vp2\s+surface\s+protein|coat\s+protein|"
              r"host\s+protective\s+antigen|putative\s+surface\s+protein)\s*$"),
    ("PVP2",  r"^\s*(vp2\s+precursor|vp2\s+protein\s+precursor|prevp2|"
              r"structural\s+polyprotein\s+vp2\s+precusor|polyprotein\s+vp2)\s*$"),
    ("PEP46", r"^\s*(pep4[0-9](\s+(peptide|protein))?|vp2a\s+protein)\s*$"),
    ("PEP7A", r"^\s*(pep7a(\s+(peptide|protein))?|vp2b\s+protein)\s*$"),
    ("PEP7B", r"^\s*(pep7b(\s+(peptide|protein))?|vp2c\s+protein)\s*$"),
    ("PEP11", r"^\s*(pep1[0-9](\s+(peptide|protein))?|vp2d\s+protein)\s*$"),
    ("X",     r"^\s*(x(\s+protein)?|protein\s+x)\s*$"),
    ("VP4",   r"^\s*(vp4|vp4\s+protein|viral\s+protein\s+4|"
              r"structural\s+protein\s+vp4|putative\s+vp4)\s*$"),
    ("VP3",   r"^\s*(vp3|vp3\s+protein|viral\s+protein\s+3|capsid\s+protein\s+vp3|"
              r"structural\s+protein\s+vp3|minor\s+structural\s+protein\s+vp3|"
              r"submajor\s+capsid\s+protein\s+vp3|putative\s+vp3)\s*$"),

    #  --- segment A, separate ORF ------------------------------------------
    ("VP5",   r"^\s*(vp5|vp5\s+protein|viral\s+protein\s+5|protein\s+vp5|"
              r"nonstructural\s+protein\s+vp5|structural\s+protein\s+vp5|"
              r"pvp5|putative\s+vp5)\s*$"),

    #  --- segment B ---------------------------------------------------------
    ("VP1",   r"^\s*(vp1|vp1\s+protein|viral\s+protein\s+1|viral\s+protein\s+vp1|"
              r"rna[\s\-]?dep[ea]nd[ae]nt\s+rna[\s\-]?polymerase(\s+\(rdrp\))?(\s+vp1)?|"
              r"rna[\s\-]?directed\s+rna\s+polymerase(\s+\(ec\s+2\.7\.7\.48\))?|"
              r"rna[\s\-]?dependant\s+rna[\s\-]?polymerase|rna\s+depended\s+rna\s+polymerase|"
              r"double-stranded\s+rna-dependent\s+rna\s+polymerase\s+vp1|"
              r"vp1\s+rna-dependent\s+rna\s+polymerase|rdrp\s+\(vp1\)|"
              r"rna\s+polymerase|viral\s+polymerase|virus\s+polymerase|polymerase|"
              r"rna-dependent\s+polymerase|possible\s+polymerase|"
              r"nonstructural\s+protein\s+vp1|putative\s+vp1|"
              r"putative\s+rna-dependent\s+rna\s+polymerase(\s+vp1)?|"
              r"putative\s+replicase)\s*$"),
]

#  A deposit whose string names one product and whose LENGTH names another is
#  dropped here with its reason recorded, never silently re-binned. Every
#  window comes from length_refs.json and carries a PMID there.
LENGTH_ROUTE = {
    "POLY":  (940, 1160,
              "polyprotein-named deposit outside the family polyprotein range "
              "(IPNV 972, IBDV 1012, BSNV 1069, TV-1 1114): 83% of IBDV records "
              "are VP2 hypervariable-region amplicons deposited under a "
              "polyprotein string, which is why the deposited median is 328 aa"),
    "VP2":   (395, 475,
              "VP2-named deposit outside MATURE VP2 length (IBDV 441, IPNV 442, "
              "BSNV 416, TV-1 450; the 452 and 465 classes are documented IBDV "
              "maturation intermediates): either a hypervariable-region amplicon "
              "below the window, or the pVP2 PRECURSOR above it, which is a "
              "different feature -- see PVP2"),
    "PVP2":  (480, 560,
              "precursor-named deposit outside pVP2 length (IBDV 1-512 PMID "
              "11836417, IPNV 1-508 PMID 10666235)"),
    "PEP46": (36, 52,  "first pVP2 C-terminal peptide outside 41-46 aa"),
    "PEP7A": (6, 10,   "second pVP2 C-terminal peptide outside 7-9 aa"),
    "PEP7B": (6, 9,    "third pVP2 C-terminal peptide outside 7 aa"),
    "PEP11": (9, 16,   "last pVP2 C-terminal peptide outside 11-13 aa"),
    "X":     (58, 120, "X-named deposit outside X length (BSNV 71, TV-1 106)"),
    "VP4":   (165, 295,
              "VP4-named deposit outside protease length (IBDV 513-755 = 243, "
              "IPNV 509-734 = 226, BSNV 234, TV-1 212)"),
    "VP3":   (215, 310,
              "VP3-named deposit outside VP3 length (IBDV 257, IPNV 238, "
              "BSNV 278, TV-1 284)"),
    "VP5":   (92, 160, "VP5-named deposit outside VP5 length (IBDV 145, IPNV 133)"),
    "VP1":   (780, 1060,
              "polymerase-named deposit outside VP1 length (IBDV 878, IPNV 845, "
              "BSNV 867, TV-1 894). The Avibirnavirus VP1 median is 255 aa for "
              "the amplicon reason above; the window is the literature's, not "
              "the deposits'. UPPER BOUND RAISED 935 -> 1060, 19 Nov 2026: the "
              "four type species are all Avibirnavirus/Aquabirnavirus-sized, "
              "and a ceiling drawn from them excludes every Entomobirnavirus "
              "and Ronavirus VP1 in the dump. Shwed et al. 2002 (Virology "
              "296:241-250, PMID 12069523) sequenced Drosophila X virus "
              "segment B and reported 'the largest birnavirus VP1' at 977 aa, "
              "which is already 42 aa above the old ceiling; the Ronavirus "
              "(Rotifer birnavirus) deposit is 1045 aa. Measured before the "
              "change: 16 polymerase-named deposits fall in 936-1060 aa, and "
              "15 of them hit the VP1 collection at e-23 to e-168 (23-38% "
              "identity over 55-85% of the query) with ZERO hits to POLY at "
              "e<1e-3 -- so the window was excluding real VP1, not admitting "
              "polyprotein. The 16th, fig|2921225.4.peg.1 (Botrytis cinerea "
              "binarnavirus 6), has no detectable homology to ANY feature in "
              "this module. It is not excluded by a rule, because its string "
              "('RNA dependent RNA polymerase') is one that real VP1 deposits "
              "use and no genome-level exclusion exists in collection_engine; "
              "it is left for audit_bins.py to flag, and the fungal "
              "binarnavirus/letobirnavirus lineage is recorded as an "
              "out-of-scope gap in NOT_CALLABLE.md"),
}

NOT_MODELLED = [
    #  --- typing amplicons: the single biggest population in this family ----
    (r"^\s*((vp2|viral\s+protein\s+2|capsid\s+protein,?\s*vp2,?|polyprotein\s+vp2|"
     r"vp2\s+protein)?\s*,?\s*(hypervariable|variable)\s+(region|domain)|"
     r"hypervariable\s+region|"
     r"vp2\s+(protein\s+)?hypervariable\s+(region|domain)|"
     r"viral\s+protein\s+2\s+variable\s+region|vp2\s+variable\s+region|"
     r"capsid\s+protein,\s*vp2,\s*variable\s+domain|"
     r"polyprotein\s+vp2\s+hypervariable\s+region)\s*$",
     "IBDV VP2 hypervariable-region typing amplicon: an internal fragment of "
     "VP2 deposited for genotyping, not a product. 83% of the VP2 bin by volume"),

    #  --- epitope sub-fragments --------------------------------------------
    (r"^\s*vp[23]\s+epitope\s+(i{1,3}|iv|v|vi)\s*$",
     "a mapped epitope is an internal fragment of VP2 or VP3, not a product"),

    #  --- a deposit measured to be wrong ------------------------------------
    (r"^\s*nonstructural\s+protein\s+vp4\s*$",
     "the sole deposit under this string is 268 aa and aligns to IPNV "
     "polyprotein 454-721 (measured 19 Nov 2026), so it begins 55 residues "
     "BEFORE the 508/509 cleavage site (PMID 10666235) and carries the three "
     "pVP2 C-terminal peptides into the VP4 bin, while stopping 13 residues "
     "short of the real C-terminus. It falls INSIDE the VP4 length window, so "
     "the window cannot catch it and it is excluded by name"),

    #  --- no product claim at all -------------------------------------------
    (r"^\s*(undefined\s+cds|unknown(\s+protein)?|hypothetical\s+protein|"
     r"unnamed\s+protein\s+product|viral\s+protein|oh|"
     r"region\s+du\s+segment\s+a|(large|small)\s+orf|orf5)\s*$",
     "no product claim at all. `ORF5` is a positional label, not a homology "
     "group: it is used here for a 1090 aa polyprotein and an 856 aa VP1"),
    (r"^\s*unnamed\s+protein\s+product;.*$",
     "a free-text submitter note, not a product name"),
]

#  Junction chemistry, MEASURED on the deposited polyproteins at the literature
#  coordinates before being applied (19 Nov 2026). VP4 is a Ser/Lys serine
#  protease (catalytic Ser633/Lys674 in IPNV, PMID 10666235) cleaving
#  (Ser/Thr/Ala)-X-Ala | (Ala/Ser):
#
#    IPNV  p1/p2  KTNA|AGGR   p2/p3  HSMA|AGGR   p3/VP4 ESWA|SGGP
#          VP4/VP3 ASTA|SGMD
#    IBDV  pep46/7a EAQA|ASGT  7a/7b TARA|ASGK   7b/11  KARA|ASGR
#          pep11/VP4 LTLA|ADKG  VP4/VP3 LAMA|ASEF
#    BSNV  VP2d/X  IPLA|SSDE   X/VP4 RPQA|ADLP   VP4/VP3 YCGA|ADEE
#
#  P1 is Ala at 15 of the 16 sites tested and P1' is Ala or Ser at every VP4
#  site. The VP2 maturation cut is a DIFFERENT chemistry -- IBDV IAGA|FGFK,
#  IPNV TSKA|WGWR -- P1 Ala but P1' aromatic, consistent with it being the
#  assembly-dependent maturation cleavage (PMID 11836416) rather than a VP4 cut.
#  That is why PEP46 admits F and W as well as A at its start.
JUNCTION = {
    "POLY":  ("M",   None),
    "VP2":   ("M",   "A"),
    "PVP2":  ("M",   "A"),
    "PEP46": ("AFW", "A"),
    "PEP7A": ("A",   "A"),
    "PEP7B": ("A",   "A"),
    "PEP11": ("A",   "A"),
    "X":     ("AS",  "A"),
    "VP4":   ("AS",  "A"),
    "VP3":   ("AS",  None),
    "VP5":   ("M",   None),
    "VP1":   ("M",   None),
}

WINDOW = {k: (v[0], v[1]) for k, v in LENGTH_ROUTE.items()}

if __name__ == "__main__":
    run(here=HERE, T="Birna", MODULE="Birnaviridae", RULES=RULES,
        NOT_MODELLED=NOT_MODELLED, LENGTH_ROUTE=LENGTH_ROUTE, WINDOW=WINDOW,
        JUNCTION=JUNCTION)
