#!/usr/bin/env python3
"""Macrobrachium rosenbergii nodavirus -- the crustacean nodaviruses.

Members: *Macrobrachium rosenbergii* nodavirus (MrNV, 81 genomes), covert
mortality nodavirus (CMNV, 69) and *Penaeus vannamei* nodavirus (PvNV, 4).

### The name

ICTV MSL41 (2025) ratifies only **Alphanodavirus** and **Betanodavirus**; these
three sit under "unclassified Nodaviridae" in NCBI taxonomy with no genus at
any rank. PMID 30346944 asks "A new genus within the Nodaviridae?" in its
title but proposes no ratified name, so *Gammanodavirus* would be an invented
name for a module key. The module is therefore named for its most-sampled
species, as `Orthopneumovirus_muris` was -- curator's rule, 22 Nov 2026.
Membership is MrNV + PvNV + CMNV, and the name understates it on purpose.

### Layout -- and why it is not the Betanodavirus layout

    RNA1   5'-[----------- protein A (RdRp) -----------]-3'
                                       [--- B2 ---]           overlapping frame
    RNA2   5'-[----------- capsid protein -----------]-3'
           ONE uncleaved chain, no gamma peptide.

**There is no B1 here.** It is not merely unannotated: no paper reports one and
none of NC_005094.1 (MrNV), NC_014978.1 (PvNV) or MT270124.1 (CMNV) annotates
one. That single missing gene is what separates this module from Betanodavirus,
which otherwise shares the two-segment / uncleaved-capsid organisation. A module
declares a protein set, so declaring B1 here would flag every genome in the
group for a gene it does not carry.

**And no subgenomic RNA3 has ever been demonstrated for any of the three.** The
alpha- and betanodavirus B2s are translated from sgRNA3; for MrNV, PvNV and
CMNV that has not been shown, so nothing in this module's annotation strings
claims it.

### What the annotation deliberately does NOT say about B2

RNAi-suppressor activity is established for alphanodavirus B2 (PMID 12016316)
and betanodavirus B2 (PMID 16352533), and is **not demonstrated for MrNV, PvNV
or CMNV**. Those ORFs are sequence predictions carried by PMIDs 23318596 and
32659308. The feature is annotated "B2 protein", never "B2 RNAi suppressor":
inheriting a function from a homologue in another genus is exactly the kind of
claim the literature-authority rule exists to stop.

### The capsid is NOT cleaved, and the one reported cut is an artefact

PMID 30346944: "Betanodaviruses and Mr_NV both have their C-termini at the
capsid exterior and do not encode a gamma-peptide"; PMID 15660790 finds the
virion capsid is "a single polypeptide of 43 kDa". PMID 34445426 did locate a
cut at Arg26 by Edman sequencing -- caused by an *E. coli* cysteine protease,
blocked by E-64, in recombinant protein. That is a heterologous-expression
artefact, not viral maturation, and no mat_peptide is built from it.

### The capsid window is genuinely bimodal, and is set from the literature

MrNV 371 aa (PMID 23318596 and PMID 30346944, independently) against CMNV
437 aa (PMID 32659308). PvNV 368 is GenBank-only. Two length classes 66
residues apart in one collection will trip a SPLIT finding in audit_bins.py;
that is species-level length variation in one protein, not two proteins in one
bin, and it is signed off rather than split.

### One mislabelled string, caught on length

Two deposits read `RDRP` at **137 aa**. The polymerase here is 1033-1045 aa,
and 137 is B2 length. The length window routes them out rather than letting a
polymerase-named 137-mer into the A collection, where it would have become a
profile that calls B2's locus under the polymerase's name.
"""
import os, sys
sys.path.insert(0, os.environ["LOWVAN_KIT"] + "/skill/scripts")
from collection_engine import run
HERE = os.path.dirname(os.path.abspath(__file__))

RULES = [
    ("A",   r"^\s*(a\s+protein|protein\s+a|"
            r"rna[\s\-]?dep[ea]nd[ae]nt\s+rna[\s\-]?polymerase|"
            r"rna[\s\-]?dependant\s+rna[\s\-]?polymerase|"
            r"rna[\s\-]?directed\s+rna\s+polymerase|rdrp|rdrp\s+polymerase|"
            r"rna\s+polymerase|polymerase|putative\s+polymerase|"
            r"putative\s+rna[\s\-]?dep[ea]nd[ae]nt\s+rna[\s\-]?polymerase|"
            r"putative\s+rna[\s\-]?dependant\s+rna[\s\-]?polymerase|"
            r"putative\s+rdrp|replicase|putative\s+replicase)\s*$"),
    ("B2",  r"^\s*(b2(\s+protein)?|protein\s+b2|b2-?like\s+protein|"
            r"b\s+protein|protein\s+b|putative\s+protein\s+b2|"
            r"non-?structural\s+protein\s+b2)\s*$"),
    ("CP",  r"^\s*((major\s+|minor\s+)?(coat|capsid|caspid)\s+protein|"
            r"(coat|capsid)\s+protein\s+precursor|(coat|capsid)|"
            r"cp(\s+protein)?|structural\s+protein|"
            r"putative\s+(coat|capsid)(\s+protein)?)\s*$"),
]

LENGTH_ROUTE = {
    "A":   (1000, 1080,
            "polymerase-named deposit outside the RdRp range (CMNV 1043, PMID "
            "32659308; MrNV 1045 and PvNV 1033, GenBank-only -- the MrNV paper "
            "PMID 23318596 reports only a 1034 aa PARTIAL). Two deposits read "
            "`RDRP` at 137 aa, which is B2 length, and are routed out here"),
    "B2":  (125, 140,
            "B2-named deposit outside B2 length (MrNV 133, PMID 23318596; CMNV "
            "132, PMID 32659308; PvNV 131, GenBank-only). TIGHTENED 120-180 "
            "-> 125-140 on 22 Nov 2026 to exclude fig|222557.10.peg.1, a "
            "174 aa deposit labelled `B2 protein` under the MrNV taxon id "
            "that is NOT MrNV B2 and NOT MrNV capsid: it hits neither "
            "collection at e<1e-3, and its arginine-rich N-terminus is a "
            "capsid RNA-binding arm. It is EXTRA SMALL VIRUS capsid CP17 -- "
            "174 aa on an 872 nt contig (AY655724) ending at nt 610, against "
            "XSV NC_043492.1 at 872 nt with CP17 174 aa at CDS 89..613. XSV "
            "is a SATELLITE virus in Sarthroviridae, genus Macronovirus "
            "(PMIDs 14993649, 30507372), so it is outside this family "
            "entirely and would need a module of its own"),
    "CP":  (355, 450,
            "capsid-named deposit outside capsid length. BIMODAL BY SPECIES: "
            "MrNV 371 (PMIDs 23318596, 30346944), PvNV 368 (GenBank-only), "
            "CMNV 437 (PMID 32659308). The window spans both classes because "
            "both are real; it is not a median"),
}

NOT_MODELLED = [
    (r"^\s*(hypothetical\s+protein(\s+[12])?|unknown(\s+protein)?|"
     r"unnamed\s+protein\s+product|viral\s+protein|orf\s*[0-9]*)\s*$",
     "no product claim at all; left to the homology rescue"),
    (r"^\s*(b1(\s+protein)?|protein\s+b1)\s*$",
     "no B1 is documented in MrNV, PvNV or CMNV, and none is annotated in "
     "NC_005094.1, NC_014978.1 or MT270124.1. A B1-named deposit here is a "
     "label carried over from the alpha-/betanodavirus literature, not a "
     "protein this group is known to encode"),
]

JUNCTION = {
    "A":  ("M", None),
    "B2": ("M", None),
    "CP": ("M", None),
}

WINDOW = {k: (v[0], v[1]) for k, v in LENGTH_ROUTE.items()}

if __name__ == "__main__":
    run(here=HERE, T="N", MODULE="Macrobrachium_rosenbergii_nodavirus",
        RULES=RULES, NOT_MODELLED=NOT_MODELLED, LENGTH_ROUTE=LENGTH_ROUTE,
        WINDOW=WINDOW, JUNCTION=JUNCTION)
