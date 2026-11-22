#!/usr/bin/env python3
"""Betanodavirus -- the fish nervous necrosis viruses (SJNNV, RGNNV, BFNNV,
TPNNV, and the AHNV/GGNNV isolates inside them).

    RNA1   5'-[------------ protein A (RdRp) ------------]-3'
                                        [--- B1 ---]            same frame as A
                                        [--- B2 ---]            +1 frame vs A
           B1 and B2 come from subgenomic RNA3, which begins at RNA1 nt 2730
           and is made by premature termination, NOT by splicing (PMID
           16186236, PMID 41146288).

    RNA2   5'-[------------ coat protein (CP) ------------]-3'
           ONE uncleaved chain. No alpha/beta/gamma.

### The capsid is NOT cleaved, and this is the whole reason the module exists

PMID 9672600 measured it directly: the betanodavirus RNA2 ORF "carries the
catalytic D-75 residue but lacks the site for autocatalytic proteolysis, the
process yielding the two capsid proteins of insect nodaviruses". The same paper
excluded alternative initiation, a -1 frameshift and precursor degradation as
explanations for the 43/45 kDa doublet. The 338-residue intact-CP crystal
structure (PMID 26491970) and T=3 VLPs assembled from a single uncleaved
polypeptide corroborate it, and the 30-year review (PMID 32050492) states it
in terms: "This mechanism was not observed in betanodavirus".

So **no BETA and no GAMMA are declared here.** The N-ARM (1-33), S domain
(52-213) and P domain (221-338) are structural regions of one chain, not
mature peptides, and are not features.

### Windows come from the literature because the deposits cannot carry them

This genus is dominated by typing amplicons: of 391 `coat protein` deposits the
median length is **181 aa against a true 338-340**, and `RNA-dependent RNA
polymerase` has a median of 300 against a true 981-983. A window anchored on
any deposited median would delete every full-length protein in the genus and
keep the fragments -- the same failure Birnaviridae's VP2 hypervariable-region
amplicons produce. Every window below is the literature's, with its PMID in
length_refs.json.

### B1 has no published start codon anywhere

PMID 19136133 cloned "the B1 gene (336 nt)" from RGNNV and showed the protein
is real; PMID 36560822 states "its open reading frame matches the 3' termini of
the RdRp reading frame". But **no paper gives a start coordinate and no RefSeq
annotates a B1 CDS** -- NC_003448, NC_003449, NC_008040 and NC_008041 carry
only protein A, "protein B" (= B2) and the coat protein. The 30-year review
even records that "there has been some controversy about the existence of B1".

BV-BRC therefore holds **2 B1 deposits in 342 genomes**. The collection is
seeded from those and expanded by extraction from the genomes, gated on the
property the literature does fix: **B1 is in frame with protein A and shares
its stop, so an accepted B1 must be an exact C-terminal substring of protein A
from its own genome.** 124 of 124 extractions that had a deposited protein A to
test against passed; the extracted lengths are 110-112 aa, against the 111 aa
(336 nt) of PMID 19136133.

B2 is likewise thin -- 26 deposits -- and is extracted the same way, but it
cannot use the substring gate because it is in the +1 frame. It is gated on
length instead, and came out at exactly 75 aa in 126 of 127 genomes, which is
the length PMID 16186236 and PMID 15109133 report for SJNNV, RGNNV and AHNV.
"""
import os, sys
sys.path.insert(0, os.environ["LOWVAN_KIT"] + "/skill/scripts")
from collection_engine import run
HERE = os.path.dirname(os.path.abspath(__file__))

RULES = [
    #  --- RNA1 -------------------------------------------------------------
    ("A",     r"^\s*(protein\s+a|protein\s+a\s*\(rdrp\)|"
              r"rna[\s\-]?dep[ea]nd[ae]nt\s+rna[\s\-]?polymerase|"
              r"rna[\s\-]?depedent\s+rna[\s\-]?polymerase|"
              r"rna[\s\-]?dependant\s+rna[\s\-]?polymerase|"
              r"rna[\s\-]?dependent[\s\-]rna[\s\-]?polymerase|"
              r"rna[\s\-]?directed\s+rna\s+polymerase|rdrp|rdrp\s+polymerase|"
              r"rna\s+polymerase|polymerase|putative\s+polymerase|"
              r"putative\s+rna[\s\-]?dep[ea]nd[ae]nt\s+rna[\s\-]?polymerase|"
              r"putative\s+rdrp|replicase|putative\s+replicase)\s*$"),
    ("B1",    r"^\s*(protein\s+b1|b1\s+protein|b1|"
              r"non-?structural\s+protein\s+b1)\s*$"),
    ("B2",    r"^\s*(protein\s+b2|b2\s+protein|b2|putative\s+protein\s+b2|"
              r"protein\s+b|non-?structural\s+protein\s+b2)\s*$"),

    #  --- RNA2 -------------------------------------------------------------
    ("CP",    r"^\s*((major\s+|minor\s+)?(coat|capsid)\s+protein|"
              r"(coat|capsid)\s+protein\s+precursor|"
              r"(coat|capsid)\s+precursor\s+protein|"
              r"cp(\s+protein)?|(coat|capsid)|"
              r"putative\s+(coat|capsid)(\s+protein)?|"
              r"structural\s+protein|caspid\s+protein)\s*$"),
]

LENGTH_ROUTE = {
    "A":     (975, 990,
              "protein-A-named deposit outside the RdRp range (SJNNV 983 PMID "
              "10580064; RGNNV/GGNNV 982 PMIDs 22270757, 11172107; AHNV 981 "
              "PMID 15109133). THE DEPOSITED MEDIAN IS 300 aa -- most "
              "polymerase-named records in this genus are partial typing "
              "fragments, so anything below this window is a fragment, not a "
              "short lineage. TIGHTENED 960-1000 -> 975-990 on 22 Nov 2026: "
              "fig|43763.22.peg.1 measured 963 aa on a 2968 nt contig, about "
              "130 nt short of a complete betanodavirus RNA1, so it is a "
              "contig-truncated deposit and not a short lineage"),
    "B1":    (100, 125,
              "B1-named deposit outside B1 length (RGNNV 111 aa from a 336 nt "
              "gene, PMID 19136133; SJNNV 112). B1 is the C-terminal fragment "
              "of protein A, so a deposit near A's length is A"),
    "B2":    (68, 85,
              "B2-named deposit outside B2 length. 75 aa exactly in SJNNV, "
              "RGNNV, AHNV, BFNNV and TPNNV (PMIDs 16186236, 15109133, "
              "36560822)"),
    "CP":    (338, 342,
              "capsid-named deposit outside coat-protein length (SJNNV and "
              "TPNNV 340, PMID 9049363; RGNNV/GGNNV/BFNNV 338, PMIDs 22270757, "
              "11172107). THE DEPOSITED MEDIAN IS 181 aa: the majority of "
              "`coat protein` records here are partial sequences deposited for "
              "genotyping, and no window drawn from them could be right. "
              "TIGHTENED 325-350 -> 338-342 on 22 Nov 2026: the two members "
              "that fell below the literature floor, fig|43762.3.peg.1 at "
              "336 aa and fig|204929.15.peg.1 at 327 aa, sit on 1009 nt and "
              "995 nt contigs, where the coding sequence fills the entire "
              "record -- so the GENOME is what is truncated, not the protein"),
}

NOT_MODELLED = [
    (r"^\s*(hypothetical\s+protein(\s+[12])?|unknown(\s+protein)?|"
     r"unnamed\s+protein\s+product|viral\s+protein|orf\s*[0-9]*|"
     r"protein\s+[cd]|core\s+protein)\s*$",
     "no product claim this module can bind. `protein C` (106 aa) and "
     "`protein D` (87 aa) are single deposits using positional letters that "
     "name no homology group -- the literature describes only A, B1, B2 and "
     "CP for this genus (PMIDs 9672600, 9049363, 11172107), so these fall to "
     "UNASSIGNED and the homology rescue arbitrates"),
]

JUNCTION = {
    "A":  ("M", None),
    "B1": ("M", None),
    "B2": ("M", None),
    "CP": ("M", None),
}

WINDOW = {k: (v[0], v[1]) for k, v in LENGTH_ROUTE.items()}

if __name__ == "__main__":
    run(here=HERE, T="N", MODULE="Betanodavirus", RULES=RULES,
        NOT_MODELLED=NOT_MODELLED, LENGTH_ROUTE=LENGTH_ROUTE, WINDOW=WINDOW,
        JUNCTION=JUNCTION)
