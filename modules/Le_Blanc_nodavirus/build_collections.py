#!/usr/bin/env python3
"""Le Blanc nodavirus -- the Caenorhabditis-infecting nodaviruses.

Members: Le Blanc nodavirus (49 genomes), Santeuil nodavirus (44), Orsay virus
(4) and Melnik nodavirus (4).

### The name

No ICTV genus covers these -- MSL41 ratifies only Alphanodavirus and
Betanodavirus, and NCBI files all four under "unclassified Nodaviridae". The
module is named for the most-sampled species, as `Orthopneumovirus_muris` was
(curator's rule, 22 Nov 2026). The margin is thin and worth recording: Le Blanc
49, Santeuil 44, Orsay 4, Melnik 4. **Orsay virus is the best-studied member by
a wide margin and carries nearly all the literature cited below**, but it has
four genomes in BV-BRC, so the rule names the module for Le Blanc.

### Layout -- three products from one 2.3 kb region on RNA2

    RNA1   5'-[-------- protein A (RdRp) --------]-3'
           ONE CDS. No sgRNA3, no B1, no B2.

    RNA2   5'-[------ CP ------]
              [------ CP ----- -delta ----------]     -1 frameshift
                              [---- delta ------]     its own AUG

### Why this is its own module, and it is not a close call

Every other nodavirus module carries B2; this group carries none. PMID 21283608
looked and says so: "All known nodavirus B2 proteins overlap with the
C-terminus of the RNA-dependent RNA polymerase and are encoded in the +1 frame
relative to the polymerase. No open reading frame with these properties was
predicted in the 3' end of the RNA1 segment." PMID 23468484 then found no RNA
silencing suppressor activity for either Orsay protein tested. And this is the
only nodavirus group with a programmed ribosomal frameshift and a non-AUG
start. Nothing about its protein set matches any other module.

### CP starts at a non-AUG codon -- verified on the genome

PMID 24503084 Edman-sequenced the mature Orsay capsid N-terminus as **NKNNT**
and reports: "The asparagine corresponds to nucleotides 89-91, an AAC codon,
which is not traditionally thought to be a suitable translation initiation
site. In the complete RNA2 segment, there was no AUG codon 5' to the codon
corresponding to the N-terminal Asn." Checked directly on HM030971 (genome
977912.4) on 22 Nov 2026: **nt 89-91 are AAC**, in the context `GCAAAC`.

Consequences carried into the JSON, not just the notes:

- `upstream_ext` is **0** for CP and CP_DELTA. The annotator's
  `scan_to_met_start` walks left in 3-nt steps until it meets a Met or a stop,
  and here there is no Met upstream to find -- the scan would run to an
  arbitrary codon or off the contig. A terminus with no start codon must not be
  scanned for one.
- `check_met_start.py` will report CP and CP_DELTA. That is the documented
  non-AUG start, not a defect.
- The JUNCTION table below expects **N**, not M, at the start of both.

The mechanism is NOT settled and is not claimed here: the 2014 paper floats an
upstream UUG at nt 56-58 and then argues against it, preferring initiation
directly at the NKNNT site "via some unknown mechanism".

### CP-delta is a -1 ribosomal frameshift, so it gets no PSSM

Slippery heptamer **G_UCA_AAC at RNA2 nt 1255-1261**, immediately 5' of the CP
stop, with a stimulatory stem-loop ~8 nt downstream and **no pseudoknot**
(PMID 24503084). The transframe peptide was sequenced by MS/MS, and CP-delta
co-purifies with virions as a head fiber (PMID 32817218, PMID 28241071).
Verified on HM030971, 22 Nov 2026: nt 1255-1261 are **GTCAAAC** and the CP stop
is the TAG at nt 1262-1264, both exactly as published.

The RefSeq `join(89..1261,1261..2385)` is therefore a **frameshift, not a
splice**. CP_DELTA is declared `special: transcript_edit` and built from a
nucleotide reference set under `Transcript-Editing/`; it has no profile, so its
collection here exists only to source those references and to measure lengths.

### Free delta is a separate initiation, and it is real

Orsay delta has its own AUG at RNA2 nt 1345 (confirmed on the genome: ATG) and
is 346 aa, equal to CP-delta residues 420-765. PMID 29743360 settled that it
exists as a free protein: an ATG->CTG mutant that still makes normal CP-delta
but no free delta has a viral-egress defect. **It is absent from RefSeq
NC_028098.1** and present only in HM030971.2, which is why the deposits are
thin.

### Two annotation traps this file exists to defuse

Both are RefSeq's, and both would put a protein in the wrong collection:

1. **Santeuil NC_015070.1 labels the FUSION `protein delta`** and gives the
   capsid `/note="protein A"`. So the string `protein delta` covers both the
   765 aa fusion and the ~346 aa free protein, and the string `protein A` is
   used on RNA2 for a capsid. Length routing separates them -- DELTA at
   320-410 and CP_DELTA at 740-840 -- and `protein A` is NOT a rule here.
2. The deposits span `alpha-delta fusion protein` and `capsid-delta fusion
   protein` for the same thing, because submitters disagree about whether to
   call the capsid alpha. Both map to CP_DELTA.
"""
import os, sys
sys.path.insert(0, os.environ["LOWVAN_KIT"] + "/skill/scripts")
from collection_engine import run
HERE = os.path.dirname(os.path.abspath(__file__))

RULES = [
    ("A",        r"^\s*(rna[\s\-]?dep[ea]nd[ae]nt\s+rna[\s\-]?polymerase|"
                 r"rna[\s\-]?dependant\s+rna[\s\-]?polymerase|"
                 r"rna[\s\-]?directed\s+rna\s+polymerase|rdrp|"
                 r"rna\s+polymerase|polymerase|putative\s+polymerase|"
                 r"putative\s+rna[\s\-]?dep[ea]nd[ae]nt\s+rna[\s\-]?polymerase|"
                 r"replicase|putative\s+replicase)\s*$"),
    ("CP_DELTA", r"^\s*((alpha|capsid|coat|cp)[\s\-]delta\s+fusion\s+protein|"
                 r"(alpha|capsid|coat)-delta\s+fusion|"
                 r"capsid-?delta|fusion\s+protein)\s*$"),
    ("CP",       r"^\s*((major\s+)?(coat|capsid)\s+protein|(coat|capsid)|"
                 r"capsid\s+protein\s+alpha|protein\s+alpha|"
                 r"cp(\s+protein)?|putative\s+(coat|capsid)(\s+protein)?)\s*$"),
    ("DELTA",    r"^\s*(protein\s+delta|delta\s+protein|delta|orf\s*delta)\s*$"),
]

LENGTH_ROUTE = {
    "A":        (930, 1000,
                 "polymerase-named deposit outside the RdRp range (Orsay 937, "
                 "Le Blanc 942, Santeuil 982; RefSeq CDS lengths, with the "
                 "genome papers PMIDs 21283608, 23043172, 24503084 describing "
                 "the single ORF). Melnik is >=745 and partial at both ends, "
                 "so its RNA1 falls out here -- MK774657-MK774660 carry no "
                 "complete CDS and no RefSeq exists. FLOOR RAISED 880 -> 930 on 22 Nov 2026. "
                 "23 Le Blanc deposits measure exactly 906 aa and looked like "
                 "a FORM by the multiplicity rule, but they are 99.89% "
                 "identical to the 942 aa protein over all 906 residues, "
                 "aligned 1-906 to 1-906, so the 942 carries 36 MORE residues "
                 "at the C-terminus. Translating the ORF on the genome "
                 "settles it: on 1241918.3 (JQ943579, 3473 nt) it runs 942 aa "
                 "to a real stop, while on 1241918.10 (MN130974, 2976 nt) the "
                 "match ends at the contig's LAST BASE with no stop at all. "
                 "The 906 class is one 3'-truncated submission batch, not 23 "
                 "independent observations of a short form"),
    "CP":       (360, 400,
                 "capsid-named deposit outside capsid length (Orsay 391, "
                 "Santeuil 389, Le Blanc 388; PMID 24503084). A `protein "
                 "delta`-named deposit at this length is the Santeuil RefSeq "
                 "trap and belongs to DELTA or CP_DELTA, not here"),
    "CP_DELTA": (740, 840,
                 "fusion-named deposit outside CP-delta length (Orsay 765, "
                 "Santeuil 765, Le Blanc 827; PMID 24503084). The -1 "
                 "frameshift product, NOT a splice: RefSeq's "
                 "join(89..1261,1261..2385) repeats nt 1261, which is the "
                 "slipped base"),
    "DELTA":    (320, 410,
                 "delta-named deposit outside FREE delta length (Orsay 346 "
                 "from the AUG at RNA2 nt 1345, Santeuil 332, Le Blanc 398; "
                 "PMIDs 21283608, 29743360). Santeuil's RefSeq uses the same "
                 "string for the 765 aa FUSION, which this window routes to "
                 "CP_DELTA instead"),
}

NOT_MODELLED = [
    (r"^\s*(hypothetical\s+protein(\s+[12])?|unknown(\s+protein)?|"
     r"unnamed\s+protein\s+product|viral\s+protein|orf\s*[0-9]*)\s*$",
     "no product claim at all; left to the homology rescue"),
    (r"^\s*(b1|b2)(\s+protein)?\s*$",
     "these viruses encode no B1 and no B2. PMID 21283608 searched RNA1's 3' "
     "end for a +1 ORF with B2's properties and found none; PMID 23468484 "
     "found no silencing-suppressor activity for either Orsay protein tested"),
    (r"^\s*protein\s+a\s*$",
     "Santeuil RefSeq NC_015070.1 puts `protein A` on an RNA2 CAPSID as a "
     "/note. The string is therefore ambiguous between the RdRp and the "
     "capsid in this group and binds nothing; the homology rescue arbitrates"),
]

#  CP and CP_DELTA begin at Asn -- the Edman-verified non-AUG start (PMID
#  24503084), confirmed as an AAC codon on HM030971. A and DELTA have ordinary
#  Met starts.
JUNCTION = {
    "A":        ("M", None),
    "CP":       ("N", None),
    "CP_DELTA": ("N", None),
    "DELTA":    ("M", None),
}

WINDOW = {k: (v[0], v[1]) for k, v in LENGTH_ROUTE.items()}

if __name__ == "__main__":
    run(here=HERE, T="N", MODULE="Le_Blanc_nodavirus", RULES=RULES,
        NOT_MODELLED=NOT_MODELLED, LENGTH_ROUTE=LENGTH_ROUTE, WINDOW=WINDOW,
        JUNCTION=JUNCTION)
