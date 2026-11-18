# Bornaviridae — build notes

Built 17–18 Nov 2026. **One module, `Bornaviridae`, eleven features, one segment.**

> **Revised 18 Nov 2026 after curator review.** Four changes: **L is now a
> `special: splice` feature** built from literature-anchored junctions instead
> of a PSSM CDS (the PSSM put up to 46 spurious residues on the orthobornavirus
> N-terminus, which the curator ruled never acceptable); **N′ (p38) is now
> modelled** as `N_PRIME`, separated from N by a measured bit cutoff; the
> question of a "G′" and of a G/L fusion was settled from the literature (there
> is no G′; the G/L fusion is mRNA-only); and the SignalP installation was
> checked and found to be missing its model weights, which is recorded in
> SKILL.md and leaves the glycoprotein signal peptide unmodelled for now.
Local build tree: `~/lowvan_scratch/borna/work`. Runtime test repo:
`~/lowvan_scratch/borna/VA_borna`.

Every literature fact and every verified PMID: **`notes/LITERATURE.md`**.
No PMID in this module was written from memory; each was resolved against
PubMed in this build and its title printed. All of them sit in
`PMID_claude_generated`; `PMID` is empty on every feature, as it must be until
a curator reads them.

## The dump

1,372 genomes, 3,424 features, 1,199 unique proteins, fetched with
`dump_api.py` in 19 s. `check_dump.py` **PASS**: every protein md5 has a
sequence and the three `uniq.*` files are line-aligned. 296 `id_md5` rows
carry a blank md5 (RNA and misc features) and are correctly absent.

Genus distribution of the genome records: Orthobornavirus 1,193,
Cultervirus 11, Carbovirus 4, Cartilovirus 1, unclassified 163.
**No `<taxon> sp.` record is present**, so the intake rule removed nothing.

## One module, and why — measured, not assumed

Gene **content** is identical in every genus: N, X, P, M, G, L. Gene **order**
is not — M and G are transposed in Carbovirus, Cultervirus and Cartilovirus
(PMID 42358456, PMID 29462190), and Cartilovirus inserts two extra ORFs
(PMID 38028148). Order does not enter a PSSM, so the partitioning question is
settled by routing, and routing says lump:

    greedy reference selection over the 218 near-complete genomes
    (>= 6,000 nt), annotator BLASTn settings (-reward 2 -penalty -3
    -word_size 11), bit >= 150:

        4 references   -> 100.0% of all 218, every genus
       16 references   -> 100.0% at a 4,000-bit floor

**16 references were shipped, not 4.** Four reach every genome at the routing
threshold, but they leave the divergent genera hanging on a single marginal
hit; selecting greedily at a 4,000-bit floor spreads one reference per species
group and costs nothing against the 25-reference budget. The shipped set is
the set that was measured.

Splitting was therefore not tried further: there is no routing gap to close,
and a Cartilovirus module would hold **one genome**.

## Collections

`build_collections.py` binds on the annotation string inside an explicit
length window per feature; `rescue_unassigned.py` then adopts by homology at
the standard 80% identity over 85% of the query.

    binned by rule        3,055 of 3,128 features = 97.7%
    rescued               27 of 35 remaining unique sequences
    left unbound          8 unique sequences (see collections/UNASSIGNED_DECISION.md)

| feature | members | window | why that window |
|---|---|---|---|
| N | 129 | 340–400 | BoDV-1 370 aa; Cult 356–383, Cart 365, Carb 370 |
| X | 134 | 80–170 | BoDV-1 87 aa; Cultervirus runs 98–156 |
| P | 124 | 160–240 | BoDV-1 201 aa; Carb 167–170, Cult 187–205 |
| M | 62 | 125–205 | BoDV-1 142 aa; one Orthobornavirus at 198 |
| G | 154 | 480–560 | BoDV-1 503 aa; observed 499–529 family-wide |
| L | 158 | 1680–1760 | spliced form 1708–1737 |
| G_N | 150 | derived | furin cut, see below |
| G_C | 150 | derived | furin cut, see below |
| VP1 | 1 | 380–460 | Cartilovirus only |
| VP2 | 1 | 190–260 | Cartilovirus only |

**Mass names are bound, but only inside the window.** This family's literature
names its proteins by mass — p40 = N, p10 = X, p24/p23 = P, p16/gp18 = M,
p57/gp94 = G — so those strings carry real information and are bound. The
window is what makes that safe: a `p24 ORF` deposit of 147 aa is not P and
does not enter the collection.

**The L window deliberately excludes the 1,608 aa form** (13 BoDV-1 deposits).
That is the L ORF read from the internal AUG without the splice, and it is an
exact C-terminal substring of the 1,711 aa protein — RefSeq NC_001607
annotates exactly it (CDS 3996..8822 inside a gene 3687..8822, and
3996 − 3687 = 309 nt = 103 codons, the measured offset). Letting it into the
collection would have put a 103-residue-truncated form into the same profile
set as the full protein.

### Four truncated P deposits, caught and dropped

`check_met_start.py` reported one P cluster whose three members begin with Gly.
Looking at the collection, **four P deposits (160–182 aa) begin 19–40 residues
inside the protein** — `TLRRERSGSP...`, `GSPRPRKVPR...`, `RPRKDPRNAL...`
against the full `MATGPSSLVDSLEDEEDPQ|TLRRERSGSP...` — all from one BoDV-1
strain series (`1714621.145`–`150`), the same submitter whose other deposits
appear in the unassigned pool as `p24 ORF` fragments. They passed the 160–240
window because they are long enough. Dropped
(`collections/NTERM_TRUNCATED_DROPPED.tsv`), P rebuilt from 124 members, and
the check is now clean: **every CDS profile begins at its protein's start.**

## The glycoprotein cleavage — the one proven cleavage in this family

Richt 1998 (**PMID 9557754**) shows furin cleaves the gp94 precursor and places
the site at **the arginine at position 249** of the 503 aa ORF-IV product.
Derived here by locating the single `R-X-[K/R]-R` match in the 40–62% window
and extending through the basic run:

    exactly one site        150 / 154 members   (4 ambiguous, left out of both products)
    P1  = Arg              150 / 150
    P1' = Asp              149 / 150   (Glu in 1)
    P4  = Arg              150 / 150
    BoDV-1 503 aa form     cuts after residue 249 -> 249 aa + 254 aa

The BoDV-1 context is `...RLKRRRR|DTQQI...`. The site is enzyme-determined, so
it re-locates itself in every genus instead of being a column number — which is
the whole reason to derive rather than to model the products from deposits
(there are none: nothing in BV-BRC is annotated as a mature bornavirus G).

`apply_signalp.py --motif` was not used: the motif is genus-variable
(`KRRRR`, `AKPRKKR`, `RRKR`, …) and a single regex would have cut the wrong
lineages, whereas "the one polybasic run in this window" generalises.
`derive_g_products.py` holds the derivation.

**The signal peptide is NOT modelled** and `G_N` therefore retains it. No
mature N-terminus has been determined for any bornavirus G, no SignalP licence
is available here, and RefSeq NC_001607 carries no `sig_peptide` or
`mat_peptide` to read one off. Shipping a predicted `_SP` would have been a
prediction wearing a measurement's clothes. This is in §07 of the audit.

### VP1 and VP2 could not be bound by their string, and were not

**Both Cartilovirus deposits carry the identical annotation `viral protein 1`.**
A string rule cannot tell them apart, and a U-number label is not a homology
group, so they were assigned from genome position and from the paper's own
description (PMID 38028148): in the single Cartilovirus genome the ORFs run
`3'-N-Vp1-Vp2-X/P-G-M-L-5'`, and the deposits sit at 1172..2419 (416 aa) and
2451..3110 (220 aa) in that order. Corroborated structurally: tblastn of the
416 aa deposit against its own genome hits the **G** locus at 31% identity and
the 220 aa deposit hits the **416 aa** locus at 42%, which is the duplication
history the paper reconstructs. They therefore appear in
`collections/UNASSIGNED_TRACKING.tsv` as two unbound `viral protein 1`
sequences; that row is the record of a decision taken by position, not a gap.

## Profiles

Built at the **documented defaults** — `-m 5 -mi 0.8 -mc 0.8` — with no
departure, for eight of the ten features. Cluster coverage at >= 5 members
and `-c 0.8` was measured first, which is what showed no departure was needed:

    mi      0.8   0.6   0.5   0.4   0.3
    N       89%   91%   91%   95%   98%
    X       82%   87%   89%   89%   89%
    P       85%   85%   85%   89%   89%
    M       79%   83%   83%   95%   95%
    G       87%   92%   92%   96%  100%
    L       87%   90%   96%  100%  100%

79–89% at the default, and the leftover pass covers the rest at the same
identity floor rather than lowering it. **57 profiles**: N 6, X 8, P 6, M 4,
G 8, G_N 8, G_C 7, L 8, VP1 1, VP2 1 — 26 first-pass, 29 from
`build_leftover_pssms.py` at a member floor of 2.

`rebuild_pssms.py`: **57 reproduce, 0 stale.**

### VP1 and VP2 are single-sequence profiles, and that is stated

The whole BV-BRC holding for Cartilovirus is **one genome** (little skate
bornavirus, BK063518), so `fasta-cluster-pssm-2.pl` produces no alignment and
the kit's `pssms_from_alignments.py` refuses a one-row MSA. The profiles were
built directly with `psiblast -in_msa` from that sequence, duplicated to meet
psiblast's two-row minimum; the alignment on disk is that two-row file, so the
profile regenerates exactly and `rebuild_pssms.py` reports it clean. They cover
100% of the Cartilovirus material that exists, which is one genome, and that
is the honest description of them.

**Their cutoffs are the measurement, not a default.** At the generated
cutoffs (60 and 30) VP1 fired on the *glycoprotein* locus of 14 panel genomes
— the expected failure for a profile built from a G paralogue at 28% identity.
Measured scores:

    VP1   on its own Cartilovirus genome   887 bits
          on Wuhan sharpbelly (Cultervirus) 109
          on BoDV-1 / parrot / python        41-60
    VP2   on its own genome                 448 bits
          on everything else               17-35

Both cutoffs were set to **200** — 1.8x above the worst foreign score and 2.2x
below the true one. The false calls went to zero and VP1/VP2 now call only
where they belong.

## QC

| check | result |
|---|---|
| `check_dump.py` | PASS |
| `audit_collection_lengths.py --require-refs` | **Unresolved: 0** |
| `audit_bins.py` | **Unresolved: 0** — every collection is one protein and its own protein |
| `qc_cross_feature.py` | no cross-feature contamination |
| `audit_contamination.py` | 1,063 sequences tested against 56 references, **clean** |
| `check_met_start.py` | every CDS profile begins at its protein's start |
| `rebuild_pssms.py` | 57 reproduce, **0 stale** |
| `install_module.py --check` | 1 installable, 0 problems |
| `qc_duplicate_calls.py` | EXACT none; the NESTED rows are declared architecture (below) |

### The nested calls are the genome's architecture, not a defect

- `G_N` and `G_C` inside `G` — the declared precursor/product relationship.
- `X` overlapping `P` — real: BoDV-1 X is 1223..1486 and P is 1272..1877, so
  the X ORF overlaps P's 5' end in a different frame. Both are called at
  **exactly** the RefSeq coordinates (below), so this is the virus's own
  arrangement.

## The held-out panel, scored against the submitters' own proteins

70 genomes drawn with `make_panel.py`, the 16 reference genomes **excluded by
accession**. Each call was compared to that genome's own deposited protein:

| feature | compared | identical | otherwise |
|---|---|---|---|
| M | 46 | **46** | — |
| P | 46 | **46** | — |
| X | 43 | **43** | — |
| N | 46 | **43** | 3 calls 13–17 aa long at the N-terminus |
| G | 33 | **21** | 12 calls 3–25 aa long at the N-terminus |
| G_C | 19 | n/a | 19 of 19 are the deposit minus its N-terminal 242–249 residues — i.e. exactly the furin product, by construction |
| L | 30 | **8** | see below |
| VP1 / VP2 | — | — | 0 calls outside Cartilovirus after the cutoff fix |

On the BoDV-1 reference genome the module reproduces **RefSeq NC_001607 base
for base**: N 54..1166, X 1223..1486, P 1272..1877, M 1893..2321,
G 2236..3747, and `GP_N` 2236..2982 = 249 codons — the Richt 1998 product to
the residue — with `GP_C` 2983..3747 running to the precursor's own stop.

`check_cleaved_ends.py`, 69 genomes:

    GP_N start vs G start          +0 nt in 35/40 (88%)
    GP_C start vs GP_N end         +0 nt in 39/40 (98%)   <- the furin cut
    GP_C end   vs G end            +0 nt in 41/41 (100%)

## L is a spliced feature, and the junction is the literature's

Walker 2000 (**PMID 10756058**) states it for strain V / NC_001607 verbatim:
transcription unit III initiates at nt 1889, **intron II spans nt 2410–3703**,
and its removal *"results in fusion of 17 nt to the 5′ end of the conserved
sequence at nt 3704 and provides an in-frame AUG (nt 2393)"*, giving the
1,711 aa, 190 kDa protein that the same paper precipitates from infected cells.

**Measured independently here, and it reproduces the paper.** Every deposited L
was required to match its own genome EXACTLY; 107 orthobornavirus junctions
round-trip with:

    canonical GT..AG donor/acceptor      107 / 107
    first exon                           17-18 nt   (published: 17)
    intron                            1,282-1,294 nt (published: 1,294)
    donor context    TCATGCGAG|GT  in 77 of 107, 16 variants in all
    acceptor context TTTCCTTACCAG|CC in 73 of 107, 22 variants in all

The product begins **MSFHAS** in strain V — exon 1 gives M-S-F-H-A plus two
bases, and the first base of exon 2 completes the sixth codon — and MSFHAR in
lineages where that base differs.

**L is collinear in the other three genera**, 11 of 11 Carbovirus, Cultervirus
and Cartilovirus deposits. One feature cannot be both spliced and collinear, so
those genera are given references with a **zero-length intron** (donor last nt
= *n*, acceptor first nt = *n*+1), which makes `get_splice_variant_features.pl`
return the collinear ORF unchanged. That keeps the family on one feature with
no duplicate calls.

**120 references ship**, every one round-tripped before writing: cut at its own
declared coordinates, translated, and required to reproduce a clean
Met-initiated ORF of the expected length. Four candidates failed that check and
were dropped rather than shipped. `splice_id` / `splice_cov` are **90/90**,
measured: at the default 95/95 the references reach 94.1% of orthobornavirus
genomes and at 90/90, 95.1%.

### What this bought, measured on the same panel

| | PSSM CDS (rejected) | spliced feature (shipped) |
|---|---|---|
| calls identical to the genome's own deposit | **8 / 30** | **29 / 31** |
| worst N-terminal error | **+46 residues** | none |
| call rate on routed near-complete exemplars | 97.2% | 91.7% |

The rate falls by two genomes and the accuracy stops being approximate. Three
lineages have no usable reference: Parrot bornavirus 3 1367 (a 6,115 nt contig
that cannot hold a 5.1 kb gene), **Canary bornavirus 1** (its only L deposit
does not match its own contig cleanly, so no junction round-trips) and
**Carlia munda bornavirus**, whose L begins at a near-cognate **TTG** codon with
residues 2–1730 collinear — a non-AUG start this module does not model.

**A trap worth recording.** A tblastn collinearity screen of all 745 collection
members against their own genomes reported *every one* as a single collinear
HSP, L included, and the first conclusion drawn from it — "nothing in this
family is spliced" — was wrong. The genomic frame upstream of the acceptor
reads `LCFLTS` where the protein reads `MSFHAS`: five residues out of 1,711,
invisible at a 95% identity threshold. Only an exact match found it.

## N′ (p38) is modelled, and separated from N by a measured cutoff

Kojima 2019 (**PMID 30541858**) uses the names plainly: *"full-length isoform N
(p40) and N-terminally truncated isoform N′ (p38), which is translated from the
second initiation codon downstream of the NLS"*, and their second-AUG knockouts
are named **M14A**. In NC_001607 the N CDS is 54..1166 and the second in-frame
AUG is nt 93–95 = **Met14**, so p38 is 357 aa and lacks the 13 residues
MPPKRRLVDDADA that carry the PKRR nuclear localization signal of
**PMID 9188580**.

**Derived by homology to that measured position, not by a "second Met" rule.**
The first attempt cut each N member at its own second in-frame Met; that
disagreed with the profiles on six genomes whose N carries an extra Met at
residue 8, and `check_met_start.py` caught it as two alignments trimmed 6
residues into their members. The collection is now built by aligning all 129 N
members with MAFFT, locating the column of BoDV's Met14 and cutting at the Met
in that column: **112 of 129 qualify**, at residue 14 in 81, 16 in 23, 18 in 4,
17 in 2, 15 in 1 and 13 in 1.

**The cutoff is the separation, and it was measured.** N′ is an exact
C-terminal substring of N, so its profiles also match the body of an N that has
no second AUG. At the generated cutoff of 60 that produced three Met-less
calls — two carboviruses and Carlia munda. Best N′ score per exemplar:

    25 genuine Met-initiated calls      510 - 727 bits
    3 false calls                       380, 117, 108
    genera lacking the isoform           65 - 79

`bit_cutoff: 440` sits in that gap, 1.16× above the worst false call and 1.16×
below the weakest true one. After it: **25 calls, every one Met-initiated and an
exact C-terminal substring of its own N, cut 12–17 residues in; zero false
calls**; and the 11 genomes with no call are the Cultervirus, Carbovirus,
Cartilovirus and Carlia munda genomes that genuinely lack the isoform.
`upstream_ext` is 0 — scanning to a Met would walk back to p40's own start.

**That 1.16× margin is the narrowest in this module** and is stated in §07.

## There is no G′, and the G/L fusion is an mRNA

Both were asked about by the curator and checked against PubMed and PMC.

**No "G′" exists.** The only G-derived proteins in the literature are the
precursor and its two furin subunits, under two naming conventions: **GP-N**
(= GP1) and **GP-C** (= GP2 = gp43). Both are experimentally real — GP-N went
undetected until **PMID 12417322** found it in infected cells at 45–50 kDa
glycosylated, 27 kDa deglycosylated — and both are already shipped here as
`G_N` and `G_C`, whose gene symbols are now `GP_N` and `GP_C` to match.

**A G/L fusion is described, as an mRNA only.** Intron III shares intron II's
donor at nt 2410 and uses an acceptor **SA3 at nt 4559**; frame 1 of that
spliced RNA encodes *"the first 58 amino acids of G protein fused to L protein
lacking one-quarter of the N terminus"* (**PMID 11070091**), which is Cubitt's
**p165** (**PMID 11172106**). It is not modelled, for two reasons: the fusion
was **never detected as a protein** (only the small 75 aa S1 ORF was expressed,
FLAG-tagged, in transfected COS-7 cells, and intron III splices at ~10%), and
**its acceptor is not conserved** — PMID 11602780 reports both donors and both
major acceptors are highly conserved across BoDV isolates *"whereas a minor
alternative splice acceptor site is not"*, and that minor site is SA3.

## The glycoprotein signal peptide: SignalP is installed but cannot run

SignalP 6.0 is installed on this machine in its own conda environment and
answers `--help`, so the earlier note that no licence was available was wrong.
It still cannot predict: `site-packages/signalp/model_weights/` contains
nothing but its README, and every prediction dies with

    FileNotFoundError: Fast mode requires model to be installed at
    .../model_weights/distilled_model_signalp6.pt

The weights ship inside the licensed tarball and are installed separately by
`signalp6_convert_models`. Until that is done the signal peptide stays
unmodelled and `G_N` remains the N-terminal furin fragment *including* its
signal peptide. The queries are already staged at
`Alignments/Bornaviridae/G_SIGNALP_QUERIES.faa` (7 consensus N-termini), so
this is one command away from being finished. SKILL.md now tells the next
build to check the weights directory before planning around SignalP.

## What is deliberately not modelled

Listed with its reason in `notes/LITERATURE.md`: the **spliced N isoforms**
N2/N2'/N3/N3'/N3C of PMID 30541858 (N2 and N2' were not detectable in infected
cells by their own authors; N3/N3C only in a fractionated ER preparation of a
recombinant virus), the **intron-III products p8.4 and p165** of PMID 11172106
and PMID 11070091 (mRNA only, ~10% splicing, non-conserved acceptor),
**Cartilovirus Vp3** and the reptilian cultervirus **P–G hybrids**
(predictions in their own papers), the **1,608 aa and 1,558 aa internal-start
forms of L** (Walker 2000 looked for the 170 kDa one in infected cells and did
not find it), the **glycoprotein signal peptide** (SignalP installed but
weightless, above), and the avian carboviruses that lack G/X/M, which are not
in this dump.

N p38 was on this list in the first build and is no longer: it is modelled, as
`N_PRIME`.

## Vocabulary

Reused from `assets/annotation-vocabulary.tsv` where one existed:
`Nucleoprotein`, `Phosphoprotein`, `Matrix protein`,
`RNA-dependent RNA polymerase`, `Uncharacterized lineage-specific protein`.
New strings proposed: `X protein`, `Glycoprotein precursor`,
`Glycoprotein N-terminal subunit`, `Glycoprotein C-terminal subunit`.
