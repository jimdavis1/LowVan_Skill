# Tobaniviridae — build notes

Built 18 November 2026 with the LowVan kit at `jimdavis1/LowVan_Skill@d0689ef`.
Four modules, one per subfamily, on branch `tobaniviridae`.

Local build tree: `~/lowvan_scratch/tobani` (local disk; Box is synced at the
end, per the work-local rule). Runtime repo used for every measurement here:
`~/lowvan_scratch/tobani/VA_tob`, which holds **this family and nothing else** —
the four modules below — as SKILL.md §11 requires of a routing repo.

---

## 1. The dump

Fetched from the BV-BRC API in batches (`dump_api.py`, `fetch_contigs.py`):

| | |
|---|---|
| genome records in family Tobaniviridae | 884 |
| features | 1,763 (1,624 protein-coding) |
| unique protein sequences | 1,115 |
| per-genome contig FASTAs, 1-50 kb | 353 |
| `check_dump.py` | **PASS** on the family dump and on all four split dumps |

139 `id_md5` rows carry a blank md5 (5'UTR, 3'UTR, misc_feature, gap,
polyA_site, regulatory). They are not protein-coding and are correctly absent.

### Out of scope at intake

125 records whose organism name is `<taxon> sp.` were dropped before any
collection was built, per the curator's intake rule: 120 `Serpentovirinae sp.`,
4 `Torovirus sp.`, 1 `Bafinivirus sp.`

**Three exceptions were kept and this is a decision to review.** BV-BRC files
three ICTV species-assigned isolates under a `Serpentovirinae sp. isolate …`
organism name while giving them real species labels — `3071291.3` *Sertovirus
cona*, `3071292.3` *Septovirus foka*, `3071293.3` *Infratovirus latu*. They are
named viruses with assigned species, which the rule is written to keep; the
`sp.` string is in the isolate name rather than in place of a species. Each is
the only genome of its genus in the dump, so dropping them would delete three
of the seven serpentovirus genera from the module.

---

## 2. Module partitioning — four modules, settled on genome organisation

Every genome was assigned to a module **from sequence**, never from the BV-BRC
genus label: the 147 near-complete genomes were grouped by replicase-protein
identity (blastx, HSP-length-weighted, single linkage at 40%), and every other
genome took the module its proteins' best blastp hits agreed on.

| module | subfamily | genomes | genome organisation |
|---|---|---|---|
| `Torovirinae` | Torovirinae (*Torovirus*) | 629 | ORF1a–ORF1b–**S–M–HE–N** |
| `Serpentovirinae` | Serpentovirinae (7 genera) | 43 | ORF1a–ORF1b–**S–ORF3–M–N–ORF5b–ORF6** |
| `Piscanivirinae` | Piscanivirinae (*Bafinivirus*, *Oncotshavirus*) | 25 | ORF1a–ORF1b–**S–M–N**, no HE, no accessory |
| `Remotovirinae` | Remotovirinae (*Bostovirus*) | 18 | ORF1a–ORF1b–**S–M–N–G2–ORF6**, no HE |

The split is by gene layout, and the layouts differ: HE is a torovirus gene,
the bostovirus G2 and 88-residue ORFs exist nowhere else, and the serpentovirus
accessories (ORF3, ORF5b, ORF6) have no counterpart outside that subfamily.
*Bafinivirus* and *Oncotshavirus* share one module because they share the
layout and their replicases are 96.7–98.8% identical.

**The BV-BRC genus labels do not survive this test**, which is why they were not
used: three genomes labelled *Torovirus* (`2116382.3`, `2116383.3`,
`2116384.3` — snake "toroviruses" from Guangdong) group with Serpentovirinae at
76–100% replicase identity and carry the serpentovirus gene layout, and 120
`Serpentovirinae sp.` records group there too.

### Tabled, with the reason

Two singleton lineages belong to no module and are **not shipped**:

- `2916852.3` *Glycera tridactyla tobanivirus* (34,177 nt) — nearest relative
  anywhere in the family is bovine nidovirus at **28.1%** replicase identity.
- `2916857.3` *Salmo salar tobanivirus* (15,461 nt) — 42.4% to a torovirus, and
  no slippery heptamer in phase 2 within 150 nt of its ORF1a stop.

One genome each, no close relative, and a module built on one genome describes
that genome. They are a documented gap, not a shortfall.

---

## 3. Routing and the rep-contig budget

Measured before building, at the project budget of 25 references, greedily
largest-cluster-first, on genomes ≥15 kb:

| module | near-complete genomes | clusters at 70% | coverage at 25 refs | refs shipped |
|---|---|---|---|---|
| Torovirinae | 51 | 3 | **100%** | 3 |
| Serpentovirinae | 30 | 16 | **100%** | 16 |
| Piscanivirinae | 21 | 8 | **100%** | 8 |
| Remotovirinae | 18 | 1 | **100%** | 1 |

Nothing was tabled for routing. `build_rep_contigs.py` reproduces the budget
figure exactly, and each module ships only the references its own coverage
needs rather than the full 25.

---

## 4. Features, with the literature behind each one

Every citation below was resolved against PubMed and the title checked against
the claim. **All of them are `PMID_claude_generated`; `PMID` is empty on every
feature**, as it must be until a human curator reads them.

### Torovirinae — 7 features

| key | annotation | symbol | type | evidence |
|---|---|---|---|---|
| `ORF1A` | ORF1a polyprotein | ORF1a | CDS | PMID 16571831 (EToV pp1a 4,569 aa), PMID 33800523 (review: 13 nsps) |
| `ORF1AB` | ORF1ab polyprotein | ORF1ab | **transcript_edit** | PMID 2388833 (−1 frameshift, 12 nt overlap, pseudoknot), PMID 29950409 (U_UUA_AAC; 27.5–29.9% efficiency by ribosome profiling) |
| `MPRO` | Main protease involved in polyprotein cleavage | Mpro | **mat_peptide**, derived | PMID 16571831, PMID 32727876 |
| `S` | Spike protein | S | CDS | PMID 2219698 (BEV peplomer, 1,581 aa), PMID 23420207 (haemagglutination) |
| `M` | Membrane glycoprotein | M | CDS | PMID 2024492 (triple-spanning, 233 aa) |
| `HE` | Hemagglutinin-esterase glycoprotein | HE | CDS | PMID 9188596, PMID 15507445, PMID 19721004, PMID 34817201 |
| `N` | Nucleocapsid protein | N | CDS | PMID 2607340 (160 aa), PMID 33177195 (localisation signals) |

### Piscanivirinae — 6 features

`ORF1A`, `ORF1AB` (transcript_edit), `MPRO` (mat_peptide, derived), `S`, `M`,
`N`. Genome organisation and all four structural proteins: PMID 16987966 (white
bream virus), PMID 22422065 (fathead minnow nidovirus). Mpro: PMID 21068254.
The Chinook salmon bafinivirus N is the only *Oncotshavirus* protein with a
paper of its own (PMID 42128078, an immunogen study).

### Remotovirinae — 7 features

`ORF1A`, `ORF1AB` (transcript_edit), `S`, `M`, `N`, `G2`, `ORF6`. One genome
paper for the whole genus: PMID 25918239 (bovine nidovirus TCH5, KM589359).
`G2` and `ORF6` carry **Uncharacterized lineage-specific protein** — G2 is a
putative accessory glycoprotein with no demonstrated function and the 88-residue
ORF has no similarity to anything. The published symbol `G2` is kept; `ORF6` is
this module's own key, since the paper gives that ORF no name.

### Serpentovirinae — 14 features

`ORF1A`, `ORF1AB` (transcript_edit), `S`, `M`, `N`, and `UNC1`–`UNC9`, all nine
of them **Uncharacterized lineage-specific protein**. Genome organisation:
PMID 25205093 (ball python nidovirus, NC_024709). No function is published for
any serpentovirus accessory ORF — that paper says so — and `UNC1` includes the
ORF reported to resemble the paramyxovirus haemagglutinin-neuraminidase
(PMID 25063552), which is a sequence similarity and not a demonstrated
function, so it keeps the uncharacterized string like the rest.

---

## 5. Not modelled, and why

| what | where | why |
|---|---|---|
| pp1ab / ORF1b as a PSSM feature | all four | ORF1b has no start codon of its own; the product is the frameshift fusion, called from nucleotide references |
| S1 / S2 subunits | Torovirinae | cleavage into two subunits is experimental (PMID 2219698) but **the scissile bond has never been mapped**; "furin site" is a motif assignment in the review literature. Literature coordinates or nothing |
| HE and S signal peptides | all | no signal-peptidase site has been mapped in any tobanivirus; the bostovirus S site (aa 18/19) is a prediction (PMID 25918239) |
| nsp1–nsp13 mature peptides | Torovirinae | only the main proteinase's own two boundaries are Edman-confirmed. The other 11 EToV sites are predicted from the FxxQ↓(S,A) consensus, and the 10 non-self PToV sites were validated only in a surrogate biosensor, not mapped in the authentic polyprotein |
| any mature peptide | Serpentovirinae, Remotovirinae | **no cleavage site of any kind is published** for these subfamilies, and no protease from either has been expressed or assayed. Source-annotation ceiling: the module stops at the polyprotein |
| U1, U2 | Torovirinae | real, conserved, efficiently translated CUG-initiated ORFs (PMID 29950409), but **CUG-initiated**: the annotator scans to a Met start. One deposit each in the whole dump |
| 30K | Torovirinae | bioinformatic proposal (PMID 19737402), largely subsumed by U1/U2 |
| single-isolate accessory clusters | Serpentovirinae | 72 deposits in clusters holding **one distinct sequence** (one isolate, usually deposited twice as GenBank + RefSeq). A profile built from them describes one genome |
| replicase fragments | Serpentovirinae, Piscanivirinae | deposits under 2,500 aa named pp1a/pp1ab/ORF1b/"helicase"/"RdRp": pieces of a polyprotein, not proteins |
| `putative polyprotein 1c`, `putative envelope protein` | Piscanivirinae | one deposit each, one isolate each; no such protein is described for bafiniviruses |

---

## 6. The −1 ribosomal frameshift, and the junction this module ships

Both halves of this were measured here and then checked against the literature.

**The slippery heptamers.** Every near-complete in-scope genome was scanned for
an `X_XXY_YYZ` heptamer in phase 2 of its ORF1a frame within 150 nt of the
ORF1a stop. The measurement reproduces the literature exactly and splits the
family the same way the gene layouts do:

| module | heptamer | genomes | literature |
|---|---|---|---|
| Torovirinae | `UUUAAAC` | 49 | PMID 2388833, PMID 29950409 |
| Piscanivirinae | `UUUAAAC` | 19 | NC_008516 `misc_feature 14549..14555` (the paper's own text was not retrievable) |
| Remotovirinae | `UUUAAAC` | 18 | PMID 25918239 |
| Serpentovirinae | `AAAAAAC` | 28 | PMID 25205093, PMID 25063552 |

**The junction.** The two representations in circulation differ by two residues.
RefSeq duplicates the **last** base of the heptamer — NC_007447
`join(859..14190,14190..21059)` with the heptamer at 14184..14190, NC_008516
`join(906..14555,14555..21523)`, NC_024709 `join(1018..18411,18411..25343)`,
KM589359 `join(487..9825,9825..15818)` — giving `…Leu-Asn…` for U_UUA_AAC.
Many bovine and porcine torovirus submitters, and MH603532's own `join`, place
it a base earlier, giving `…Phe-Lys…`.

**This module ships the RefSeq construction**, and the reason is mechanistic
rather than a preference between conventions. −1 frameshifting at a shifty
heptamer is tandem slippage of the P-site and A-site tRNAs, so the residues
incorporated are the ones decoded from the **0-frame** P and A codons — UUA
(Leu) and AAC (Asn). The direct evidence is in SARS-CoV, not in a tobanivirus:
Baranov *et al.* 2005 (PMID 15680415) analysed a tagged SARS-CoV frameshift
product by mass spectrometry and concluded that frameshifting "occurs only due
to tandem slippage of the P-site tRNA-Leu and A-site tRNA-Asn … at the standard
U_UUA_AAC slippery site", and mutating either codon reduced frameshifting to
background. Brierley & Dos Ramos 2006 (PMID 16310880) state the same
assignment for the coronavirus ORF1a/1b junction. **No torovirus, bafinivirus,
bostovirus or serpentovirus pp1ab product has ever been sequenced or analysed
by mass spectrometry** — the junction here is inferred from the shared
mechanism at the *same* U_UUA_AAC heptamer, and that is the evidence level.

Corroboration that the construction, not just the mechanism, is right: applying
it to SARS-CoV-2 (NC_045512, heptamer 13462..13468) reproduces the nsp12
N-terminus `SADAQSFLNRVCGVSAARL` used in every nsp12 structure; the alternative
construction yields `SADAQSFFKRVCG`, which matches nothing.

**Serpentovirinae is a separate argument**, because its heptamer is
`A_AAA_AAC`, which Baranov did not test. There is no product-level evidence for
it in any serpentovirus. What carries it is the general tandem-slippage rule
established with product sequencing in other systems — Jacks *et al.* 1988
(PMID 2447506) localised the RSV/HIV-type site by amino-acid sequencing to the
0-frame P-site codon, and Brierley *et al.* 1989/1992 (PMID 2720781,
PMID 1404364) fixed the IBV site at UUUAAAC and showed the A-site codon
(positions 5–7) governs efficiency. Applied to A_AAA_AAC that gives P-site AAA
(Lys) and A-site AAC (Asn), i.e. `…Lys-Asn…`, which is exactly what NC_024709's
own `join` encodes. Evidence level: **mechanism-level inference plus the RefSeq
record, with no product sequenced for any serpentovirus.**

The records that place the junction a base earlier are using a different join
convention, not a different biology; where a genome's deposited pp1ab differs
from this module's, it differs by those two residues and nothing else.

**Reference sets.** `<Module>/Transcript-Editing/<Module>/ORF1AB.fasta`, built
from each in-scope near-complete genome: Torovirinae 44, Serpentovirinae 20,
Remotovirinae 15, Piscanivirinae 11. Rejected rather than guessed: 23 genomes
carrying an ambiguous base anywhere in the construct, 4 with no slip site in
the window, 3 with no ORF ≥2,500 codons, 1 whose ORF1b part was under 1,500
codons. The ORF1a start the builder picked agrees with the submitter's own
coordinate on **97 of 101** genomes.

---

## 7. The mature peptides: Mpro, derived from published cut sites

Nothing in BV-BRC is annotated as a tobanivirus mature peptide, so this feature
exists only because the cleavage sites are published. Every coordinate below
was re-read from the paper and then **checked against the actual deposited
sequence**; all matched with no offset.

**Torovirinae.** Equine torovirus (PMID 16571831, GenBank DQ310701, numbering
from Met1 of pp1a): N terminus **Gln3251↓Ser3252** in `SNFSFQ↓SVFSKA`,
C terminus **Gln3543↓Ser3544** in `QPFKKQ↓SVSNNV`, both by Edman degradation of
the cleavage products; mature protein **292 aa**. Porcine torovirus
(PMID 32727876, GenBank MH603532): **Gln3095↓Ser3096** in `SNFSFQ↓SVFSKV` and
**Gln3382↓Ser3383** in `QPFKKQ↓SVSTNV`, Edman plus Q→A mutagenesis; mature
3CLpro **pp1a 3096–3382 = 287 aa**. (An earlier draft of these notes wrote the
C-terminal site as "Gln3382|Ser3382"; it is 3382↓3383.) A third site,
**Gln3397↓Ala3398** in `SEFATQ↓AWQT`, cuts EToV Mpro in half in *E. coli* and
was not detected in infected cells, so it is not a feature here.

**Piscanivirinae.** White bream virus (PMID 21068254, NC_008516): all four
autoprocessing sites mapped by Edman — N1 **Gln3442↓Ala3443** (`ITRQ↓ALQR`,
predominant), N2 **Gln3449↓Ala3450**, C1 **Gln3709↓Thr3710** (`VTQQ↓TSVT`),
C2 **Gln3716↓Thr3717**; catalytic triad Ser3589-His3492-Asp3518. The paper's own
N1–C1 marker protein is pp1a **3443–3709**.

**Junction chemistry, measured before cutting** (the rule: measure, then cut):

| module | derived | P1 | P1′ | P1 (C-term) | P1′ (C-term) | lengths |
|---|---|---|---|---|---|---|
| Torovirinae | 32 of 32 pp1a | Gln 100% | Ser 100% | Gln 100% | Ser 100% | 287–292 aa |
| Piscanivirinae | 5 of 6 pp1a | Gln 100% | Ser 80% / Ala 20% | Gln 100% | Ser 80% / Thr 20% | 263–275 aa |

Torovirinae cuts at the published motif directly. Piscanivirinae cannot: the
white bream virus N1 context is not conserved in *Oncotshavirus*, so the
derivation anchors on the protein whose ends were sequenced (pp1a 3443–3709,
located by blastp) and snaps each boundary to the nearest Gln↓(Ala,Thr,Ser),
that protease's published consensus. The sixth sequence is **refused, not
guessed** — no such junction within 12 residues of the projected boundary.

**And the called products were checked on the panel**, which is the test that
matters: Torovirinae Mpro, 24 of 24 calls begin with Ser and end with Gln,
287–307 aa. Piscanivirinae, 13 calls, 77% begin Ser and 15% Ala (both are
published P1′ residues), 92% end in Gln; one call drifts at both ends.

`check_cleaved_ends.py`'s tiling test does not apply to a single **internal**
mature peptide — both of its offsets from the precursor are large by
construction — so the chemical test above is the one reported. The offset
itself is a spike: Mpro's C terminus sits at +3078 nt from ORF1a's in 21 of 24
genomes (88%), and within ±3 of that in 92%.

---

## 8. Splice variants

**None, in any genus.** Tobaniviruses are cytoplasmic +ssRNA viruses that make
a 3′-coterminal nested set of subgenomic mRNAs (PMID 2293666); toroviruses are
unusual in making only the S sgRNA discontinuously (PMID 12456663,
PMID 15956573). There is no spliced mRNA to model, and the only `join()` CDS in
any reference record is pp1ab via ribosomal slippage. No `special: splice`
feature is declared.

---

## 9. What the measurements say

(filled in from the coverage run — see §10 of the coverage audit page)

---

## 10. Open questions for the curator

1. **The three kept `sp.` isolates** (§1). They are the only genomes of
   *Sertovirus*, *Septovirus* and *Infratovirus latu* in the dump.
2. **`PMID` is empty on every feature.** Every citation is in
   `PMID_claude_generated` and needs a human to read it before it moves.
3. **Single-sequence profiles.** Most Serpentovirinae and Piscanivirinae
   features rest on one-sequence profiles (see §9); they are honest models of
   one sequence each, and the alternative was no profile at all.
