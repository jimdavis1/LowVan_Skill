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
(Lys) and A-site AAC (Asn), i.e. `…Lys-Asn…`. Evidence level: **mechanism-level
inference from product sequencing in other viruses, with no serpentovirus
frameshift measured and no serpentovirus pp1ab product sequenced.** NC_024709's
own `join` happens to encode the same two residues; that is agreement, not
evidence, and it is not what this construction rests on.

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

### Profiles

| module | features | profiles | single-sequence | collections |
|---|---|---|---|---|
| Torovirinae | 7 (1 `transcript_edit`) | 57 | 42 | ORF1A 32, S 62, M 37, HE 87, N 52, MPRO 32 |
| Piscanivirinae | 6 (1) | 24 | 19 | ORF1A 6, S 13, M 9, N 7, MPRO 5 |
| Remotovirinae | 7 (1) | 10 | 5 | ORF1A 16, S 12, M 8, N 5, G2 10, ORF6 9 |
| Serpentovirinae | 14 (1) | 85 | 83 | S 17, M 17, N 14, ORF1A 12, UNC1-9 2-7 each |

`-m` was lowered to 2, recorded in every `BUILD_PARAMS`, because at the
documented floor of 5 most of these features build **zero** profiles: 6 to 17
unique sequences per protein spread over several genera never reach five per
cluster. `-mi` was **not** lowered — the identity floor is what makes a profile
a model of one protein. The cost of holding it is that most sequences end up in
a cluster of one and get a single-sequence profile, built by duplicating the
sequence into a two-row alignment, which is what psiblast will accept. A
single-sequence profile is a lookup table with a bit score: it calls a close
relative and does not generalise.

### Held-out panels

Drawn from the build's own contigs with the reference genomes excluded by
accession. Rates are computed from the feature tables, by gene symbol — NOT
from `run_gto_eval.py`'s per-feature figure, which counts a feature as called
when its annotation STRING appears in the GTO and so reports one identical rate
for the nine Serpentovirinae features that share "Uncharacterized
lineage-specific protein".

**Torovirinae** (24 panel genomes)

```N          100.0%  (24/24)
   S          100.0%  (24/24)
   Mpro       100.0%  (24/24)
   M          100.0%  (24/24)
   ORF1a      100.0%  (24/24)
   HE          95.8%  (23/24)
   ORF1ab      87.5%  (21/24)
                8.3%  (2/24)
```

**Piscanivirinae** (13 panel genomes)

```N          100.0%  (13/13)
   S          100.0%  (13/13)
   M          100.0%  (13/13)
   Mpro       100.0%  (13/13)
   ORF1a       84.6%  (11/13)
   ORF1ab      61.5%  (8/13)
               23.1%  (3/13)
```

**Remotovirinae** (8 panel genomes)

```ORF1a      100.0%  (8/8)
   S          100.0%  (8/8)
   ORF6       100.0%  (8/8)
   M          100.0%  (8/8)
   N          100.0%  (8/8)
   ORF1ab      87.5%  (7/8)
   G2          25.0%  (2/8)
               12.5%  (1/8)
```

**Serpentovirinae** (14 panel genomes)

```S          100.0%  (14/14)
   M          100.0%  (14/14)
   ORF1a       92.9%  (13/14)
   N           92.9%  (13/14)
   ORF1ab      85.7%  (12/14)
   UNC1        57.1%  (8/14)
   UNC2        50.0%  (7/14)
   UNC3        42.9%  (6/14)
   UNC4        14.3%  (2/14)
   UNC6         7.1%  (1/14)
   UNC5         7.1%  (1/14)
   UNC7         7.1%  (1/14)
   UNC8         7.1%  (1/14)
   UNC9         7.1%  (1/14)
```


### Bit cutoffs, reset from measured scores

The generator derives a cutoff from the profile, and for a short protein that
lands at 30 bits. On Berne virus a second *Membrane glycoprotein* was called at
position 11,322 — inside the replicase ORF — scoring **30.34 against a cutoff
of 30**, while the genuine M on that genome scored 476.8. Every feature's
cutoff was reset to half the lowest score its own panel calls reached, rounded
down, never below the generator's value, and **capped at 500**: half of
Piscanivirinae ORF1A's lowest panel score is 3,000, and at 3,000 the divergent
fish tobaniviruses whose pp1a matches partially stop being called at all
(measured: 70 per cent falling to 33). Features with fewer than five observed calls — UNC4 to
UNC9 — keep the generator's value, because one observation says nothing about
how a divergent member of the same orthogroup would score.

The measured score distributions behind those cutoffs:

```
module           feature      n       min    median       max
Piscanivirinae   M           13     167.4     501.8     512.6
Piscanivirinae   MPRO        13     110.2     540.1     548.6
Piscanivirinae   N           13      66.7     324.4     367.1
Piscanivirinae   ORF1A       11    6960.5    9254.0    9747.1
Piscanivirinae   S           13    1687.2    2435.2    2525.3
Remotovirinae    G2           7     907.8     974.8    1004.8
Remotovirinae    M            8     497.5     502.9     502.9
Remotovirinae    N            8     355.5     357.5     357.5
Remotovirinae    ORF1A        8    6077.9    6670.0    6693.5
Remotovirinae    ORF6         8     177.9     180.6     185.6
Remotovirinae    S            8    1163.5    1168.9    1182.0
Serpentovirinae  M           14     441.7     474.7     567.7
Serpentovirinae  N           13     132.1     301.2     489.6
Serpentovirinae  ORF1A       13    1085.1   11039.4   13270.5
Serpentovirinae  S           14    1852.7    1967.5    2388.2
Serpentovirinae  UNC1         8     986.9    1072.7    1073.0
Serpentovirinae  UNC2         8      63.5     694.8     872.4
Serpentovirinae  UNC3         6     471.4     500.0     518.4
Serpentovirinae  UNC4         2     970.6    1015.2    1015.2
Serpentovirinae  UNC5         1     551.6     551.6     551.6
Serpentovirinae  UNC6         1     783.5     783.5     783.5
Serpentovirinae  UNC7         1     684.8     684.8     684.8
Serpentovirinae  UNC8         1     474.1     474.1     474.1
Serpentovirinae  UNC9         1    1454.4    1454.4    1454.4
Torovirinae      HE          23     841.6     887.1     916.7
Torovirinae      M           24     468.3     481.0     500.7
Torovirinae      MPRO        24     463.7     613.9     632.0
Torovirinae      N           24     288.2     312.8     426.5
Torovirinae      ORF1A       24    6703.5    9378.0    9761.6
Torovirinae      S           24    3145.9    3197.5    3322.6
```

After the change the false call is gone, Berne annotates ORF1a, Mpro, S, M and
N (correctly no HE: its HE is a pseudogene), and every panel rate above is
unchanged.

### Whole-taxon run

Every Tobaniviridae record BV-BRC holds between 1 and 50 kb (353 genomes),
clustered at 95% nucleotide identity, one exemplar per cluster, routed against
a repository holding this family and nothing else.

**Routing.** 91 of the 117 exemplars were annotated through the GTO pipeline
and **every one of them routed**. The other 26 never reached a module:
`rast-create-genome` returned HTTP 403 from the BV-BRC ID server, which the
coverage report then counts as UNROUTED — an infrastructure failure wearing a
routing failure's name. Re-run through the flat annotator, which needs no ID
server:

```
26 of 26 ID-server-failed exemplars annotated with the flat annotator
routed: 24/26
modules they route to: {'Torovirinae': 5, 'Piscanivirinae': 2, 'Serpentovirinae': 17}
```

So 115 of 117 distinct genome types route, 98.3%. The two that do not are a
1,625 nt and a 1,316 nt fragment.

**Quality, with an explicit completeness cutoff.** "Near-complete" is measured
against the genome each module expects, not a family-wide constant: an
exemplar counts as near-complete when it reaches 90% of the median length of
its own module's reference genomes (Torovirinae 28,511 nt, Piscanivirinae
27,004, Remotovirinae 20,386, Serpentovirinae 30,104). 49 of the 117 clear
that; the other 68 have a median length of 4,731 nt and **90.5% of them are
correctly flagged incomplete**.

The 49 were re-run against the module as it now ships:

```
near-complete exemplars re-run against the shipped module: 36
  routed 36 (100.0%), good 30 (83.3%), poor 6
  Serpentovirinae  16 genome(s): ORF1A 100%, M 100%, N 100%, S 100%, UNC2 69%, UNC1 69%, UNC3 56%, ORF1AB 44%, UNC4 12%, UNC8 6%, UNC9 6%
  Torovirinae      12 genome(s): ORF1A 100%, M 100%, S 100%, MPRO 100%, N 92%, HE 92%, ORF1AB 75%
  Piscanivirinae   4 genome(s): ORF1A 100%, M 100%, N 100%, S 100%, MPRO 100%, ORF1AB 100%
  Remotovirinae    4 genome(s): ORF1A 100%, M 100%, N 100%, S 100%, ORF6 100%, G2 50%, ORF1AB 50%

  poor, with the reason:
    2116384.3    Guangdong red-banded snake torovirus LPS ORF1a polyprotein: Feature is too short
    2116385.3    Hainan hebius popei torovirus LPSC33749  ORF1a polyprotein: Feature is too short
    2529475.6    Torovirus sp. 15211_TORO_len_26822_cov_3 Genome is missing essential feature: Nucleocapsid protein | ORF1a polyprotein: Feature is too short
    2661817.62   Serpentovirinae sp. H0-2                 ORF1a polyprotein: Feature is too short
    2816848.3    Bangali torovirus Bangali/H.rufipes/2018 ORF1a polyprotein: Feature is too long
    3403327.3    Sinibotia superciliaris tobanivirus 1    Nucleocapsid protein: Feature is too short
```

**The six poor ones are all profile coverage, and all visible.** Four carry
`ORF1a polyprotein: Feature is too short` — a divergent genome whose
polyprotein matches partially, which is what a thin profile set does; one is
missing N; and Bangali torovirus carries `Feature is too long` on ORF1a. None
is a silent failure, and the fix for all of them is more sequences rather than
different thresholds.

### Where ORF1AB misses, classified

| cause | n | what it means |
|---|---|---|
| ambiguous bases in the genome | 7 | the annotator demotes to `partial_cds` by design, because the gap-fill draws its bases from the subject; the same ambiguity is why those genomes are not references |
| no pp1a on the contig | 2 | 5'-truncated submissions |
| no reference inside the 95% gate | 2 | divergent, or a partial best HSP |

`check_te_leak.py` states the other half plainly: the reference sets are built
from every genome that yields a clean construct, so most panel genomes are
also references and ORF1AB's rate on them is self-recall. On the genomes that
are NOT references it is 2 of 13 — and 11 of those 13 are exactly the genomes
rejected from the reference set for ambiguity or partiality. The two
populations coincide, which is the honest reading.

---

## 10. Open questions for the curator

1. **The three kept `sp.` isolates** (§1). They are the only genomes of
   *Sertovirus*, *Septovirus* and *Infratovirus latu* in the dump.
2. **`PMID` is empty on every feature.** Every citation is in
   `PMID_claude_generated` and needs a human to read it before it moves.
3. **Single-sequence profiles.** Most Serpentovirinae and Piscanivirinae
   features rest on one-sequence profiles (see §9); they are honest models of
   one sequence each, and the alternative was no profile at all.
4. **The bit-cutoff rule is mine, not the kit's.** Half the lowest panel score,
   capped at 500, with features under five observed calls left alone. It was
   derived from one measured false positive and one measured over-correction
   (§9). If the project has a standing rule for this, use that instead.
5. **`ORF1AB` carries no `copy_num` in any module**, although it measures
   85.7–87.5% on the panels. The misses are genomes with ambiguous bases, which
   the annotator demotes by design; marking the feature essential would flag
   those genomes twice for one defect. Reversible in one line if you disagree.
6. **The Serpentovirinae `UNC` features are an orthogroup split, not a naming
   claim.** Nine features, all `Uncharacterized lineage-specific protein`,
   separated by homology with position used only to tell one group from
   another. If the literature later names any of them, the keys stay and the
   annotation changes.

## 11. Known limits of this build

- `check_rep_contigs.py` was not run: it fetches by accession from NCBI, and
  routing is already measured twice here — by the budget script before the
  build and by the whole-taxon run after it (§9).
- The whole-taxon pass is incomplete through no fault of the modules: the
  BV-BRC ID server returned HTTP 403 for 26 of 117 exemplars on the first pass
  and 13 of 49 on the near-complete re-run. Routing for the first 26 was
  recovered with the flat annotator; the quality figures rest on the 36
  near-complete exemplars that did annotate.
- `archive_workdir.py` was not run. The build tree is at
  `~/lowvan_scratch/tobani` on this machine, and the alignments behind every
  shipped profile are committed under `modules/<Module>/PSSM-Alignments/`.
