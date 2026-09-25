# Alsuviricetes module partitioning

Settled 17 September 2026, before any collection was built, from three sources:

1. **The dump** — 87,179 BV-BRC genome records for the class outside Togaviridae
   and Matonaviridae, 142,697 PATRIC features, 69,231 distinct protein sequences.
2. **Coordinates** — every CDS start/end, so gene order is measured, not assumed
   (`Alsu.layouts.txt`, built from the 764 RefSeq reference isolates).
3. **ICTV** — the current report chapters, for what the genome *is* rather than
   what was submitted.

The rule applied throughout is the kit's: **a module is a group of taxa that
share a genome organisation, not a clade.** Segment count and accessory gene set
split a module; phylogeny on its own does not.

## Already built, not revisited

| Module | Where |
|---|---|
| Togaviridae | `../Togaviridae` |
| Matonaviridae | `../Matonaviridae` |

## Segments are stored one record per genome

Every BV-BRC record in this class has `contigs = 1`, including the tripartite
Bromoviridae and the four-segment Blunervirus. A three-segment ilarvirus is
three genome records sharing one `genome_name`. So `merge_segments.py` is not
optional here the way it was for a monopartite family, and `"segments": N` in
each module block is load-bearing.

## The modules

### Hepelivirales

| Module | Taxa | Segments | Layout | Depth (uniq seq) |
|---|---|---|---|---|
| `Hepeviridae` | Paslahepevirus, Orthohepevirus, Rocahepevirus, Avihepevirus, Chirohepevirus, Piscihepevirus, Hepevirus | 1 | ORF1 nonstructural polyprotein, ORF2 capsid, ORF3 phosphoprotein (+ORF4 in Rocahepevirus and HEV-1) | 12,299 |
| `Benyviridae` | Benyvirus + unclassified | 4 (2-5) | RNA1 replicase; RNA2 CP + CP-readthrough + TGB1/2/3 + cys-rich; RNA3 p25; RNA4 p31 | 1,094 |
| `Alphatetraviridae` | Betatetravirus (1 seg), Omegatetravirus (2 seg) | 1 / 2 | RdRp + capsid | 78 |
| `Mycoalphaviridae` | Alphasclernavirus, Betasclernavirus | 1 | polyprotein (+CP in Betasclernavirus) | 14 |

### Martellivirales

| Module | Taxa | Segments | Layout | Depth |
|---|---|---|---|---|
| `Closterovirus` | Closterovirus | 1 |
| `Ampelovirus` | Ampelovirus | 1 |
| `Velarivirus` | Velarivirus, Olivavirus, Bluvavirus, Menthavirus | 1 |

**Closteroviridae was split three ways on 24 September 2026**, on the same
measured grounds as Quinvirinae. All six genera share the organisation
(ORF1a L-Pro/MTR/HEL, ORF1b RdRp by -1 frameshift, p6, HSP70h, HSP90h/p61,
CPm, CP), so the organisation rule alone would keep them together.

| grouping | genomes | 25 refs |
|---|---|---|
| lumped | 947 | **77.0%** |
| two-way: Closterovirus+small / Ampelovirus | 509 / 438 | 88.6% / 89.5% |
| **three-way: Closterovirus** | 333 | **95.8%** |
| **three-way: Ampelovirus** | 438 | **89.5%** |
| **three-way: Velarivirus + Olivavirus + Bluvavirus + Menthavirus** | 176 | **99.4%** |

Three beats two for every group and a fourth module adds nothing over the
three-way. Old row, for reference: | `Closteroviridae` | ... | 1 | ORF1a (L-Pro/MTR/HEL), ORF1b RdRp (-1 frameshift), p6, HSP70h, HSP90h/p61, CPm, CP, 3' accessories | 11,645 |
| `Crinivirus` | Crinivirus | 2 | RNA1 ORF1a/1b + p26 + p6; RNA2 HSP70h, p60, CP, CPm, p26 | 1,968 |
| `Bromoviridae` | Cucumovirus, Ilarvirus, Alfamovirus, Bromovirus, Anulavirus, Oleavirus | 3 | RNA1 1a (MTR+HEL); RNA2 2a RdRp (+2b); RNA3 3a MP + CP | 7,084 |
| `Tobamovirus` | Tobamovirus | 1 | 126K replicase / 183K **readthrough**, MP, CP | 3,068 |
| `Tobravirus` | Tobravirus | 2 | RNA1 134K/194K **readthrough**, MP, 16K cys-rich; RNA2 CP | 401 |
| `Furovirus` | Furovirus | 2 | RNA1 replicase + **readthrough** + MP; RNA2 CP + **CP-readthrough** + cys-rich | 240 |
| `Pomovirus` | Pomovirus | 3 | RNA1 replicase + **readthrough**; RNA2 CP + **CP-readthrough**; RNA3 TGB1/2/3 + cys-rich | 384 |
| `Hordeivirus` | Hordeivirus | 3 | alpha replicase; beta CP + TGB1/2/3; gamma RdRp-**readthrough** + cys-rich | 120 |
| `Peclu_Goravirus` | Pecluvirus, Goravirus | 2 | RNA1 replicase (+readthrough); RNA2 CP + TGB1/2/3 + cys-rich | 77 |
| `Endornaviridae` | Alphaendornavirus, Betaendornavirus | 1 | one polyprotein, no CP, no MP | 671 |
| `Cilevirus` | Cilevirus | 2 | RNA1 replicase + p29; RNA2 p15, p61 glycoprotein, MP, p24 | 438 |
| `Higrevirus` | Higrevirus | 3 | RNA1 polyprotein; RNA2 p50/p39/p9/p6; RNA3 p33/p29/p23 | 46 |
| `Blunervirus` | Blunervirus | 4 | RNA1 MTR-HEL; RNA2 HEL-RdRp; RNA3 accessories + virion membrane protein; RNA4 MP | 333 |
| `Idaeovirus` | Idaeovirus | 2 | RNA1 replicase; RNA2 MP + CP | 210 |
| `Pteridovirus` | Pteridovirus | 2 | RNA1 replicase (+p12); RNA2 MP + ORF2 + CP | 28 |

### Tymovirales

| Module | Taxa | Segments | Layout | Depth |
|---|---|---|---|---|
| `Alphaflexiviridae` | Potexvirus, Allexivirus, Mandarivirus, Lolavirus | 1 | replicase, TGB1/2/3, CP (+40K, +NABP in Allexivirus/Mandarivirus/Lolavirus) | 6,614 |
| `Botrexvirus_Platypuvirus_Sclerodarnavirus` | Botrexvirus, Platypuvirus, Sclerodarnavirus | 1 | replicase, CP, no TGB; fungal/non-virion | 42 |
| `Carlavirus` | Carlavirus | 1 | replicase, TGB1/2/3, CP, NABP | see note |
| `Quinvirinae` | Foveavirus, Robigovirus, Banmivirus, Ravavirus, Sustrivirus | 1 | replicase, TGB1/2/3, CP, NABP | see note |

**Quinvirinae was split on 24 September 2026, on measurement rather than on
organisation** — all six genera share the replicase/TGB1-3/CP/NABP layout, so
the organisation rule alone would keep them together. Lumped, 25 references
reach 81.0%. The cost is not confined to the divergent genus:

| genus | total | uncovered lumped | coverage lumped | coverage split |
|---|---|---|---|---|
| Foveavirus | 794 | 43 | 94.6% | **100%** |
| Carlavirus | 714 | 246 | 65.5% | **76.3%** |
| Robigovirus | 121 | 22 | 81.8% | **100%** |
| Banmivirus | 61 | 7 | 88.5% | (with Foveavirus) |
| Ravavirus, Sustrivirus | 4 | 4 | **0%** | (with Foveavirus) |

The mechanism is budget competition. Greedy largest-cluster-first gives **13 of
the 25 references to Carlavirus**, which has 102 clusters over 714 genomes and
still reaches only 65.5%, while Foveavirus — 21 clusters, 92% on five
references — gets 7 and falls from 100% to 94.6%. Splitting improves every
genus including Carlavirus, which gains 11 points purely by no longer sharing
a budget it cannot spare.

Same precedent as Tymoviridae (74.6% lumped; 100%, 91.2% and 100% per genus)
and Allexivirus.
| `Trivirinae` | Trichovirus, Vitivirus, Citrivirus, Prunevirus, Tepovirus, Chordovirus, Divavirus, Wamavirus | 1 | replicase, MP, CP, NABP | 6,688 |
| `Capillovirus` | Capillovirus | 1 | ORF1 polyprotein with **CP fused at its C-terminus**, ORF2 MP | 1,308 |
| `Tymovirus` | Tymovirus | 1 | ORF1 MP (overlapping), ORF2 replicase polyprotein, ORF3 CP | 296 |
| `Marafivirus` | Marafivirus | 1 | ORF1 polyprotein with **CP at its C-terminus**, ORF2 CP | 586 |
| `Maculavirus` | Maculavirus | 1 | replicase, CP, two small ORFs | 355 |
| `Deltaflexiviridae` | Deltaflexivirus | 1 | polyprotein + 3-4 small ORFs, no CP | 145 |
| `Gammaflexiviridae` | Mycoflexivirus, Gammaflexivirus, Xylavirus | 1 | replicase, (MP), CP | 59 |

## Why these splits and not others

**Collapsed, on the Rhabdoviridae precedent.** `Alphaflexiviridae` puts four
genera in one module because all four are monopartite replicase-TGB1/2/3-CP.
Allexivirus's 40K and the NABP are accessories, and `bit_cutoff` keeps them from
firing on a potexvirus; routing does not have to. Same argument for the six
genera in `Quinvirinae` and the six in `Closteroviridae`.

**Split on segment count**, which is not negotiable — the module block declares
`segments`, and merging is gated on it. Crinivirus leaves Closteroviridae for
this reason alone; Virgaviridae breaks into six modules for this reason mostly.

**Split on a genuinely different architecture.** `Capillovirus` is the clearest
case and it was measured, not assumed: BLASTp of Trichovirus/Tepovirus coat
proteins against Capillovirus ORF1 hits residues 1909-2100 of 2105
(28-33% identity, e-values to 5e-22). Capillovirus does not encode a separate CP
CDS, so it cannot share `Trivirinae`, whose CP feature is one.
`Marafivirus` is the same shape in Tymoviridae and splits from `Maculavirus` and
`Tymovirus` for the same reason.

## Declared but expected to be thin

`Alphatetraviridae` (78 unique sequences, and its two genera differ in segment
count), `Mycoalphaviridae` (14), `Gammaflexiviridae` (59), `Pteridovirus` (28),
`Botrexvirus_Platypuvirus_Sclerodarnavirus` (42) and `Higrevirus` (46) are all below or near the
kit's ~5-genome-per-genus heuristic once segments are merged. They are listed
here so the gap is on the record; whether each ships is a data question settled
at step 3, not now.

## Not modelled

371 genome records sit in the class with no family assignment — "martelli-like",
"tymo-like" and "hepe-like" singletons from metagenomic surveys, one species per
record, 528 unique proteins between them. No group of them reaches the member
floor. They go to `NOT_MODELLED.tsv`.

## Special features to settle before the JSON (step 6b)

These are the ones the coordinate data already flags. Each is checked with the
tblastn two-HSP test before anything is declared.

| Taxon | Symptom | Expected treatment |
|---|---|---|
| Tobamovirus, Tobravirus, Furovirus, Pomovirus, Hordeivirus, Pecluvirus, Benyvirus | replicase pairs 126K/183K, 134K/194K etc. differing by a fixed amount | leaky amber **readthrough** -> `internal_stop: 1` on the long form |
| Furovirus, Pomovirus, Benyvirus, Hordeivirus | CP and a much longer CP-readthrough at the same start | `internal_stop: 1` on the readthrough form |
| Closteroviridae, Crinivirus | ORF1a and ORF1b annotated separately, and as an "ORF 1a/1b fusion polyprotein" | **-1 ribosomal frameshift**; test, then `special: transcript_edit` or two collinear CDS |
| Capillovirus, Marafivirus | CP inside the ORF1 polyprotein | mature peptide, `cleave=` |
| Hepeviridae | RefSeq annotates MeT / Y / PCP / PPR / X / Hel / RdRp on ORF1 | domains, not established mature peptides -- see below |

## Settled: Hepeviridae ORF1 is not cleaved

Asked and answered from the primary literature, 17 September 2026, not from the
source annotation.

RefSeq annotates seven domains on the 1693-aa ORF1 of NC_001434
(methyltransferase, Y, papain-like protease, poly-proline hinge, X, helicase,
RdRp). That annotation descends from **Koonin et al. 1992, PNAS
(PMID 1518855)**, which is titled "*Computer-assisted* assignment of functional
domains" and describes them throughout as "*putative* functional domains". It is
a sequence-comparison prediction and was never a demonstration of processing.

The experimental record since then does not support cleavage:

- **Perttilä, Spuul & Ahola 2013, J Gen Virol (PMID 23255617)** expressed
  full-length pORF1 and six truncations in HeLa and Huh-7 by several vector
  systems, with antisera against the MT, HEL and POL domains. "The weight of
  evidence supports the proposition that pORF1 is not subjected to specific
  proteolytic processing, which is unusual among animal positive-strand RNA
  viruses but common for plant viruses."
- **LeDesma et al. 2023, eLife (PMID 36852909)** scanned the putative PCP domain
  by mutagenesis: "ORF1 operates as a multifunctional protein, which is not
  subject to proteolytic processing." The six essential cysteines coordinate
  divalent metal ions and hold the fold together; the domain rescues replication
  only in the context of full-length ORF1, never as an isolated subdomain.

So the module declares **ORF1 whole**, and the seven domain boundaries are
recorded in the build notes with coordinates rather than shipped as features.
Since these three papers were read, they go in `PMID`, not
`PMID_claude_generated`.

A useful cross-check fell out of Koonin 1992: by polymerase and helicase, HEV
groups with **rubella virus and beet necrotic yellow vein virus** -- Matonaviridae
and Benyviridae. The module neighbourhood this build assumes is the one the
founding paper describes.
