# Kolmioviridae — build notes

Built 25–26 November 2026. **Two modules**, `Deltavirus` and `Kolmioviridae`.

Cross-check done first, at Jim's request: **no existing LowVan module covered
delta hepatitis virus**, despite its former classification outside the family.

| module | features | genomes | refs | PSSMs |
|---|---|---|---|---|
| `Deltavirus` | `S-HDAg` 685, `L-HDAg` 852 | 1,176 | 6 | 20 |
| `Kolmioviridae` | `DAg` 29 | 41 | 17 | 16 |

## Why two modules

`DAg` and `S-HDAg` are the same protein, in different lineages, on the same
locus. Inside one module the annotator called both on **601 of 660 exemplars**
— on `12475.1049` all three features landed on one strand within 36 nt.

**No bit_cutoff can separate them.** Chosen-profile scores are fully
interleaved: `DAg` runs 35–442 bits on HDV against 157–432 on its own genomes.
A length or coverage threshold trades duplicates against genuine calls roughly
one-for-one, so that is the wrong instrument too.

**The per-genome winner was right 660/660**, and rep-contig routing reproduces
that split **659/659** under the annotator's own blastn
(`-reward 2 -penalty -3 -word_size 11`, best bit wins). So the family is
partitioned by reference set, which needs no change to the shared annotator.

An earlier note of mine said the divergent genera "all route to the Deltavirus
references". That was not a test of separability: with one module there was only
one destination. Re-measured after building the second reference set, they
separate cleanly.

**Partition rule: the genome name, not the taxon id.**
`^hepatitis\s+(delta|d)\s+virus` → `Deltavirus`, else `Kolmioviridae`. Taxon
prefix `12475` is wrong — HDV isolates also sit under `510xxx` and
`2847173/4`. Under the name rule the collections were already clean on both
sides, so the split was a pure relabelling with no re-binning.

**BV-BRC's genus column is stale here.** 13 divergent agents (termite ×5, bat
×2, rodent, fish, avian, *Lonchura*, *Epsilon*) are still filed under genus
`Deltavirus`; ICTV moved them. Scores and routing both agree with ICTV.

## The edited codon, and `transl_except`

The small antigen stops at an amber codon; ADAR1 edits it to UGG and the large
antigen runs 19 residues further to the CXXQ farnesylation box at Cys211. Of
the 852-member large-antigen collection, **508 carry TAG, 392 the IUPAC TRG,
143 TGG**. TRG is physical evidence of a mixed population; it translates as `X`
and is **not a stop**, so `internal_stop` never sees it. Hence a second,
independent mechanism, new in this build:

```json
"transl_except": {
  "replace_residue": ["X", "*"],
  "insert_residue":  "W",
  "offset":          [19, 20],
  "from":            "C-term"
}
```

Offsets are measured on the **emitted protein**, never the PSSM query, so the
JSON never names a profile. The window is 19–20 because genotype III's
extension is genuinely one residue longer (PMID 12097551). All 852 members
resolve to exactly one position; widening by 1 would make 2 ambiguous, by 2
would make 4. More than one candidate in the window → declined and logged, not
guessed.

**`S-HDAg` must keep `downstream_ext: 0`.** It ends AT the editing site, and on
an edited or ambiguous genome no stop codon marks that end. The S/L pair is
declared in `protected_pairs.json` — a strict prefix by biology.

## Cutoffs, measured from the annotator's trace

At the generated cutoff of 30, 29 genomes carried a duplicate call.

| feature | genuine, 1st pct | secondary, median | cutoff | cost |
|---|---|---|---|---|
| `S-HDAg` | 191.5 (min 100.5) | 33.2 | 100 | 0 genomes; 92/99 secondaries removed |
| `L-HDAg` | 215.3 | 35.1 | 100 | 1 fragment; 175/192 removed |
| `DAg` | 156.8 (min 70.9) | none observed | 60 | 0 genomes |

Duplicates fell 29 → 8. The 8 that remain are deposit properties: four
(AF01–AF04) carry the antigen **at full length on both strands**, 481 vs 468
bits. `DAg`'s cutoff had to be measured because **13 of its 16 profiles are
built from a single deposit**; HDV is kept away from them by routing, not by
that number.

## What the modules refuse, and why that is right

**8 exemplars are declined.** They route at 1,159–1,756 bits, but seven carry
an antigen ORF with **2–4 internal stop codons** — the translation is plainly a
delta antigen (`MSRSESRKNRGGREDILEQWVSGRKRL…`) but coverage before the first
stop falls to 0.07–0.49, under the 0.65 floor. The eighth is 27% N.

**108 of 631 Deltavirus genomes score Poor**, none of it a module defect:
85 contigs with >1% ambiguous bases (expected — the editing site itself is
deposited as `TRG` in 392 genomes), 30 features short, 12 HSP-count or
missing-essential, 7 long.

## Fixes made to the shared kit (committed as `f1495200`)

- `gen_module_json.py` dropped `transl_except` silently on regeneration.
- `evaluate_coverage.py` reached `fold_circular()` only from the special-feature
  branch, so a circular module with no special features counted every feature
  twice.
- `annotator_regression.py` parsed the `no_features_called` sentinel as a call
  and died on `int('')` after both passes had finished.

Regression: **45 genomes, 755 features, 0 differ** against `VA_hepadna`.

## Fixes made during this build

- `install_module.py` preserves source filenames, so `Deltavirus` was serving
  profiles named `Kolmioviridae.*`. Renamed at source.
- 13 single-sequence `DAg` profiles had no two-row MSA in `corrected_alis`, so
  `rebuild_pssms.py` reported FAILED rather than stale. Fixed; 36/36 reproduce.
- `protected_pairs.json` was in a shape and location `audit_bins.py` never
  reads (it wants a bare list of pairs beside the module's collections).
- The white-tailed deer deltavirus (*Perideltavirus* ov-mh1) routed to the HDV
  module by 12 bits. Added as `Kolmioviridae.17`; now 2,684 vs 489.
- The synonym table was family-wide, so both collapse pages showed HDV strings.
  Rebuilt per module from each module's own genomes.

## Final coverage (660 exemplars, 95% nucleotide identity)

| | `Deltavirus` | `Kolmioviridae` |
|---|---|---|
| routed | 631 | 21 |
| routed, called nothing | 0 | 0 |
| called a foreign module's feature | 0 | 0 |
| `S-HDAg` / `L-HDAg` / `DAg` | 630 / 628 | 21 |
| genomes Good | 523 (82.9%) | 21 (100%) |

Every profile recovers **100%** of its own collection: 852/852, 685/685, 29/29.

## Open

- The 8 deposits with disrupted ORFs are left uncalled by design. If BV-BRC
  ever corrects them they will annotate without any change here.
- `PMID_claude_generated` throughout. **No citation here has been curator-checked.**

---

## RESOLVED — ORF2 built, 27 Nov 2026

The defect below is closed. `ORF2` is now a declared feature of the
`Kolmioviridae` module: 34 collection members, 17 profiles, called on 36 of the
41 genomes. The `Deltavirus` module is untouched.

### What the paused investigation had wrong

**The ORF survey in the defect note below assigned DAg by position, not by
homology, and got it backwards in several genomes.** Re-running with DAg located
by blastp against the module's own `DAg` collection:

- the "196 aa ORF at nt 1, shared by 11 Rodent deltaviruses plus Tome's
  spiny-rat virus", cited below as the strongest sequence evidence for ORF2,
  **is DAg itself**;
- the rule used to dismiss candidates -- *antisense-nested inside DAg = shadow
  ORF, not real* -- **is the opposite of the truth**. Kishimoto et al. state
  verbatim that ORF2 is "present on the genomic-sense strand in the region
  complementary to the DAg ORF". Half of the 36 ORF2s found here lie wholly
  inside the antigen gene antisense, and that is the architecture, not an
  artefact. The HDV plus-strand ORF dismissed below on the same reasoning was
  dismissed for the wrong reason -- though the conclusion stands, since HDV's
  two antigens are the only proteins it is known to translate and the preprint
  does not claim ambisense expression in HDV.

### The literature, read

| source | what it establishes |
|---|---|
| PMID 42842300, J Gen Virol Oct 2026, full text read | SwSCV-1 has "two partially
  overlapping >500 base ORFs, the DAg ORF in the antigenomic orientation and the
  ORF2 in the genomic orientation"; 177 aa. The knockout is real: replication,
  persistence and particle formation unaffected, no ORF2 mRNA. **But the mutation
  reverted to wild-type Met during >1 yr persistence (92.3% G at that position)**,
  and the authors say they "cannot rule out that ORF2 would be beneficial". |
| bioRxiv 2025.01.08.631871 v3, Kishimoto, Imai, Sono, Igarashi, Sassa-O'Brien,
  Rubbenstroth, Horie | ORF2 of *Serinus canaria* kolmiovirus is expressed by
  ambisense transcription and promotes helper-dependent transmission with avian
  bornavirus helpers. **ORF2a 135 aa unspliced, ORF2b 149-150 aa spliced**, three
  genotypes (Supplementary Table S3, parsed locally from the media-2.xlsx
  workbook). Five further KoV lineages have genomic-sense ORFs with poly-A
  signals, incl. *Indirana chiravasi* deltavirus and DrDV-B. No PMID; cited by
  DOI only. |

### How the collection was built

ORF2 is **absent from BV-BRC entirely** -- 0 of 3,368 unique sequences. Every
member was translated from the module's own genomes, into
`Extracted/Kolmioviridae.ORF2.faa`.

Rule: locate DAg by blastp; call ORFs on the OPPOSITE strand, circular-aware,
>=60 aa; take the longest whose span is >=40% complementary to the DAg ORF.
36 of 41 genomes carry one; 65-212 aa, median **149** -- the scKoV ORF2b length.

**Three corrections after the first annotation pass, homology over position:**

| genome | rule picked | profiles called | action |
|---|---|---|---|
| 2767006.6 DrDV-B CAJ1 | 120 aa @1298-1660 | 73 aa @1018-1240 | replaced; the 73 aa
  ORF is the orthologue of the one in both sibling DrDV-B genomes |
| 2751480.5 Rodent deltavirus 1056 | 102 aa (unrelated ORF) | 204 aa, no Met start |
  dropped from training; the ORF2 ATG is inside a run of Ns |
| 2751480.7 Rodent deltavirus 0354 | 91 aa (unrelated ORF) | 197 aa, no Met start |
  dropped from training; same cause |

Both rodent genomes are already flagged "Contig has too many ambiguous bases".
The locus is still called in both. Collection: 36 -> 34.

### ORF2 is eleven proteins sharing a position, not one family

Masking the DAg-antisense portion of each ORF2 and blasting the remainder against
itself leaves matches only WITHIN a genus: snake/snake 100% over 98 aa, the three
bird agents 97-100% over 30 aa, and across genera nothing better than snake vs
marmot at 45% over 22 aa, e=0.002. The cross-lineage similarity in the unmasked
protein is DAg's antisense signal. Both papers say the same in words. Hence:
`Uncharacterized lineage-specific ORF2 protein`, gene_symbol `ORF2`, **17
profiles**, 12 single-sequence, **no copy_num** (uncharacterized AND
lineage-specific AND absent from 5 genomes).

### Cutoff: 100, and it is set by the OTHER module

```
genuine ORF2 calls (36)      min 140.5   5th pct 148.4   median 313.2
the 5 genomes with no ORF2   best ORF2 profile score      0.0
60 HDV genomes, 6 frames,    best ORF2 profile score     86.6   <- the hazard
  circularly doubled           (ORF2.8 86.6, ORF2.1 85.8, ORF2.7 72.7)
```

ORF2 profiles carry DAg's antisense signal, so they light up on HDV's antigen
locus. 100 sits 13 bits above everything foreign and 40 below the lowest genuine
call, and costs nothing -- the call rate is identical at 40 and at 100. Routing
already separates HDV: 40 HDV genomes run end to end routed to `Deltavirus`
40/40 and produced **zero** ORF2 calls. The cutoff no longer depends on that.

### Measured result

| | |
|---|---|
| genomes scored | 41 |
| routed to `Kolmioviridae` | 41 (100%) |
| `DAg` called | 41/41, copy_num 1 |
| `ORF2` called | 36/41 (87.8%) |
| ORF2 on the strand opposite DAg | **36 of 36** |
| ORF2 called twice in one genome | 0 |
| duplicate calls | 0 (ORF2:DAg declared expected architecture) |
| genomes Good | 38/41 (92.7%); the 3 are ambiguous-base flags |
| self-recall | DAg 29/29, ORF2 34/34 |
| `rebuild_pssms.py` | 53 reproduce, 0 stale |
| `check_met_start.py` | every CDS profile begins at its protein's start |
| `qc_cross_feature.py` | 0 Kolmioviridae rows |

### Two stale artifacts found while doing this

- **`Rep-Contigs/` in scratch was the PRE-SPLIT set**: 22 references, five of them
  hepatitis delta virus genomes (12475.2500, .1517, .3659, .2105, .2131, .2664).
  Regenerating the JSON from it re-imports HDV into this module and undoes the
  split. Moved to `Rep-Contigs.PRESPLIT-STALE/` with a README; the live set is
  restored from the shipped module, 17 references, no HDV.
- **`gen_module_json.py` writes `replicon_geometry: linear`** when no `--segments`
  file is passed, silently turning a circular module linear and disabling the
  fold. Now driven from `segments.Kolmioviridae.json`.

### Kit fix made here

`annotation_rarefaction.py` counted the annotator's `circular_contig` bookkeeping
row as a controlled-vocabulary string, so this two-feature module reported three.
It reaches **every circular taxon** -- Kolmioviridae and Hepadnaviridae so far.
Fixed by skipping non-feature rows.

### Still not built

- **The spliced ORF2b isoform.** Junctions appear only in figures of the preprint
  (Figs 1C, 1E, 2C, 2D); bioRxiv rate-limited every attempt to read them. Splice
  junctions come from published coordinates or not at all.
- **The demonstrating genome.** BV-BRC holds scKoV only as a 761 nt partial
  (BR001664). The 16 complete genomes are DDBJ **BR002460-BR002472**, plus
  BR002539 (ichiDV) and BR002541-BR002543 (anatid KoVs) -- all unreleased as of
  27 Nov 2026 (NCBI, ENA and DDBJ all return nothing). Fetch them when they land
  and rebuild; the canary and frog lineages are the two with direct evidence and
  neither contributes a profile today.

---

## OPEN DEFECT — ORF2 is missing from the Kolmioviridae module (26 Nov 2026)

**Status: the Kolmioviridae module is incomplete. Do not treat it as done.**

Raised by Jim from the genome map: one called gene on a 1.7 kb genome, half the
genome unoccupied.

**My error.** I declared ORF2 known-absent across the whole family on the
strength of ONE negative experiment — the snake kolmiovirus ORF2 initiation-codon
knockout (J Gen Virol 10.1099/jgv.0.002327; "no evidence of translation from
ORF2", boa constrictor cells). That result is real and stands FOR THE SNAKE.
Generalising it to every genus was unsupported. I never ran a literature search
for ORF2 being translated anywhere else.

**The literature I missed.** bioRxiv 2025.01.08.631871 v3, *Ambisense
transcription uncovers hidden coding capacity in compact circular RNA satellite
viruses*. ORF2 lies on the GENOMIC-sense strand and is expressed by ambisense
transcription; alternative splicing yields two ORF2 isoforms; the authors report
comparative evidence for ambisense expression in further kolmiovirus lineages.
v1 (Jan 2025) was bioinformatic, v3 adds experimental validation. NOT YET READ IN
FULL — the agent reading it was stopped at the pause. Re-read before building.
No peer-reviewed version found as of 26 Nov 2026.

**Sequence evidence, measured here (all 41 module genomes, ORFs >=70 aa,
circular-aware).** 38 of 41 carry a large genomic-sense ORF. Non-overlapping
with DAg — i.e. a genuine second gene in the empty half:

| genome | DAg (-) | plus-strand ORF | overlap |
|---|---|---|---|
| Chinese fire belly newt virus 1 | 1385-1663 | 211 aa @ 73-705 | 0% |
| rhinotermitid virus 1 (2848056.3) | 555-1034 | 184 aa @ 1-552 | 0% (abuts at 553) |
| Kalotermes flavicollis agent | 959-1378 | 183 aa @ 1-549 | 0% |
| Chusan Island toad virus 1 | 334-654 | 172 aa @ 794-1309 | 0% |
| Epsilon virus (2854842.8) | 526-888 | 172 aa @ 1097-1612 | 0% |
| ray-finned fish virus 1 | 675-1214 | 114 aa @ 1351-86 | 0% |
| Bat deltavirus AYA11 / LMA6 | 1018-1605 | 92 / 86 aa | 0% |
| Neotermes cubanus agent | 1150-16 | 190 aa @ 1-570 | 3% |

Conserved, not chance: all 11 Rodent deltaviruses plus Tome's spiny-rat virus 1
give 196 aa starting at nt 1 — same length, same coordinate, 12 independent
genomes.

**Deltavirus (HDV) looks defensible.** The exemplar's 192 aa plus-strand ORF runs
937-1512, essentially antisense-nested inside the DAg gene — a shadow ORF whose
apparent conservation merely follows DAg's. Across 281 full-length low-N HDV
genomes the longest plus-strand ORF has no consistent length (209 x35, 192 x13,
202 x7, 185 x6, 336 x4 ...). CAUTION: 260 such ORFs cluster into only 20 groups
at 40% id / 60% cov, but that is the same antisense confounder, NOT evidence of
coding. Do not read it as support for an HDV ORF2.

**The demonstrating genome is not in the module.** BV-BRC holds the canary virus
only as `2760986.3 Serinus canaria-associated deltavirus sc-mh1`, 761 nt partial
(BR001664), so the length filter dropped it. Fetch the complete genome from the
preprint's accessions before building ORF2.

**To do, in order.**
1. Re-read bioRxiv 2025.01.08.631871 v3 in full: which viruses, ORF2 length and
   coordinates per virus, what was demonstrated wet-lab vs predicted, the splice
   junctions for the two isoforms, proposed function/motifs, and how it treats
   the snake negative.
2. Decide the feature set: one ORF2, or ORF2 + a spliced isoform (if spliced,
   this needs `special: splice`, with junctions from the paper — never from
   GenBank joins).
3. Decide which lineages get it. The snake is experimentally negative; declaring
   ORF2 for Daletvirus would assert a product shown not to exist. Per-lineage
   call, documented.
4. Pull the complete canary genome and any other ambisense exemplars missing
   from the dump.
5. Rebuild collections, cutoffs, copy_num, the four pages and the genome map.
   `copy_num` must be OMITTED on ORF2 unless it is present in every genome of
   the module, or every other genus is flagged.
6. All new PMIDs go in `PMID_claude_generated` only.

The shipped module, its PSSMs and its ten published pages remain valid for DAg.
Only the claim of completeness is withdrawn.
