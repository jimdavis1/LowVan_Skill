# Caliciviridae — four modules, complete

Built 13–15 October 2026. Four modules, one per genus, partitioned on genome
organisation rather than phylogeny (`PARTITIONING.md`). Cleavage maps taken
from the literature before any profile was built (`CLEAVAGE_SITES.md`).

| module | features | PSSMs | refs | routed | panel | tiling |
|---|---|---|---|---|---|---|
| Norovirus | 9 | 176 | 25 | 411/411 | 39/40 | exact |
| Sapovirus | 9 | 125 | 25 | 220/223 | 35/35 | exact |
| Lagovirus | 10 | 15 | 2 | **69/69** | **40/40** | exact |
| Vesivirus | 11 | 32 | 19 | **276/276** | 39/40 | exact |

Self-recall is 100% on every feature of every module except Norovirus VP2
(2,794 of 2,795). String collapse 17.6× / 12.2× / 7.2× / 8.7×; vocabulary
collapse at equal depth 6.1× / 6.4× / 3.2× / 4.9×.

**"Tiling exact" means the products sum to their precursor minus one codon** —
the stop, which a CDS carries in its coordinates and a mature peptide never
does. Vesivirus ORF1 reaches p05 = p95 = −1 codon: every genome, zero spread.

## Why four modules

    Norovirus   ORF1 (six NS)                  | ORF2 VP1 | ORF3 VP2   3 ORFs
    Vesivirus   ORF1 (six NS)                  | ORF2 LC+VP1 | ORF3 VP2 3 ORFs
    Sapovirus   ORF1 (six NS + VP1 IN FRAME)   | ORF2 VP2             2 ORFs
    Lagovirus   ORF1 (six NS + VP1 IN FRAME)   | ORF2 VP2             2 ORFs

Genome organisation, not the tree. Best blastp identity between genera is
30–47% and mostly low 30s, far under the 80% this project requires before
anything is shared across an alignment. Lagovirus cleaves Pro from Pol where
Sapovirus and Vesivirus keep a single Pro-Pol, so it declares one more ORF1
product; Vesivirus cleaves ORF2 into a leader plus the mature capsid, which no
other genus does, so it declares eleven features — the most in the family.

## The vocabulary cross-check, which earned its place

For every genus, the **modal deposited length** of the records carrying each
mature-peptide name was compared against the literature cleavage map. It
agreed every time, and it did three further jobs no one asked of it:

- **Filled two undocumented Lagovirus boundaries.** p23/p37 and p29/VPg are
  not in the literature. The medians fixed them, and the same measurement
  reproduced all five documented RHDV sites without being given them. The
  eight products sum to 2,343 against an ORF1 of 2,344.
- **Corrected two errors in this project's own notes.** The Vesivirus ORF2 cut
  was recorded as E157/S158 while the same line described a 14.2 kDa leader
  (14.2 kDa is ~124 residues); and P1′ is **A** in 275 of 323 precursors, not
  S in 39. The cut is **E124/A125**.
- **Settled the U-number trap by measurement.** Sapovirus `ORF2` and
  `ORF3 protein` both sit at a median of 166 residues: submitters disagree
  about whether the VP1-containing frame counts as one ORF or two, and only
  the length distinguishes them. Both name VP2.

## What is NOT built, and the numbers that justify it

Seven genera and the records carrying no genus at all. The floor is what
Seadornavirus cleared: **46 complete genomes, 1,202 features, 696 sequences in
collections.** Nothing below is close.

| genus | genomes | complete ≥6.5 kb | features | verdict |
|---|---|---|---|---|
| Nebovirus | 285 | **34** | 361 | closest; still under the floor on both counts |
| Bavovirus | 56 | **28** | 87 | a quarter of Seadornavirus's features |
| Recovirus | 72 | 5 | 103 | 5 complete genomes |
| Valovirus | 9 | 5 | 19 | — |
| Salovirus | 6 | 5 | 12 | — |
| Nacovirus | 4 | 4 | 10 | — |
| Minovirus | 2 | 2 | 8 | — |
| **(no genus on the record)** | **634** | not measured | 1,066 | **worth triaging — see below** |

**1,068 unbuilt genomes against 86,614 built**, so the family is 98.8% covered
by genome count. Nebovirus is the one to revisit first if it grows: it shares
the Sapovirus/Lagovirus two-ORF layout, so it would join that partition rather
than need a new one.

**The 634 records with no genus are the real loose end.** They carry 1,066
protein features — more than Nebovirus, Recovirus and Bavovirus combined — and
nothing has been done to find out what they are. Some will be assignable by
homology to a module that already exists, which would cost nothing but a
rescue pass; the rest should be declared.

## Open items, in the order I would take them

1. **Triage the 634 unclassified genomes.** Largest unbuilt population in the
   family and the only one that might be free.
2. **Norovirus ORF4 / VF1 is not declared.** 39 features carry the name. It is
   murine norovirus virulence factor 1, a real protein and not a synonym for
   anything the module models. It was never evaluated as a feature; it should
   be, and either declared or tabled with a measurement.
3. **Sapovirus: 147 bindable sequences tabled.** 22 annotation strings map
   unambiguously onto a declared feature but are not matched by the rules —
   `VP-1` hyphenated, bare `ORF1`, `small basic protein`, the typo
   `major capid protein`. Binding all of them and keeping only what clears the
   length window adds 147 sequences: POLY +4.2%, everything else under 1.5%.
   Measured and tabled rather than guessed at, because re-running
   `build_collections.py` discards the homology rescue. `UNASSIGNED_DECISION.md`
   has the spellings.
4. **Vesivirus ORF1 products are FCV-only.** 210 of 212 feline genomes carry
   all six; **none** of the 58 non-FCV genomes carries any. The product names
   (`p5.6`, `p32`, `p39`, `p30`, `Pro-Pol`) are FCV nomenclature, so the
   collections are FCV sequences and the profiles do not reach a canine or mink
   vesivirus. All seven affected features are declared **non-essential**, so
   those genomes are not flagged poor. Closing it needs collections for the
   non-FCV products — either their deposits under other names, or derivation
   from the polyprotein by motif as VP1 was. It does **not** close by lowering
   a cutoff: self-recall is already 100%, so the profiles are not failing to
   find what they were built from; they were built from the wrong half of the
   genus.

## What the four builds cost, and where

The collections and the literature maps were right on the first pass every
time. **Every defect came from tooling**, and each produced a permanent check
that then found something the previous checks could not:

| defect | found on | caught by | now |
|---|---|---|---|
| extension flags discarded by the generator | Norovirus | measuring annotations | `gen_module_json` honours them; mat_peptide defaults to no Met scan |
| N-terminal trimming of CDS profiles | Norovirus | `check_met_start` (new) | gate; `filter_nterm` |
| bit cutoffs clipping their own training data | Sapovirus | self-recall | `derive_cutoffs` (new), margin 0.9 |
| the gate skipped *leading* mature peptides | Lagovirus | extending that gate | tests `upstream_ext: 1` mat_peptides |
| `filter_nterm` believed two over-annotated records | Vesivirus | regression test | extensions must be corroborated |
| self-recall measured with the wrong search | Vesivirus | NS1/NS5 reading 66%/36% | `collect_registry` uses `-comp_based_stats 0` |
| registry pages printing a Python dict | all four | reading a published page | generator renders mappings; lint catches it |

The last two are worth remembering together: one was caught only because a
**published page was read before republishing it**, and the other only because
the module with the shortest features finally made a long-standing measurement
error visible.
