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

## What is NOT built — and two that SHOULD be

**Corrected 15 October 2026.** The first version of this section declared all
seven remaining genera unbuilt, on the grounds that their complete-genome
counts fall below the 46 that Seadornavirus cleared. **That was the wrong
test.** Complete-genome count is a proxy for buildability, not a measure of
it: what decides whether a module can be built is whether the sequences
cluster into groups that clear the member floor. Measured properly, two of
them do.

| genus | genomes | complete | ORF1 full-length | clusters @80% | largest | verdict |
|---|---|---|---|---|---|---|
| **Nebovirus** | 285 | 34 | **27** | 3 | **19** | **BUILD** |
| **Bavovirus** | 56 | 28 | **20** | 2 | **12** | **BUILD** |
| Recovirus | 72 | 5 | not measured | — | — | below floor |
| Valovirus | 9 | 5 | — | — | — | below floor |
| Salovirus | 6 | 5 | — | — | — | below floor |
| Nacovirus | 4 | 4 | — | — | — | below floor |
| Minovirus | 2 | 2 | — | — | — | below floor |
| (no genus on the record) | 634 | — | — | — | — | triage first |

Clusters of 19, 19, 12 and 14 clear `-m 5` outright, with no concession.

### These 341 genomes are currently annotated with NOTHING

61 near-complete Nebovirus and Bavovirus genomes were run through the shipped
four-module repo. **Zero features were called on any of them** — not a poor
annotation, no annotation at all. That is the strongest case for building a
module that this project has measured.

### Nebovirus is fully worked out and ready to build

Two ORFs, the Sapovirus/Lagovirus layout. ORF1 is **2,210 aa** in the
literature and the deposits give a median of exactly 2,210; ORF2 is 225 aa and
the deposits give exactly 225.

The published map (Newbury-1, PMC7111791 and J Virol 76:10089) is
**NH2-P34-P37(2C)-P30-P7.6(VPg)-P20(3C)-P53(3D)-P57(VP1)-COOH**, with cleavage
at `[EQ]↓[GA]`. The deposits do NOT name the mature peptides, so the vocabulary
cross-check that pinned Lagovirus's two undocumented boundaries is unavailable
here — but three independent lines converge anyway:

- **catalytic anchors land inside their predicted products**: `GPPGHGKS`
  (2C Walker A) uniquely at 456, `GYCG` (3C) at 1101, `GLPSG` at 1440 and
  `YGDD` at 1485 (3D);
- **conserved cleavage sites sit where the kDa map predicts**, every one of the
  six junctions, with a majority of the 27 sequences agreeing:

```
P34/P37   302  PKQE|GG   17/27      VPg/3C   990  LGRE|GN   20/27
P37/P30   645  PDRE|GP   19/27      3C/3D   1174  MKEE|AA   13/27
P30/VPg   925  NDRQ|GW   19/27      3D/VP1  1659  LIRE|GN   20/27
```

- the shared **P2 arginine** across four of the six sites is a real signal, and
  makes per-junction motifs specific.

So Nebovirus is a full module — POLY, seven ORF1 products derived at measured
motifs, VP2 from its own 28 deposits — using the technique already proven
twice in this family (Vesivirus VP1, the Norovirus conserved cores).

### Bavovirus is probably the same shape, and is not yet worked out

Chicken calicivirus. 20 full-length ORF1 in 2 clusters, 20 VP2 in 2 clusters,
so it builds. Its literature has **not** been searched, and it is far less
studied than Nebovirus, so its mature-peptide set may have to be scoped to
whatever is documented — or projected from Nebovirus, which would need the 80%
identity floor checked first.

## Open items, in the order I would take them

1. **Triage the 634 unclassified genomes.** Largest unbuilt population in the
   family and the only one that might be free.
2. **Norovirus ORF4 / VF1 is not declared.** 39 features carry the name. It is
   murine norovirus virulence factor 1, a real protein and not a synonym for
   anything the module models. It was never evaluated as a feature; it should
   be, and either declared or tabled with a measurement.
3. **Vesivirus ORF1 products are FCV-only.** 210 of 212 feline genomes carry
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

## Settled, not open

**Sapovirus's 147 bindable sequences: not worth binding.** Decided by the
curator, 15 October 2026. The measurement stands — 22 annotation strings map
unambiguously onto a declared feature, and binding every one of them while
keeping only what clears the length window adds 147 sequences: POLY +4.2%,
everything else under 1.5%. Against that, re-running `build_collections.py`
discards the homology rescue and costs a full rebuild of a module measuring
411/411 routed. **The gain does not pay for the rebuild.**

The spellings are still recorded in `UNASSIGNED_DECISION.md`, and they were
worth recording: Lagovirus and Vesivirus bound their equivalents correctly on
the first pass because of them. This closes the question for Sapovirus only —
a future rebuild of that module for some other reason should pick them up on
the way past.

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
