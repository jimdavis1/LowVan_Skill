# PICO_LDR_VP4 build notes

Built 5 October 2026. Six genera sharing one genome organisation:
**Bopivirus, Erbovirus, Mischivirus, Hunnivirus, Tottorivirus, Anativirus.**
Leader present; capsid cleaved to VP4+VP2+VP3+VP1; one VPg.

396 genomes (254 >= 5 kb). No genus here holds enough sequence to model
alone - Tottorivirus has 10 genomes, Anativirus 38 - so they are grouped and
their diversity becomes extra profiles inside each feature.

15 features, 95 PSSMs, 25 rep contigs routing 218/250 = **87.2%**
(28 clusters holding 32 genomes uncovered).

## Two slots held two different proteins

A slot name is a homology group INSIDE a genus and not across one. Checking
per-genus medians inside the group before building caught two cases where a
single feature would have blurred unrelated proteins:

| slot | split off | lengths |
|---|---|---|
| 2A | **Bopivirus** -> `2A_LONG` | 212 aa vs 12-55 for the StopGo peptides of the rest |
| Ldr | **Erbovirus** -> `Ldr_PRO` | 219 aa vs 51-84 for the rest |

`2A` keeps the vocabulary string "Protein 2A, StopGo peptide" and `2A_LONG`
takes "Protein 2A". `Ldr` is "Leader protein"; `Ldr_PRO` is "Leader protease
(Lpro)" - ERBV's L protein has C-terminal processing activity (Hinton 2002
J Gen Virol 83:3111, PMID 12466488, model-proposed).

Two features in one module must never share an annotation string: the quality
tool keys %essential by the string, so the second read silently overwrites the
first's length window. make_group_json.py refuses to emit a module with a
duplicate.

## Known gaps

- 87.2% routing is the lowest of the four grouped modules and the 32 unrouted
  genomes are a long tail of singletons across six genera.
- Every citation is model-proposed and sits in PMID_claude_generated.

## Mature-peptide projection rebuilt at a defensible identity floor (6 October 2026)

The first projection ran at `--min-id 0.50 --min-seg-id 0.35`. Both are
far too low: below ~80% the global alignment stops being trustworthy
residue by residue and a transferred cut lands wherever the gaps fell.
754 of 13,884 projected sequences family-wide had come across below 80%.

Rebuilt from the pre-projection collections at **`--min-id 0.80`,
`--min-seg-id 0.80`**, plus a cleavage-site test: for every feature the
P1 residue is measured across all the references that place it, and the
feature is projected only if P1 is >=90% conserved **and every projected
cut lands on that residue**. Conservation is the criterion, not the
identity of the residue, because the polyprotein carries three kinds of
cut -- 3C protease (Q|G canonical, also Q|S, Q|R, E|G, E|Q; the S1 pocket
reading P1 is the most conserved part of the protease across the family,
Meng 2022 *J Virol* 96:e0073622, PMID 35727031), the StopGo skip at the
end of the D(V/I)ExNPG|P motif, and the autocatalytic VP0 -> VP4+VP2
maturation cut.

**202 sequences added**, every one at >=80% identity to a reference that
annotates the cut itself.

| key | projected | new | was in collection |
|---|---|---|---|
| 2A | 5 | 2 | 16 |
| 2B | 52 | 22 | 41 |
| 2C | 49 | 21 | 35 |
| 3A | 40 | 17 | 28 |
| 3C | 44 | 27 | 36 |
| 3D | 63 | 27 | 62 |
| Ldr | 59 | 26 | 33 |
| Ldr_PRO | 6 | 0 | 6 |
| VP1 | 42 | 14 | 57 |
| VP2 | 46 | 16 | 49 |
| VP3 | 41 | 14 | 48 |
| VP4 | 18 | 6 | 36 |
| VPg | 26 | 10 | 27 |

### Features the cleavage-site test refused

Not projected, because the cut cannot be placed with confidence:

- **Anativirus** `VP1 P1 K 50% (n=2) P1' G2 NOT CONSERVED (50%)`
- **Bopivirus** `2A_LONG P1 D 88% (n=25) P1' G23 L1 D1 NOT CONSERVED (88%)`
- **Bopivirus** `VP4 P1 T 88% (n=25) P1' G24 M1 NOT CONSERVED (88%)`
- **Erbovirus** `2A P1 S 75% (n=4) P1' N4 NOT CONSERVED (75%)`
- **Erbovirus** `3C P1 E 86% (n=7) P1' N6 E1 NOT CONSERVED (86%)`
- **Erbovirus** `VPG P1 S 86% (n=7) P1' R7 NOT CONSERVED (86%)`
- **Mischivirus** `2A P1 D 50% (n=8) P1' G6 N2 NOT CONSERVED (50%)`
- **Mischivirus** `2B P1 G 50% (n=8) P1' P4 G4 NOT CONSERVED (50%)`
- **Mischivirus** `2C P1 Q 86% (n=7) P1' G5 P2 NOT CONSERVED (86%)`
- **Mischivirus** `3A P1 Q 80% (n=5) P1' G4 A1 NOT CONSERVED (80%)`
- **Mischivirus** `3C P1 Q 88% (n=8) P1' S4 G2 A2 NOT CONSERVED (88%)`
- **Mischivirus** `3D P1 Q 50% (n=8) P1' G8 NOT CONSERVED (50%)`
- **Mischivirus** `VP1 P1 Q 75% (n=8) P1' G6 T2 NOT CONSERVED (75%)`
- **Mischivirus** `VP2 P1 K 67% (n=6) P1' S4 D2 NOT CONSERVED (67%)`
- **Mischivirus** `VP3 P1 Q 88% (n=8) P1' G4 M2 Q1 NOT CONSERVED (88%)`
- **Mischivirus** `VP4 P1 C 57% (n=7) P1' G7 NOT CONSERVED (57%)`
- **Mischivirus** `VPG P1 E 57% (n=7) P1' G6 S1 NOT CONSERVED (57%)`
- **Tottorivirus** `2A P1 T 50% (n=2) P1' S1 G1 NOT CONSERVED (50%)`
- **Tottorivirus** `3C P1 N 50% (n=2) P1' G2 NOT CONSERVED (50%)`

### What it cost

Coverage over the same 201 exemplars: **2072 calls -> 2381 (+309)**.
No feature lost a single call.

### Per genus, where the calls actually moved

| genus | feature | before | after | exemplars |
|---|---|---|---|---|
| Hunnivirus | `2B` | 0 | **76** | 76 |
| Hunnivirus | `3A` | 0 | **71** | 76 |
| Hunnivirus | `Ldr` | 0 | **54** | 76 |
| Hunnivirus | `VP1` | 0 | **75** | 76 |
| Mischivirus | `2C` | 36 | **38** | 42 |
| Mischivirus | `Ldr` | 12 | **21** | 42 |
| Mischivirus | `VP2` | 34 | **35** | 42 |
| Tottorivirus | `2B` | 0 | **6** | 6 |
| Tottorivirus | `VP1` | 0 | **6** | 6 |

