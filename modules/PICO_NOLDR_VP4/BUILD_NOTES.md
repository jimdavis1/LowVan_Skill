# PICO_NOLDR_VP4 build notes

Built 5 October 2026. **Cosavirus, Dicipivirus, Tremovirus.** No leader;
capsid cleaved to VP4+VP2+VP3+VP1; one VPg. Cosavirus is a human enteric
virus, Tremovirus is avian encephalomyelitis virus.

181 genomes (134 >= 5 kb). 13 features, 68 PSSMs, 25 rep contigs routing
133/134 = **99.3%**.

## One slot held two proteins

| slot | split off | lengths |
|---|---|---|
| 2A | **Tremovirus, Dicipivirus** -> `2A_LONG` | 141-164 aa vs 31 for Cosavirus |

Cosavirus 2A keeps "Protein 2A, StopGo peptide"; the long form takes
"Protein 2A".

## Known limitation

Dicipivirus has a **dicistronic** genome - two ORFs with an intergenic IRES
rather than one polyprotein. This module represents its polyprotein products
only and does not model that arrangement. The genus was grouped here on its
capsid and replicase layout; whether it deserves its own module on genome
organisation grounds is an open question and 30 genomes is thin either way.

Every citation is model-proposed and sits in PMID_claude_generated.

## Tremovirus is mis-partitioned, and `VP4` is the symptom

`VP4` is called on 54/54 Cosavirus and 16/19 Dicipivirus exemplars and on
**0/10 Tremovirus**, unchanged by the projection rebuild. The cause is a size
difference, not a profile gap:

| genus | VP4 median | n |
|---|---|---|
| Cosavirus | 68 | 25 |
| Dicipivirus | 44 | - |
| Tremovirus | **20** | 4 (range 19-23) |

The feature declares 54-82 with `bit_cutoff` 31, computed from the pooled
median of 68. A 20-residue Tremovirus VP4 fails the window, and would fail the
cutoff even if the window were widened -- a profile that short frequently
cannot score its own training sequences above 31 bits.

**A second VP4 feature is not the fix.** The two forms are the same protein, so
it would need a second annotation string for `Capsid protein VP4 (1A)`, and one
annotation string means one feature. The real statement is that a ~20 aa VP4 is
a different genome organisation: it is the *hepatovirus* form (Hepatovirus VP4
median 23, declared 18-28 in its own module), not the cosavirus one. By the
rule that a module is defined by genome organisation rather than phylogeny,
Tremovirus does not belong in this module.

Moving it is a partitioning change and is left for a curator. Tremovirus routes
here and gets its other twelve features; only VP4 is lost.

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

**68 sequences added**, every one at >=80% identity to a reference that
annotates the cut itself.

| key | projected | new | was in collection |
|---|---|---|---|
| 2B | 31 | 9 | 39 |
| 2C | 29 | 12 | 45 |
| 3A | 28 | 6 | 42 |
| 3C | 45 | 15 | 44 |
| 3D | 32 | 9 | 58 |
| VP2 | 31 | 4 | 59 |
| VP3 | 44 | 9 | 60 |
| VP4 | 29 | 4 | 37 |

### Features the cleavage-site test refused

Not projected, because the cut cannot be placed with confidence:

- **Cosavirus** `2A P1 G 30% (n=27) P1' S10 N4 I4 NOT CONSERVED (30%)`
- **Cosavirus** `VP1 P1 Q 70% (n=27) P1' G9 T5 N5 NOT CONSERVED (70%)`
- **Cosavirus** `VPG P1 A 56% (n=27) P1' G27 NOT CONSERVED (56%)`
- **Tremovirus** `2A_LONG P1 D 50% (n=2) P1' E1 A1 NOT CONSERVED (50%)`
- **Tremovirus** `2C P1 Q 86% (n=7) P1' G5 S1 K1 NOT CONSERVED (86%)`
- **Tremovirus** `3A P1 Q 50% (n=2) P1' M1 G1 NOT CONSERVED (50%)`
- **Tremovirus** `3D P1 Q 86% (n=7) P1' C5 R2 NOT CONSERVED (86%)`
- **Tremovirus** `VP1 P1 Q 86% (n=7) P1' G5 A1 T1 NOT CONSERVED (86%)`
- **Tremovirus** `VP2 P1 N 71% (n=7) P1' D6 N1 NOT CONSERVED (71%)`
- **Tremovirus** `VPG P1 Q 86% (n=7) P1' S6 G1 NOT CONSERVED (86%)`

### What it cost

Coverage over the same 83 exemplars: **901 calls -> 901 (+0)**.
No feature lost a single call.

