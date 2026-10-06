# PICO_LDR_VP0 build notes

Built 5 October 2026. **Sicinivirus, Salivirus, Oscivirus, Passerivirus,
Shanbavirus.** Leader present; VP0 left uncleaved (so VP0-VP3-VP1, three
capsid proteins not four); one VPg. Salivirus causes human gastroenteritis;
the rest are bird viruses.

281 genomes (219 >= 5 kb). 13 features, 67 PSSMs, 25 rep contigs routing
196/218 = **89.9%**.

## One slot held two proteins

| slot | split off | lengths |
|---|---|---|
| Ldr | **Sicinivirus** -> `Ldr_LONG` | 462 aa vs 64-174 for the other four |

`Ldr` is "Leader protein". `Ldr_LONG` is declared **"Uncharacterized
lineage-specific protein"** with gene_symbol `Ldr-long`: a PubMed search found
no functional description for the Sicinivirus leader, and the skill forbids
inventing a designation. The protein is modelled and called; it is simply not
named beyond its lineage. Revisit if a description appears.

## Known gaps

- 22 genomes unrouted (22 clusters), a long tail of singletons.
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

**177 sequences added**, every one at >=80% identity to a reference that
annotates the cut itself.

| key | projected | new | was in collection |
|---|---|---|---|
| 2A | 7 | 0 | 24 |
| 2B | 6 | 0 | 23 |
| 2C | 50 | 35 | 27 |
| 3A | 45 | 31 | 26 |
| 3C | 7 | 0 | 26 |
| 3D | 51 | 37 | 38 |
| Ldr | 30 | 19 | 13 |
| Ldr_LONG | 37 | 33 | 8 |
| VP0 | 8 | 0 | 27 |
| VP1 | 22 | 7 | 27 |
| VP3 | 8 | 0 | 28 |
| VPg | 27 | 15 | 19 |

### Features the cleavage-site test refused

Not projected, because the cut cannot be placed with confidence:

- **Salivirus** `2A P1 Y 50% (n=8) P1' S4 G2 L1 NOT CONSERVED (50%)`
- **Salivirus** `2B P1 Q 86% (n=7) P1' G6 D1 NOT CONSERVED (86%)`
- **Salivirus** `2C P1 Q 86% (n=7) P1' G6 D1 NOT CONSERVED (86%)`
- **Salivirus** `3A P1 Q 83% (n=6) P1' G5 A1 NOT CONSERVED (83%)`
- **Salivirus** `3C P1 Q 71% (n=7) P1' G5 P1 A1 NOT CONSERVED (71%)`
- **Salivirus** `3D P1 Q 80% (n=5) P1' S4 H1 NOT CONSERVED (80%)`
- **Salivirus** `VP0 P1 Q 67% (n=6) P1' G4 M1 I1 NOT CONSERVED (67%)`
- **Salivirus** `VP1 P1 Q 86% (n=7) P1' S6 D1 NOT CONSERVED (86%)`
- **Salivirus** `VP3 P1 P 62% (n=8) P1' Q5 H2 K1 NOT CONSERVED (62%)`
- **Salivirus** `VPG P1 Q 80% (n=5) P1' G4 N1 NOT CONSERVED (80%)`
- **Sicinivirus** `2A P1 Q 60% (n=5) P1' G4 S1 NOT CONSERVED (60%)`
- **Sicinivirus** `2B P1 Q 86% (n=7) P1' A7 NOT CONSERVED (86%)`
- **Sicinivirus** `3C P1 Q 71% (n=7) P1' I4 C2 L1 NOT CONSERVED (71%)`
- **Sicinivirus** `VP0 P1 Q 57% (n=7) P1' T4 G3 NOT CONSERVED (57%)`
- **Sicinivirus** `VP3 P1 Q 86% (n=7) P1' G4 A2 H1 NOT CONSERVED (86%)`

### What it cost

Coverage over the same 136 exemplars: **1121 calls -> 1122 (+1)**.
Features that lost calls: `3A` -10.

### Per genus, where the calls actually moved

| genus | feature | before | after | exemplars |
|---|---|---|---|---|
| Oscivirus | `Ldr` | 4 | **2** | 4 |
| Salivirus | `3A` | 13 | **4** | 13 |
| Shanbavirus | `3A` | 5 | **4** | 53 |
| Shanbavirus | `Ldr` | 0 | **2** | 53 |
| Shanbavirus | `VPg` | 21 | **25** | 53 |

**Salivirus `3A` 13 -> 4 is a re-clustering effect, not bad data.** The 3A
collection went 26 -> 57 sequences, all of the new ones at >=80% identity to a
reference that annotates the cut. With the larger set the main pass clusters
differently -- 2 profiles and 20 leftovers where it previously covered
Salivirus -- and the profile that used to call those nine genomes no longer
exists. The length window is not the cause: it was held at 120-202 by the
"a window may widen but never narrow" rule, and widening it alone recovers
only three of the lost calls.

The module is net +1, and VPg gains 11. Recovering Salivirus 3A means splitting
that feature by genus or forcing the leftover pass to cover it, and is worth
doing before this module is considered finished.
