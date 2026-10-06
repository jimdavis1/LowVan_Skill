# PICO_NOLDR_VP0 build notes

Built 5 October 2026. **Pasivirus, Avisivirus, Limnipivirus, Potamipivirus,
Aquamavirus** - fish, bird and seal picornaviruses. No leader; VP0 left
uncleaved; one VPg.

The smallest of the four grouped modules: 119 genomes, 68 above 5 kb.
11 features, 28 PSSMs, **14 rep contigs route all 68 = 100%**.

## Aquamavirus fragments dropped

Aquamavirus contributed 2C records of 15 aa and VP3 records of 40 aa against
330 and 230 for every other genus in the group. Those are **fragments, not
short forms**, and left in they would have dragged the 2C window floor to 12
and the VP3 floor to 32 - wide enough to admit any fragment of anything.
Dropped from both collections rather than modelled.

Consequence, recorded rather than hidden: Aquamavirus is represented in the
capsid and polymerase features but NOT in 2C or VP3, because the only records
BV-BRC offers for those are fragments.

Every citation is model-proposed and sits in PMID_claude_generated.

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

**119 sequences added**, every one at >=80% identity to a reference that
annotates the cut itself.

| key | projected | new | was in collection |
|---|---|---|---|
| 2A | 32 | 10 | 23 |
| 2B | 40 | 13 | 32 |
| 2C | 38 | 12 | 33 |
| 3A | 36 | 14 | 25 |
| 3C | 30 | 10 | 29 |
| 3D | 40 | 16 | 30 |
| VP0 | 47 | 15 | 32 |
| VP1 | 37 | 11 | 38 |
| VP3 | 38 | 13 | 32 |
| VPg | 15 | 5 | 18 |

### Features the cleavage-site test refused

Not projected, because the cut cannot be placed with confidence:

- **Avisivirus** `3C P1 E 67% (n=3) P1' G2 S1 NOT CONSERVED (67%)`
- **Avisivirus** `VPG P1 S 67% (n=3) P1' N1 S1 Q1 NOT CONSERVED (67%)`
- **Limnipivirus** `2C P1 E 67% (n=3) P1' G2 D1 NOT CONSERVED (67%)`
- **Limnipivirus** `3A P1 Q 67% (n=3) P1' A2 G1 NOT CONSERVED (67%)`
- **Limnipivirus** `3C P1 E 67% (n=3) P1' G1 A1 S1 NOT CONSERVED (67%)`
- **Limnipivirus** `VP1 P1 Q 67% (n=3) P1' A2 G1 NOT CONSERVED (67%)`
- **Limnipivirus** `VP3 P1 E 67% (n=3) P1' G2 S1 NOT CONSERVED (67%)`
- **Limnipivirus** `VPG P1 Q 50% (n=2) P1' R1 A1 NOT CONSERVED (50%)`
- **Potamipivirus** `2B P1 T 67% (n=3) P1' G2 Q1 NOT CONSERVED (67%)`
- **Potamipivirus** `2C P1 Q 33% (n=3) P1' A1 F1 G1 NOT CONSERVED (33%)`
- **Potamipivirus** `3C P1 Q 33% (n=3) P1' G2 A1 NOT CONSERVED (33%)`
- **Potamipivirus** `3D P1 Q 33% (n=3) P1' G2 V1 NOT CONSERVED (33%)`
- **Potamipivirus** `VP1 P1 Q 67% (n=3) P1' M1 E1 N1 NOT CONSERVED (67%)`
- **Potamipivirus** `VP3 P1 Q 67% (n=3) P1' G2 H1 NOT CONSERVED (67%)`
- **Potamipivirus** `VPG P1 E 33% (n=3) P1' R2 Y1 NOT CONSERVED (33%)`

### What it cost

Coverage over the same 53 exemplars: **379 calls -> 438 (+59)**.
No feature lost a single call.

### Per genus, where the calls actually moved

| genus | feature | before | after | exemplars |
|---|---|---|---|---|
| Aquamavirus | `2C` | 1 | **2** | 2 |
| Aquamavirus | `3C` | 0 | **2** | 2 |
| Aquamavirus | `3D` | 0 | **2** | 2 |
| Aquamavirus | `VP0` | 0 | **2** | 2 |
| Avisivirus | `2B` | 0 | **8** | 8 |
| Avisivirus | `3A` | 0 | **8** | 8 |
| Avisivirus | `3D` | 5 | **8** | 8 |
| Avisivirus | `VP0` | 0 | **8** | 8 |
| Avisivirus | `VP1` | 0 | **8** | 8 |
| Avisivirus | `VP3` | 0 | **8** | 8 |
| Limnipivirus | `2B` | 0 | **4** | 11 |
| Limnipivirus | `2C` | 0 | **4** | 11 |
| Limnipivirus | `3C` | 1 | **0** | 11 |
| Potamipivirus | `VP0` | 1 | **2** | 3 |
| Potamipivirus | `VPg` | 0 | **1** | 3 |

