# Hepatovirus build notes

## Module boundary

One module for the genus: VP4-VP2-VP3-VP1 (VP0 cleaved), 2A, 2B, 2C, 3A,
one VPg, 3C, 3D. No leader. Same organisation as Enterovirus. 25 rep
contigs. Built 2 October 2026; REBUILT 5 October.

1,718 genomes but only 919 unique protein sequences, and 37 species of
which most contribute one or two genomes.

## Rebuild, 5 October 2026 -- the thinnest-trained module of the ten

Before: 42 PSSMs over 12 features from ~340 annotated mature peptides.
VPg had 9 training sequences, 2B had 13, VP4 had 20.

Collections enlarged by projection (cut coordinates carried from annotated
genomes across global alignments; nothing predicted; kept only if the
polyprotein is >=50% identical to its reference (SUPERSEDED: 0.80, see the
6 October section below), the segment is inside the
JSON window and >=35% identical to the reference segment):

3D 40->203, 2C 61->142, VP1 77->138, 2A 31->77, 3C 19->76, VP3 26->74,
3A 23->72, VP2 24->68, 2B 13->49, VP4 20->44, VPG 9->37.

86 PSSMs. Every feature at the documented defaults.

## Measured 5 October 2026

- Coverage: 1,207 genomes >= 5 kb, 116 exemplars, 113 routed, 89 complete
  (was 84). Lowest feature VP4 96/113; everything else >= 99.
- Held-out panel: 39 scored, 36 clean.
- Self-recall 1318/1318.

## Known gaps

- **Three unrouted, all bird-associated**: 2107573.39 and .41 (goose faeces
  Hepatovirus sp.) and 2832213.15 (mute swan hepatovirus 3). Never
  annotated at all.
- **VP4 in the non-primate species.** VP4 is ~23 aa, the shortest capsid
  protein, so divergence costs most there. Misses: 6/12 Hepatovirus sp.,
  2/2 bar-headed goose, 1/1 Loch Leven, 1/1 and 1/1 mute swan - against
  2/59 for hepatitis A virus itself.
- Projection could not reach these lineages either: it needs a reference
  polyprotein >=50% identical (SUPERSEDED: 0.80) and they are further away. The rebuild
  improved the module mostly where it was already adequate.

Hepatitis A virus itself is called completely.

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

**344 sequences added**, every one at >=80% identity to a reference that
annotates the cut itself.

| key | projected | new | was in collection |
|---|---|---|---|
| 2A | 17 | 16 | 31 |
| 2C | 74 | 54 | 61 |
| 3A | 47 | 28 | 23 |
| 3C | 47 | 35 | 19 |
| 3D | 156 | 130 | 40 |
| VP1 | 54 | 33 | 77 |
| VP2 | 33 | 17 | 24 |
| VP3 | 39 | 26 | 26 |
| VP4 | 19 | 5 | 20 |

### Features the cleavage-site test refused

Not projected, because the cut cannot be placed with confidence:

- `2B P1 Q 58% (n=26) P1' G20 S3 A3 NOT CONSERVED (58%)`
- `VPG P1 E 79% (n=29) P1' G28 V1 NOT CONSERVED (79%)`

### What it cost

Coverage over the same 116 exemplars: **1316 calls -> 1303 (-13)**.
Features that lost calls: `2A` -1, `2B` -5, `3A` -6, `VPg` -1.

