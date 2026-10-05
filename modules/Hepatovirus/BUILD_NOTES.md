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
polyprotein is >=50% identical to its reference, the segment is inside the
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
  polyprotein >=50% identical and they are further away. The rebuild
  improved the module mostly where it was already adequate.

Hepatitis A virus itself is called completely.
