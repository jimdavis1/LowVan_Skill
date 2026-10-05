# Kobuvirus build notes

## Module boundary

One module for the genus: Ldr, VP0 (uncleaved), VP3, VP1, 2A, 2B, 2C, 3A,
one VPg, 3C, 3D. VP0 uncleaved groups this genus with Parechovirus and
Avihepatovirus rather than Enterovirus. 23 rep contigs. Built 2 October
2026; leader renamed L -> Ldr 3 October; REBUILT 5 October.

Most species-rich module in the family: 45 species, 956 genomes.

## Rebuild, 5 October 2026

Collections enlarged by projection (as for the other genera; nothing
predicted): VP1 222->408, 3D 222->403, VP0 182->351, VP3 116->275,
3C 100->258, 2A 96->249, 2C 103->244, Ldr ->268, 3A 93->224, 2B 85->209,
VPG 43->92. 97 -> 115 PSSMs, all at the documented defaults.

**The rebuild did not change the coverage numbers**: 278 complete genomes
before and after, Ldr 287/313 both times. This is informative. Projection
converts UNANNOTATED genomes into training data by carrying cuts from a
close relative, so it helps enormously where a feature is annotated in
almost no genome (Sapelovirus: 11 leaders -> 6/209 calls became 170/208)
and very little where the feature is already represented and the gap is a
lineage no reference comes within 50% identity of. Kobuvirus is the second
case.

## Measured 5 October 2026

- Coverage: 515 genomes >= 5 kb, 313 exemplars, 313 routed (100%), 278
  complete. Lowest feature Ldr 287/313; everything else >= 305.
- Held-out panel: 40 genomes, 32 clean. The eight flags are ambiguous-base
  contigs and doubled polyprotein calls - deposit quality, and more common
  here because much of this genus comes from metagenomic surveys.
- Self-recall 3331/3331.

## Known gap

The leader in the one- and two-genome species: 3/3 K. ecuni, 3/3
K. femyomini, 2/2 grey squirrel kobuvirus, 5/15 Kobuvirus sp. - against
2/146 K. cebes and 2/49 K. aichi. Every other protein is called on those
same genomes, so the module reaches them; the leader is simply the least
conserved feature. A sampling limit in BV-BRC, demonstrated by the rebuild
leaving the number unmoved.
