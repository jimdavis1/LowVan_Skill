# Sapelovirus build notes

## Module boundary

One module for the genus: Ldr, VP4-VP2-VP3-VP1, 2A (~226 aa, a full protein
not a StopGo peptide), 2B, 2C, 3A, one VPg, 3C, 3D. 25 rep contigs.
Built 2 October 2026; leader renamed L -> Ldr 3 October; REBUILT 5 October.

## Rebuild, 5 October 2026 -- the leader

First full-taxon run: Ldr called on 6 of 209 exemplars, and only 5 exemplars
carried all thirteen features. The Ldr collection held 11 sequences.

Projection (coordinates carried from annotated genomes across a global
alignment) could not help: not one of the 11 genomes with a leader annotation
also places VP4, so no reference positioned both.

**The leader was DERIVED instead.** In an L+ picornavirus the polyprotein
starts with the leader, so Ldr = polyprotein start -> start of VP4. Of 170
genomes placing both POLY and VP4, **166 give exactly 84 residues and the
junction reads ...NKPQ|GAYN - a canonical 3C Q|G site.** The boundary lands
on the protease's own recognition signature in 98% of genomes, independently
of anything we imposed.

(The 11 source-annotated leaders are 88 aa, 4 longer, and come from a
different species group; both forms are kept in the collection.)

Now a general `--first-product KEY:NEXTKEY` option in project_matpeptides.py.

| | before | after |
|---|---|---|
| Ldr collection | 11 | 75 |
| Ldr called | 6/209 | 170/208 |
| exemplars with all 13 features | 5/209 | 167/208 |

63 PSSMs. No feature needed a departure from the documented defaults.

## Measured 5 October 2026

- Coverage: 350 genomes >= 5 kb, 209 exemplars, 208 routed, 167 complete.
- Held-out panel: 40 genomes, 34 clean. Three of the six flagged carry
  contigs >1% ambiguous bases; one has POLY called twice. Lowest clean rate
  of the ten modules, tracking deposit quality rather than the module.
- Self-recall 1862/1862.

## Known gap

38 exemplars still miss Ldr and they are entirely divergent host lineages:
4/4 Pteropodidae bat, 2/2 bat, 2/2 and 2/2 mute swan feces associated,
2/2 coypu sapelovirus 2, 12/21 `Sapelovirus sp.` All route and get their
other twelve features. The derivation needs a reference polyprotein >=50%
identical and these are further away. Needs annotated genomes from those
lineages, not a threshold change.
