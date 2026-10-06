# Senecavirus build notes

## Module boundary

One module for the genus: Ldr, VP4-VP2-VP3-VP1, 2A (9 aa StopGo), 2B, 2C,
3A, one VPg, 3C, 3D. 3 rep contigs route all 993 genomes.
Built 2 October 2026; REBUILT 5 October.

993 genomes collapse to **10 exemplars** at 95% nucleotide identity. Every
coverage percentage here has a denominator of 10.

## Rebuild, 5 October 2026

Collections enlarged by projection (cut coordinates carried from annotated
genomes across global alignments; nothing predicted):
3D 35->371, 3A 29->215, VP3 28->230, VP2 13->189, 2B 13->118, 2C 21->146,
3C 19->166, VP1 25->189, Ldr 13->98, VP4 4->31, VPG 3->23.

2A did not grow: it is 9 residues (SGDVETNPG, the StopGo NPG motif) and the
genus holds 8 distinct ones.

## DEFECT FOUND AND FIXED: 2A shipped with no profile

The first rebuild used the thin-collection default `-m 2` at `-mi 0.8`.
Eight 9-residue peptides produce no cluster at that identity, so no profile
was built - but 2A stayed declared in the JSON. A declared feature with no
profile is never called and is silent at runtime; the annotator reports it
missing on every genome.

`install_module.py --check` DID report it (`1 problem`), but the wrapper
script grepped only for success lines so the warning was never seen. The
2 October build had used `-m 1 -mi 0.6` and been correct.

Both halves fixed:
 - rebuild_genus now picks `-m 1` for any collection with <10 sequences or a
   median length <20 residues. **The MEMBER floor only - never -mi.** The
   first fix also dropped the identity floor to 0.6, which was wrong and
   unnecessary: 2A builds at the documented `-mi 0.8` with `-m 1` alone.
   Lowering -mi redefines "same protein" and merges unrelated sequences
 - a hard gate refuses to ship if any non-special feature has no PSSM

2A: 0 calls -> 10/10.

Departures from the documented defaults: POLY `-m 2`, 2A `-m 1`. 22 PSSMs.

## Measured 5 October 2026

- Coverage: 956 genomes >= 5 kb, 10 exemplars, 10 routed, 6 complete.
  VPG 7/10, Ldr 8/10, 3A 8/10; everything else 10/10.
- Held-out panel: 40 genomes, 39 clean; the one flag is a missing POLY.
- Self-recall 2511/2511.

## Known limitation

Not a gap but a denominator. Ten distinct genome types is all BV-BRC holds,
so the full-taxon run is barely larger than the held-out panel. The module
is sound for Seneca Valley virus as deposited; nothing here says how it
behaves on a divergent senecavirus, because none exists in the data.

The 9 aa 2A works because this genus is homogeneous enough that nine
residues are near-invariant. That would not transfer: the same StopGo slot
in Teschovirus and Aphthovirus holds 21 and 18 aa peptides.
