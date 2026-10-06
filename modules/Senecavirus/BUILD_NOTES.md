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

### 2A `-m 1` is a deliberate, approved exception

The documented pass (`-m 5`) and `-m 2` both produce no cluster and no output
directory for a set of eight 9-residue peptides, so the leftover pass - the
sanctioned route to a 2-member profile - has no `Leftover_Seqs.aa` to read.
`-m 1` yields one 2-member cluster and that profile calls 2A on 10/10
exemplars. Approved 5 Oct 2026 and recorded in the feature's BUILD_PARAMS.

**The identity floor is NOT lowered**: `-mi` stays at the documented 0.8.
Two-member clusters are fine; a member floor of 1 is not, and this is the
only feature in the family that uses one. Do not copy it anywhere that can
be built by a sanctioned route.

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

**1559 sequences added**, every one at >=80% identity to a reference that
annotates the cut itself.

| key | projected | new | was in collection |
|---|---|---|---|
| 2A | 3 | 1 | 4 |
| 2B | 113 | 104 | 13 |
| 2C | 136 | 124 | 21 |
| 3A | 205 | 185 | 29 |
| 3C | 160 | 145 | 19 |
| 3D | 358 | 336 | 35 |
| Ldr | 93 | 85 | 13 |
| VP1 | 181 | 162 | 25 |
| VP2 | 182 | 174 | 13 |
| VP3 | 219 | 199 | 28 |
| VP4 | 28 | 25 | 4 |
| VPg | 21 | 19 | 3 |

### What it cost

Coverage over the same 10 exemplars: **126 calls -> 126 (+0)**.
No feature lost a single call.

