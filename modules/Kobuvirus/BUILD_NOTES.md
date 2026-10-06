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

**1207 sequences added**, every one at >=80% identity to a reference that
annotates the cut itself.

| key | projected | new | was in collection |
|---|---|---|---|
| 2B | 163 | 107 | 85 |
| 2C | 209 | 131 | 103 |
| 3A | 182 | 112 | 93 |
| 3C | 212 | 143 | 100 |
| 3D | 266 | 172 | 222 |
| Ldr | 234 | 191 | 60 |
| VP0 | 258 | 161 | 182 |
| VP1 | 269 | 172 | 222 |
| VPg | 51 | 18 | 43 |

### Features the cleavage-site test refused

Not projected, because the cut cannot be placed with confidence:

- `2A P1 Q 48% (n=118) P1' V53 C41 G12 NOT CONSERVED (48%)`
- `VP3 P1 P 53% (n=116) P1' Q61 H55 NOT CONSERVED (53%)`

### What it cost

Coverage over the same 313 exemplars: **3734 calls -> 3738 (+4)**.
No feature lost a single call.

### Why `2A` and `VP3` were refused: two incompatible annotation conventions

This is not divergence. The source `mat_peptide` records place the *same* cut
at two different residues, and the site test is what exposed it. Counting the
six residues before each annotated start and the five after:

```
VP3    LRYVAP|QHWKT   61      VP3 begins AT the P1 glutamine
       RHVTKQ|HWKTR   46      VP3 begins AFTER it -- the real Q|H 3C site

2A     FIVKVQ|RPTYV    3      2A begins at the 3C site
       AQRPTY|VHWAL   55      2A begins five residues downstream
```

Both forms describe one genomic position; they disagree by one residue for VP3
and by five for 2A. Pooled, P1 reads "P 53%" and "Q 48%" and the feature looks
unconserved, which is why neither is projected -- a projected cut would inherit
whichever convention the nearest reference happened to use, and the collection
would carry both.

The minority form is the correct one in each case: 3C cleaves *after* P1, so
VP3 starts at H and 2A starts at R. Normalising the majority to the motif would
make both projectable, but it changes boundaries on 61 and 55 source-annotated
genomes, so it is a curation decision and is left for a curator to make
deliberately rather than taken here. Both features keep their full collections
(VP3 116, 2A 96) and are unaffected at runtime.
