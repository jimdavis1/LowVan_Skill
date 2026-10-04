# Avihepatovirus build notes

## Module boundary

One module for the genus. L-less, VP0 left uncleaved, three 2A proteins
(2A1 StopGo, 2A2 AIG1-like, 2A3 parechovirus-like), one VPg. 3C is the only
protease. 9 rep contigs route 100% of genomes at the budget measurement.

## Rebuild, 4 October 2026

The first build (2 October) had 2 to 5 sequences per mature-peptide
collection: 196 unique polyproteins, but only five genomes carry mat_peptide
records. Full-taxon coverage called VPG on 68% of exemplars and VP1 on 86%.

**Projection.** `project_matpeptides.py` takes each reference genome whose
mat_peptides are exact substrings of its own polyprotein, globally aligns
every other polyprotein to its nearest reference that places the feature,
and carries the cut coordinates across. Accepted only when the polyprotein is
>= 50% identical to the reference, the segment is inside the JSON window and
>= 35% identical to the reference segment. Nothing is predicted.

| feature | before | after |
|---|---|---|
| 2A1 | 2 | 18 |
| 2A2 | 3 (285 aa, lumped) | 94 (161 aa) |
| 2A3 | none | 62 |
| 2B | 3 | 51 |
| 2C | 5 | 92 |
| 3A | 4 | 61 |
| 3C | 2 | 55 |
| 3D | 5 | 117 |
| VP0 | 5 | 100 |
| VP1 | 5 | 120 |
| VP3 | 5 | 96 |
| VPG | 4 | 25 |

Cut sites are crisp: every projected 3C, 2C and VP0 ends in Q; lengths are
constant per feature except VP1 (237-240).

**2A2 / 2A3.** RefSeq and BV-BRC annotate a single 285-aa "2A2". The DHAV-3
genome 1006063.235 annotates 2A2 (161 aa) and 2A3 (124 aa) separately, 2A3
starting at the conserved HLPR motif, which is present at the same position in
every DHAV-1 and -3 sequence. Split with `--split '2A2:2A3:HLPR[LIV]'`.
Organisation: Tseng, Knowles & Tsai 2007, Virus Res 123:190 (PMID 17067712);
Tseng & Tsai 2007, Virus Res 126:19 (PMID 17292992). Both model-proposed,
in PMID_claude_generated as coordinates.

Windows 0.8-1.2x and bit cutoffs max(20, 0.45x) of the new collection medians.
All features built at the documented defaults except 2A1 (-m 2, 18 seqs).
25 PSSMs.

## Measured after the rebuild

- Coverage: 22 exemplars (95% nt clusters, >= 5 kb), 22 routed, 12 with every
  feature. 2A3 called on all 22.
- Held-out panel: 40 genomes, rep-contig genomes excluded, 39 clean. The one
  flagged genome, 885378.3, has a premature stop inside 2A3 in the deposited
  sequence; POLY and 2A3 are correctly withheld.

## Known gaps

- **Divergent lineage.** Duck egg-reducing syndrome virus, a mute-swan virus
  and four goose/crane *Avihepatovirus sp.*: 13 of 196 polyproteins, 7 of 22
  exemplars, 8.3-9.0 kb genomes, < 40% identical to any reference, no
  mat_peptide records. VPG and 2A2 are missed on all seven. A projection test
  at 25% identity drifted on 2B, 3A, VP3 and VPG ends, so it was not shipped.
  Needs annotated genomes from the lineage.
- **3C starts one codon late.** `check_cleaved_ends.py` over 196 genomes: every
  junction is exactly adjacent in 91-100% of genomes, except 3C start vs VPG
  end, which sits at +3 nt in 94%. The 3C call loses its first residue; with
  both extensions at 0 the annotator does not walk back to the cut. Systematic,
  one codon, not drift.
