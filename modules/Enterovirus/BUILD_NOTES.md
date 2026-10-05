# Enterovirus build notes

(2 October 2026 build; measured 5 October. See git history for the original
triage notes, which this file supersedes.)

## Module boundary

One module for the genus. L-less, VP4-VP2-VP3-VP1, chymotrypsin-like 2A
protease, 2B, 2C, 3A, one VPg, 3C, 3D. 12 features, 701 PSSMs, 25 rep
contigs. Largest module in the project: 26,364 genomes, 70,226 unique
protein sequences, 53,036 mature peptides.

Two independent lines of evidence for lumping the genus:
 - Genome organisation is identical across all seven species groups. 3C is
   182-183 aa in every one, 2C 322-330, 3D 457-462.
 - Lumped routing reaches 97.3% of genomes at 25 references (5:73%, 10:87%,
   15:93%, 20:95%). Splitting by species group reaches 99-100% each but
   costs seven reference budgets to buy two points.

The EV_* species subsets exist for measurement only; they are not modules.

## Triage

Rules assign 99.0% of 55,534 mat_peptides. A string matching two different
features is refused, not arbitrated - which caught `coat protein VP4 (P1B)`
and `1D protein, VP4`, both self-contradictory.

425 proteins annotated only `protease` are 2A or 3C and the string cannot
say which. Placed by blastp against the collections the unambiguous strings
had already built, at >=80% identity over >=85% query coverage:
  protease   408 -> 3C 408 (100%)
  proteinase  11 -> 2A 8, 3C 3
The same floor correctly REFUSED the precursors (VP4-2C x36, P1 x30,
P1-P2 x25) - a precursor cannot cover 85% of a single mature peptide. They
are held in PRECURSORS.fasta, not deleted. 71 of 70,226 proteins unplaced.

## Length windows

Set from per-species-group medians, checked against published sizes:

| | VP4 | VP2 | VP3 | VP1 | 2A | 2B | 2C | 3A | VPg | 3C | 3D |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median | 69 | 261 | 238 | 296 | 150 | 99 | 329 | 86 | 22 | 183 | 461 |
| published | 69 | 272 | 238 | 302 | 149 | 98 | 329 | 87 | 22 | 183 | 461 |

Every median within 4% of the literature, from rules written on annotation
strings alone. 6,294 sequences fall outside the windows; 98.4% are short
(fragments from partial genome records). Only 100 are long.

## VPg is 22 aa and was measured callable, not assumed

Against 300 full-length polyproteins over 5,052 profile x subject pairs:
at the VPg locus 10th percentile 21.2 bits, median 27.0; anywhere else on
the same polyproteins 99th percentile 17.6, max 19.3. `bit_cutoff: 20`
calls it and is silent across the other 2,160 residues.

NOTE the JSON key is VPG and the gene_symbol is VPg. The feature table
reports the SYMBOL. Counting by key returns 0 - this bit the audit
generator on 5 October.

## Measured 5 October 2026

- Coverage: 17,755 genomes >= 5 kb, 2,856 exemplars at 95% nt identity,
  2,856 routed (100%), 2,690 with all twelve features. Lowest feature VPg
  2,784/2,856 = 97.5%. No feature trails the best by 20 points.
- Held-out panel: 40 genomes (rep-contig genomes excluded), 36 clean.
  Three of the four flags are a polyprotein called twice (frameshifting
  indels in the deposit); one is a contig >1% ambiguous.
- Self-recall 63,428/63,527.
- No rebuild needed. Collections were already large enough that the
  projection technique used on the thin genera had nothing to add.

## The eight genomes that route and call nothing

All eight are broken deposits, in two forms:

| genome | name | verdict |
|---|---|---|
| 12080.3096 | Human poliovirus 1 | **0.0% C** |
| 12090.79 | Enterovirus D70 | **0.0% C** |
| 12090.87 | Enterovirus D70 | **0.0% C** |
| 138951.213 | Enterovirus D DW-RAT-200 | 94% ambiguous |
| 138951.219 | Enterovirus D DW-RAT-206 | 88% ambiguous |
| 147711.7818 | Rhinovirus A DW-RAT-152 | 92% ambiguous |
| 147711.7823 | Rhinovirus A DW-RAT-136 | 92% ambiguous |
| 147711.8186 | Rhinovirus A DW-RAT-319 | 92% ambiguous |

The first three are full length with ZERO ambiguous bases and zero
cytosine: every C has become T. 5' end reads TTAAAATAGTTTTGGGG... where an
enterovirus reads TTAAAACAGCTCTGGGG... These are bisulfite-converted
sequences deposited as genomes. The conversion destroys every reading frame
(longest ORF in the poliovirus genome: 158 aa against ~2,200 for a real
polyprotein), so no profile can match and none should.

No measured coverage gap in this module.
