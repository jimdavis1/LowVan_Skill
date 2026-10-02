# Enterovirus build notes

## Module boundary

One module for the genus. Two independent lines of evidence:

- **Genome organisation is identical across all seven species groups.** No L
  protein; L-less P1-P2-P3 with 11 mature peptides. The per-group median
  lengths agree to within a few residues — 3C is 182–183 aa in every one of
  the seven groups, 2C is 322–330, 3D is 457–462.
- **Routing fits the budget without splitting.** Lumped Enterovirus reaches
  **97.3% of genomes at 25 references** (5:73%, 10:87%, 15:93%, 20:95%).
  Splitting by species group reaches 99–100% each but costs seven reference
  budgets to buy two points.

So the EV_* species subsets exist for measurement only; they are not modules.

## Triage

Picornavirus naming is positional (1A–1D, 2A–2C, 3A–3D) *and* functional
(VP4/VP2/VP3/VP1, VPg, picornain, polymerase), often both in one string.
Unlike the U-number trap these positional names **are** homology groups — but
only inside a genus. Enterovirus 2A is a chymotrypsin-like protease,
cardiovirus 2A is not one, and aphthovirus 2A is an 18 aa StopGo peptide. None
of these rules transfer to another genus.

Rules assign **99.0% of 55,534 mat_peptides**. A string matching two different
features is refused, not arbitrated — which caught `coat protein VP4 (P1B)`
and `1D protein, VP4`, both self-contradictory on their face.

### The ambiguous strings were placed by homology, not by guess

425 proteins are annotated `protease` with no number. That is 2A or 3C and the
string cannot say which. Rather than split them on a prior, they went to a pool
and were placed by blastp against the collections the unambiguous strings had
already built, at the project floor of ≥80% identity over ≥85% query coverage:

```
protease                      408  ->  3C   408 (100%)
proteinase                     11  ->  2A     8 (73%),  3C 3 (27%)
capsid protein                 34  ->  VP1   27 (79%),  VP4 5, VP3 2
VP1 / 3D / RdRp (as CDS)      251  ->  VP1 / 3D        (100%)
```

`protease` alone means 3C in this genus, unanimously over 408 sequences. That
is a measurement, and it is not transferable either.

The same floor correctly **refused** the precursors — `VP4-2C` ×36, `P1` ×30,
`P1-P2` ×25 — because a precursor cannot cover 85% of a single mature peptide.
They are held out in `PRECURSORS.fasta`, not deleted: they are real proteins,
just not the product any feature models. 71 of 70,226 proteins end up unplaced.

## Length windows

Set from the **per-species-group medians**, not from the observed extremes, and
checked against the published sizes rather than fitted to them:

| | VP4 | VP2 | VP3 | VP1 | 2A | 2B | 2C | 3A | VPg | 3C | 3D |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median | 69 | 261 | 238 | 296 | 150 | 99 | 329 | 86 | 22 | 183 | 461 |
| published | 69 | 272 | 238 | 302 | 149 | 98 | 329 | 87 | 22 | 183 | 461 |

Every median is within 4% of the literature value and most are exact, from
rules written on annotation strings alone — an independent check that the
triage is right, since nothing in the rules knows about length.

6,294 sequences fall outside the windows. **98.4% of them are short** — a 17 aa
"VP4", a 7 aa "VPg polypeptide", a 470 aa "polyprotein" — i.e. fragments from
partial genome records. Only 100 are long. The windows are excluding
truncations, which is what they are for.

## Literature

Every PMID was checked against PubMed before being written in. Three of four
recalled from memory were wrong, pointing at papers on *E. coli* sex factor F,
HIV epidemiology and AML stromal layers. Recalled PMIDs are not evidence.

| PMID | | category |
|---|---|---|
| 6264310 | Kitamura 1981 *Nature* 291:547–553 — poliovirus primary structure and gene organisation | coordinates |
| 6202874 | Toyoda 1984 *J Mol Biol* 174:561–585 — all three serotypes, cleavage map | coordinates |
| 3011278 | Toyoda 1986 *Cell* 45:761–770 — 2A is a second, distinct proteinase | function |
| 2252396 | Palmenberg 1990 *Annu Rev Microbiol* 44:603–623 — proteolytic processing | function |
| 217003  | Rothberg 1978 *PNAS* 75:4868–4872 — VPg–pUp tyrosine linkage | function |

## VPg

22 aa, and callable — measured, not assumed. Against 300 full-length
polyproteins over 5,052 profile × subject pairs: at the VPg locus 10th
percentile 21.2 bits, median 27.0; anywhere else 99th percentile 17.6, max
19.3. `bit_cutoff: 20` calls it and is silent across the other 2,160 residues.

## Held-out evaluation

Panel drawn with `make_panel.py` from the build's own BV-BRC contigs,
stratified across species, **with all 25 reference genomes excluded by
accession**. 120 genomes requested; 89 scored — `make_gto.py` lost 31 to
BV-BRC ID-server 403s, which is rate limiting, not a module property.

```
89 genomes | 83 clean (93.3%) | 85 of 89 called all 12 features
```

| | POLY | VP4 | VP2 | VP3 | VP1 | 2A | 2B | 2C | 3A | VPg | 3C | 3D |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| called | 89 | 87 | 89 | 88 | 89 | 88 | 88 | 89 | 88 | **87** | 88 | 88 |
| % | 100 | 97.8 | 100 | 98.9 | 100 | 98.9 | 98.9 | 100 | 98.9 | **97.8** | 98.9 | 98.9 |

Every genome that fell short is flagged for ambiguous bases — the four not
calling exactly 12 features are 147711.7824 (6), 147711.9528 (10),
1193974.973 (11), and 12089.2262 (14, a doubled polyprotein). Only two flags
appear across the whole panel: `Contig has too many ambiguous bases` ×6 and
`too many HSPs for: Genome polyprotein` ×2.

**VPg at 97.8% is the result worth keeping.** A 22 aa feature was expected to
be uncallable on a length prior; measuring it instead gave a cutoff of 20 bits
that separates cleanly, and it now calls on 87 of 89 unseen genomes. The two
misses are ambiguous-base contigs, not threshold failures.

The two doubled polyprotein calls are the one open defect. `qc_cross_feature.py`
is the next step there.
