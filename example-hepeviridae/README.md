# Hepeviridae — a worked example

The first module of the **Alsuviricetes** class build, and the only family in
that class that infects humans: *Paslahepevirus balayani* is hepatitis E virus,
around 20 million infections a year.

Built 17 September 2026. Four features, 63 PSSMs, 24 rep contigs.

```
self-recall (in-range)   97.5 / 99.4 / 99.4 / 100%
held-out panel           42 genomes, 115 calls, 0 duplicates
whole taxon              669 exemplars, 622 routed
protein distribution     ORF1 100.8%  ORF2 99.4%  ORF3 94.9%  ORF4 10.6%
genome quality           555 good (89.2%), 67 poor, ~2% attributable to the module
string collapse          157 BV-BRC strings -> 4 annotations, 39x
vocabulary saturation    free text knee at 533 genomes, controlled at 13
```

## The programs

| file | what it decides |
|---|---|
| `build_collections.py` | which annotation strings mean which protein, plus the ORF4 homology gate |
| `build_hepeviridae_json.py` | what the taxon has, how long, how confident a match must be, and every citation's provenance |
| `features.json` | the per-feature annotation strings handed to the pipeline |

Everything else is the kit. Eight programs this build needed that are **not**
taxon-specific were promoted into `skill/scripts/` instead of living here:
`fetch_features.py`, `fetch_seqs.py`, `fetch_contigs.py`, `subset_dump.py`,
`build_features.py`, `rescue_unassigned.py`, `reroute_outliers.py` and
`qc_truncation_symmetric.py`.

## Three decisions worth reading

**ORF1 is declared whole.** RefSeq annotates seven domains on it. That
annotation descends from Koonin 1992 (PMID 1518855), whose title says
*computer-assisted* and which calls the domains *putative*; Perttilä 2013
(PMID 23255617) and LeDesma 2023 (PMID 36852909) both looked for processing
experimentally and found none. Declaring seven mature peptides would assert
seven cleavage events the field has not agreed happen. The boundaries are
recorded in the build notes so the decision can be revisited.

**ORF4 is two different proteins.** Paslahepevirus ORF4 and Rocahepevirus ORF4
share a positional label and nothing else — all-vs-all BLASTp over every
ORF4-annotated sequence in the dump gives a best cross-genus alignment of **5
residues at e-value 138**. Only the Rocahepevirus protein is declared; the
Paslahepevirus one has 2 unique sequences and is retired.

The rule is therefore **not scoped by genus**, because genus is unreliable here:
an audit of every ORF1 against every other found 16 of 1,417 genomes whose
BV-BRC genus label contradicts its own sequence — rat and ferret HEV deposited
twice under two taxon ids, appearing as both Paslahepevirus and Rocahepevirus at
96–99% identity. A `HOMOLOGY_GATE` against a verified seed does the
discrimination instead, and 13 of the 45 sequences it admitted are labelled
Paslahepevirus and are really rat HEV.

**ORF3's `bit_cutoff` is 90, measured on both sides.** The noise ceiling is 78.6
— a second HSP of the chosen profile landing in ORF2's frame, which the
extension flags then grew to ORF2's exact span and which was the only duplicate
in the whole evaluation. The weakest genuine call is 99.0. The cost is recorded
rather than hidden: ORF3 is no longer called on three Chirohepevirus genomes
that scored 62.6, 77.2 and 84.9. That is a coverage hole, not a threshold
problem — the ORF3 collection holds 863 Paslahepevirus sequences against 2
Chirohepevirus and 0 Piscihepevirus.

## Citations

`PMID` is **absent on all four features**; every id is in
`PMID_claude_generated`, and the `annotation-vocabulary.tsv` rows carry an empty
PubMed column. Every id was checked to resolve to the paper named. None has been
read and confirmed by a curator, which is the only act that clears the flag.
