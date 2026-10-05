# PICO_NOLDR_VP4 build notes

Built 5 October 2026. **Cosavirus, Dicipivirus, Tremovirus.** No leader;
capsid cleaved to VP4+VP2+VP3+VP1; one VPg. Cosavirus is a human enteric
virus, Tremovirus is avian encephalomyelitis virus.

181 genomes (134 >= 5 kb). 13 features, 68 PSSMs, 25 rep contigs routing
133/134 = **99.3%**.

## One slot held two proteins

| slot | split off | lengths |
|---|---|---|
| 2A | **Tremovirus, Dicipivirus** -> `2A_LONG` | 141-164 aa vs 31 for Cosavirus |

Cosavirus 2A keeps "Protein 2A, StopGo peptide"; the long form takes
"Protein 2A".

## Known limitation

Dicipivirus has a **dicistronic** genome - two ORFs with an intergenic IRES
rather than one polyprotein. This module represents its polyprotein products
only and does not model that arrangement. The genus was grouped here on its
capsid and replicase layout; whether it deserves its own module on genome
organisation grounds is an open question and 30 genomes is thin either way.

Every citation is model-proposed and sits in PMID_claude_generated.
