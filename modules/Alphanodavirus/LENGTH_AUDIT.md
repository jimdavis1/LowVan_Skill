# Length audit -- Alphanodavirus

Collections: `./collections/Alphanodavirus`

| feature | n | min | median | max | frag | giant | split | lit | reference |
|---|---|---|---|---|---|---|---|---|---|
| `A` | 11 | 973 | 998 | 1043 |  |  |  |  | PMID 10799587 |
| `ALPHA` | 8 | 399 | 407 | 407 |  |  |  |  | PMID 6548308 |
| `B1` | 4 | 102 | 116 | 131 |  |  |  |  | PMID 15163762 |
| `B2` | 12 | 90 | 106 | 143 |  |  | 1 |  | PMID 15163762 |
| `BETA` | 8 | 355 | 363 | 363 |  |  |  |  | PMID 9765417 |
| `GAMMA` | 8 | 40 | 44 | 44 |  |  |  |  | PMID 9765417 |

## Findings

**Unresolved: 0.**

## Signed off (1)

Checked and accepted, with the reason. These still describe the
collections -- they are not suppressed, only answered.

| kind | feature | member | measurement | why it is correct |
|---|---|---|---|---|
| SPLIT | `B2` | fig | 7x 90-106 aa + 5x 137-143 aa | Two length classes, 90-106 aa and 137-143 aa, are the documented biology of this genus and not a mixed bin. PMID 15163762 states in one sentence "FHV, 106 amino acids; NoV, 134 and 137 amino acids", and PMID 12573589 detected both Nodamura forms (B2-134 and B2-137) by MALDI-TOF, arising from the first and second AUG of subgenomic RNA3. Pariacoto B2 is 90 aa (PMID 10799587). So the short mode is FHV/BBV/Boolarra/Pariacoto and the long mode is Nodamura and its close relatives -- species-level length variation in one protein, with every member in the +1 frame relative to protein A. Splitting the bin would build two profiles for one gene; the alternative-initiation pair within Nodamura cannot be separated by any profile anyway, since B2-134 is B2-137 minus three N-terminal residues. |
