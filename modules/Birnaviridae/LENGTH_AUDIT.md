# Length audit -- Birnaviridae

Collections: `Birnaviridae/collections/Birnaviridae`

| feature | n | min | median | max | frag | giant | split | lit | reference |
|---|---|---|---|---|---|---|---|---|---|
| `PEP11` | 11 | 11 | 13 | 13 |  |  |  |  | PMID 11836417 |
| `PEP46` | 42 | 41 | 46 | 46 |  |  |  |  | PMID 11836417 |
| `PEP7A` | 14 | 7 | 7 | 9 |  |  |  |  | PMID 11836417 |
| `PEP7B` | 7 | 7 | 7 | 7 |  |  |  |  | PMID 11836417 |
| `POLY` | 403 | 940 | 1012 | 1139 |  |  |  |  | PMID 11836417 |
| `PVP2` | 307 | 486 | 512 | 545 |  |  |  |  | PMID 11836417 |
| `VP1` | 379 | 785 | 879 | 1045 |  |  |  |  | PMID ['12477876', '30484762', '12069523'] |
| `VP2` | 398 | 417 | 441 | 467 |  |  |  |  | PMID 11836417 |
| `VP3` | 227 | 215 | 257 | 308 |  |  |  |  | PMID 11836417 |
| `VP4` | 188 | 212 | 243 | 314 |  |  | 1 |  | PMID 11836417 |
| `VP5` | 192 | 105 | 145 | 154 |  |  |  |  | PMID 12477876 |
| `X` | 2 | 71 | 88 | 106 |  |  |  |  | PMID 12477876 |

## Findings

**Unresolved: 0.**

## Signed off (1)

Checked and accepted, with the reason. These still describe the
collections -- they are not suppressed, only answered.

| kind | feature | member | measurement | why it is correct |
|---|---|---|---|---|
| SPLIT | `VP4` | fig | 15x 212-243 aa + 2x 314-314 aa | The short mode (212-243 aa) is Avibirnavirus/Aquabirnavirus/TV-1 and is correct. The long mode is now a single genome, Largemouth bass birnavirus at 314 aa (2 identical deposits). Corrected 21 Nov 2026: this row previously also covered Blotched snakehead virus at 305 aa and asserted that both modes were 'validated at both ends'. That was wrong. qc_cross_feature found the BSNV projected member 100% identical over its first 71 residues to fig|311176.4.mat_peptide.6, the deposited X peptide -- the projection had started at the pVP2|X junction instead of the X|VP4 junction, because in BSNV the polyprotein is pVP2-X-VP4-VP3 (PMID 12477876) and the projection skipped X. The 71 aa prefix was trimmed; the member is now 234 aa, which is both the published BSNV VP4 size and the length implied by X's own C-terminal A|A -- two independent confirmations. It had its own single-sequence cluster and therefore its own profile, so uncorrected it would have called VP4 on X's locus in every BSNV-lineage genome. Largemouth bass birnavirus is left at 314 aa: it deposits no X peptide, its one internal A|[AS] would imply a 49 aa X rather than a 71 aa one, and no publication places a cut there. Projecting one would be invention, so the length stands and is recorded as a known uncertainty. |
