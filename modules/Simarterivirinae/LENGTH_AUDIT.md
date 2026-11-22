# Length audit -- Simarterivirinae

Collections: `./collections/Simarterivirinae`

| feature | n | min | median | max | frag | giant | split | lit | reference |
|---|---|---|---|---|---|---|---|---|---|
| `E` | 169 | 75 | 80 | 81 |  |  |  |  | PMID 25036340 |
| `EP` | 23 | 52 | 94 | 94 |  |  | 1 |  | PMID 33917085 |
| `GP2` | 171 | 158 | 205 | 235 |  |  |  |  | PMID 25036340 |
| `GP2P` | 171 | 214 | 232 | 281 |  |  |  |  | PMID 9420301 |
| `GP3` | 166 | 120 | 179 | 270 |  |  | 1 |  | PMID 25036340 |
| `GP3P` | 171 | 198 | 204 | 220 |  |  |  |  | PMID 9420301 |
| `GP4` | 168 | 152 | 173 | 182 |  |  |  |  | PMID 25036340 |
| `GP4P` | 171 | 122 | 194 | 205 |  |  |  |  | PMID 9420301 |
| `GP5` | 171 | 230 | 245 | 278 |  |  |  |  | PMID 25036340 |
| `M` | 171 | 159 | 161 | 173 |  |  |  |  | PMID 25036340 |
| `N` | 171 | 96 | 111 | 132 |  |  |  |  | PMID 25036340 |
| `ORF1A` | 110 | 2048 | 2052 | 2105 |  |  |  |  | PMID 21544192 |
| `ORF5A` | 165 | 45 | 64 | 68 |  |  |  |  | PMID 21307223 |

## Findings

**Unresolved: 0.**

## Signed off (2)

Checked and accepted, with the reason. These still describe the
collections -- they are not suppressed, only answered.

| kind | feature | member | measurement | why it is correct |
|---|---|---|---|---|
| SPLIT | `EP` | fig | 10x 52-59 aa + 13x 94-94 aa | CONFIRMED lineage biology: E-prime genuinely comes in two size classes -- 94 aa in SHFV and 110 in Free State vervet virus, 52 aa in Pebjah and Southwest baboon virus -- and it is absent or truncated by a premature stop codon in four further simian arteriviruses (PMID 25187550). This is the labile member of the duplicated cassette and the reason it must not carry copy_num; Cai 2021 (PMID 33917085) shows ORF2b' is not required for replication in vitro. |
| SPLIT | `GP3` | fig | 164x 120-204 aa + 2x 270-270 aa | CONFIRMED record defect, not a bin defect, and it is documented upstream: the two 270 aa members are Zambian malbrouck virus 1 (NC_048209), the one simian arterivirus record that annotates no GP4 at all and whose GP3/GP4 boundary is unresolved. Every other GP3 in the subfamily is 136-204 aa. Kept rather than dropped because the sequence is a real GP3-region product; the module calls it as GP3 and the register is the submitter's. |
