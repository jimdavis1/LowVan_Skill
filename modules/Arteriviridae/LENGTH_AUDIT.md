# Length audit -- Arteriviridae

Collections: `./collections/Arteriviridae`

| feature | n | min | median | max | frag | giant | split | lit | reference |
|---|---|---|---|---|---|---|---|---|---|
| `E` | 587 | 45 | 73 | 84 |  |  |  |  | PMID 10400725 |
| `GP2` | 1852 | 187 | 256 | 304 |  |  |  |  | PMID 1328669 |
| `GP3` | 1981 | 163 | 254 | 269 |  |  |  |  | PMID 12743278 |
| `GP4` | 1777 | 81 | 178 | 185 | 1 |  | 1 |  | PMID 12743278 |
| `GP5` | 2375 | 161 | 200 | 255 |  |  |  |  | PMID 1328669 |
| `M` | 1650 | 158 | 174 | 180 |  |  |  |  | PMID 1328669 |
| `N` | 1287 | 100 | 123 | 134 |  |  |  |  | PMID 14604534 |
| `NSP1` | 6 | 260 | 260 | 260 |  |  |  |  | PMID 7769711 |
| `NSP10` | 97 | 413 | 441 | 467 |  |  |  |  | PMID 9971783 |
| `NSP11` | 40 | 219 | 223 | 225 |  |  |  |  | PMID 9971783 |
| `NSP12` | 39 | 118 | 143 | 164 |  |  |  |  | PMID 9971783 |
| `NSP1A` | 137 | 180 | 180 | 180 |  |  |  |  | PMID 20006994 |
| `NSP1B` | 140 | 203 | 205 | 205 |  |  |  |  | PMID 20006994 |
| `NSP2` | 964 | 571 | 1059 | 1251 |  |  | 1 |  | PMID 20696193 |
| `NSP3` | 543 | 226 | 230 | 231 |  |  |  |  | PMID 20696193 |
| `NSP4` | 7 | 203 | 204 | 204 |  |  |  |  | PMID 9371590 |
| `NSP5` | 10 | 170 | 170 | 170 |  |  |  |  | PMID 9371590 |
| `NSP6` | 5 | 16 | 16 | 22 |  |  |  |  | PMID 9371590 |
| `NSP7A` | 14 | 123 | 123 | 149 |  |  |  |  | PMID 17098961 |
| `NSP7B` | 164 | 102 | 120 | 120 |  |  |  |  | PMID 17098961 |
| `NSP8` | 9 | 45 | 45 | 50 |  |  |  |  | PMID 30353315 |
| `NSP9` | 11 | 640 | 685 | 693 |  |  |  |  | PMID 9971783 |
| `ORF1A` | 1805 | 1706 | 2391 | 2539 |  |  |  |  | PMID 20696193 |
| `ORF5A` | 375 | 43 | 46 | 60 |  |  |  |  | PMID 21307223 |

## Findings

**Unresolved: 0.**

## Signed off (3)

Checked and accepted, with the reason. These still describe the
collections -- they are not suppressed, only answered.

| kind | feature | member | measurement | why it is correct |
|---|---|---|---|---|
| FRAGMENT | `GP4` | fig|1118369.9.peg.6 | 83 aa vs fig median 178 | Same finding as the GP4 SPLIT above, reported per member: this IS the 83 aa wobbly possum disease virus GP4, measured against a module median of 178 that is set by PRRSV. Not a fragment. |
| SPLIT | `GP4` | fig | 4x 81-83 aa + 1772x 148-185 aa | CONFIRMED lineage biology: wobbly possum disease virus GP4 is 83 aa against 148-185 aa elsewhere (PRRSV-1 183, EAV 152). It is genuinely short, not a fragment -- the 83 aa form is what NC_026811.2 annotates and it is the only GP4 that genome has. GP4 forms the minor-glycoprotein heterotrimer with GP2 and GP3 (PMID 12743278). |
| SPLIT | `NSP2` | fig | 22x 571-576 aa + 942x 751-1251 aa | CONFIRMED lineage biology and the reason this module keeps one nsp2 feature: the 571-576 aa class is Alphaarterivirus (EAV nsp2 = pp1ab 261-831 = 571 aa, PMID 7622476 for the PLP2 protease and UniProt P19811 for the Gly831 bond) and the 830-1196 aa class is Betaarterivirus (PRRSV-1 1078, PRRSV-2 1196, Fang & Snijder Table 1, PMID 20696193). One orthogroup, two lineage length classes. Every member of both classes begins at the measured nsp1beta/nsp2 cleavage residue; the 279 internal typing amplicons that do not were removed by the junction filter. |
