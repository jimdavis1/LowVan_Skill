# Do the 36 genera with no module of their own route anyway? (6 October 2026)

404 genomes >=6 kb from the 36 Picornaviridae genera that no module claims,
fetched and run against the full 14-module runtime. `--min-complete 6000`: a
picornavirus genome is ~7.5 kb and the default floor for genomes of unknown
module is 10 kb, which once produced a wrong "94% unrouted" on the no-genus set.

**271 exemplars, 144 routed (53%), 127 unrouted.
Genome-level: 197 of 404 routed — 0.5% of the family, already being annotated
by modules built for other genera.**

## 18 genera route, and several are annotated nearly completely

| genus | exemplars | routed | mean features | lands in |
|---|---|---|---|---|
| Rosavirus | 48 | 48 | 7.1 | PICO_NOLDR_VP4 |
| Boosepivirus | 15 | 15 | **13.9** | PICO_LDR_VP4 |
| Rabovirus | 12 | 12 | 4.1 | Kobuvirus 9, Sapelovirus 3 |
| Parabovirus | 12 | 12 | 9.8 | Sapelovirus |
| Gallivirus | 11 | 11 | 7.9 | PICO_LDR_VP0 9, Kobuvirus 2 |
| Ludopivirus | 7 | 7 | 10.1 | Kobuvirus |
| Mupivirus | 8 | 8 | 7.0 | Cardiovirus |
| Aalivirus / Grusopivirus | 5 / 5 | 5 / 5 | 9.8 | Avihepatovirus |
| Sakobuvirus | 5 | 5 | 9.8 | Kobuvirus 4, PICO_LDR_VP0 1 |
| Ailurivirus | 5 | 5 | 2.8 | PICO_LDR_VP0 |
| Gruhelivirus / Crahelivirus | 3 / 1 | 3 / 1 | 8.3 / 9.0 | Hepatovirus |
| Malagasivirus, Mosavirus, Felipivirus, Livupivirus, Pygoscepivirus | 1-2 each | all | 6-11 | various |

**Rosavirus and Boosepivirus have ZERO mat_peptides of their own** and are
nonetheless annotated at 7.1 and 13.9 features per genome, entirely from
profiles trained on other genera. This is the strongest evidence in the family
that module partitioning by genome organisation — not phylogeny — is right.

## What does not route is mostly one genus

| genus | unrouted exemplars |
|---|---|
| **Megrivirus** | **64** |
| Mosavirus | 10 |
| Kunsagivirus | 8 |
| Crohivirus | 7 |
| Orivirus, Fipivirus | 6 each |
| Torchivirus | 5 |
| 15 further genera | 1-3 each |

Megrivirus is half the unrouted set on its own, and it does not route because
its genome is genuinely a different size: **median 9,142 nt against ~7,500 for
a typical picornavirus**, 138 genomes, 18 species.

## Megrivirus is below the budget bar

Measured on its own redundancy clustering, 109 genomes >=6 kb in 64 clusters,
**48 of them singletons**:

| greedy budget | genomes reached |
|---|---|
| 10 refs | 49/109 (45%) |
| **25 refs (project budget)** | **70/109 (64%)** |
| 40 refs | 85/109 (78%) |

Closteroviridae and Tymoviridae were tabled at 74.6%. At the project budget
Megrivirus reaches 64%, so by the standard already applied it is **tabled, not
built** — unless the budget is raised for it deliberately. It also has only 12
mat_peptide records across 194 proteins, so most features would come from
projection off one or two reference genomes.

## Conclusion

The remaining gap is smaller than the genus count suggests. Of 36 uncovered
genera, 18 already route and get 3-14 features with no further work. The real
remainder is Megrivirus plus a long singleton tail across ~20 genera of 1-10
genomes each — the documented tail the 25-reference budget cannot reach, not a
defect in the modules.
