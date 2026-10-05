# PICO_LDR_VP0 build notes

Built 5 October 2026. **Sicinivirus, Salivirus, Oscivirus, Passerivirus,
Shanbavirus.** Leader present; VP0 left uncleaved (so VP0-VP3-VP1, three
capsid proteins not four); one VPg. Salivirus causes human gastroenteritis;
the rest are bird viruses.

281 genomes (219 >= 5 kb). 13 features, 67 PSSMs, 25 rep contigs routing
196/218 = **89.9%**.

## One slot held two proteins

| slot | split off | lengths |
|---|---|---|
| Ldr | **Sicinivirus** -> `Ldr_LONG` | 462 aa vs 64-174 for the other four |

`Ldr` is "Leader protein". `Ldr_LONG` is declared **"Uncharacterized
lineage-specific protein"** with gene_symbol `Ldr-long`: a PubMed search found
no functional description for the Sicinivirus leader, and the skill forbids
inventing a designation. The protein is modelled and called; it is simply not
named beyond its lineage. Revisit if a description appears.

## Known gaps

- 22 genomes unrouted (22 clusters), a long tail of singletons.
- Every citation is model-proposed and sits in PMID_claude_generated.
