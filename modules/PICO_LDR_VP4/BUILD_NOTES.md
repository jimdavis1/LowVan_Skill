# PICO_LDR_VP4 build notes

Built 5 October 2026. Six genera sharing one genome organisation:
**Bopivirus, Erbovirus, Mischivirus, Hunnivirus, Tottorivirus, Anativirus.**
Leader present; capsid cleaved to VP4+VP2+VP3+VP1; one VPg.

396 genomes (254 >= 5 kb). No genus here holds enough sequence to model
alone - Tottorivirus has 10 genomes, Anativirus 38 - so they are grouped and
their diversity becomes extra profiles inside each feature.

15 features, 95 PSSMs, 25 rep contigs routing 218/250 = **87.2%**
(28 clusters holding 32 genomes uncovered).

## Two slots held two different proteins

A slot name is a homology group INSIDE a genus and not across one. Checking
per-genus medians inside the group before building caught two cases where a
single feature would have blurred unrelated proteins:

| slot | split off | lengths |
|---|---|---|
| 2A | **Bopivirus** -> `2A_LONG` | 212 aa vs 12-55 for the StopGo peptides of the rest |
| Ldr | **Erbovirus** -> `Ldr_PRO` | 219 aa vs 51-84 for the rest |

`2A` keeps the vocabulary string "Protein 2A, StopGo peptide" and `2A_LONG`
takes "Protein 2A". `Ldr` is "Leader protein"; `Ldr_PRO` is "Leader protease
(Lpro)" - ERBV's L protein has C-terminal processing activity (Hinton 2002
J Gen Virol 83:3111, PMID 12466488, model-proposed).

Two features in one module must never share an annotation string: the quality
tool keys %essential by the string, so the second read silently overwrites the
first's length window. make_group_json.py refuses to emit a module with a
duplicate.

## Known gaps

- 87.2% routing is the lowest of the four grouped modules and the 32 unrouted
  genomes are a long tail of singletons across six genera.
- Every citation is model-proposed and sits in PMID_claude_generated.
