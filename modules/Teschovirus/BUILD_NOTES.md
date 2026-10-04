# Teschovirus build notes

## Module boundary

One module for the genus: leader (Ldr, 86 aa), VP4-VP2-VP3-VP1, 2A StopGo
(~20 aa, NPG|P), 2B, 2C, 3A, one VPg, 3C, 3D. 7 rep contigs, 100% routing at
the budget measurement. Built 2 October 2026 at the documented defaults for
every feature (49 PSSMs); leader key renamed L -> Ldr on 3 October.

## Measured 4 October 2026

- Full-taxon coverage: 250 genomes >= 5 kb, 160 exemplars at 95% nt identity,
  160 routed, 148 with every feature. Every mature peptide >= 158/160 except
  Ldr 151/160. No feature trails by 20 points; no rebuild.
- Held-out panel: 40 genomes (rep-contig genomes excluded), 39 clean. The
  flagged one, 118140.487 (6,389 nt), has its polyprotein split into three
  HSPs by frameshifting indels in the deposited sequence.
- Self-recall 1508/1509.

## Known gaps

- Bat teschovirus leader: 7 of the 9 Ldr misses are Rousettus, Eidolon and
  Pteropodidae bat teschoviruses. No leader annotation exists for that lineage.
- 23 unique source leaders are 118 aa: the 86-aa leader plus 32 residues from
  an upstream in-frame AUG (93% identical over the shared span). Modelled as
  the 86-aa form; the window (68-103) excludes the long form by design.
