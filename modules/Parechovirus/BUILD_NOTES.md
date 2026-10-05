# Parechovirus build notes

## Module boundary

One module for the genus: VP0 (uncleaved), VP3, VP1, 2A, 2B, 2C, 3A, VPg,
3C, 3D, POLY. No leader; VP0 is never cleaved to VP4+VP2, so the capsid has
three proteins not four. 25 rep contigs. Built 2 October 2026 at the
documented defaults for all 11 features (66 PSSMs).

## Measured 5 October 2026

- Full-taxon coverage: 667 genomes >= 5 kb, 245 exemplars at 95% nt
  identity, 244 routed, 207 with every feature.
- Held-out panel: 40 genomes (rep-contig genomes excluded), 39/39 scored,
  all clean - no genome, contig or feature flag.
- Self-recall 2,931 sequences, 99%.

## Known gaps

- **3A on the bat parechoviruses.** 3A is called on 218/244. The misses are
  8/8 Scotophilus kuhlii, 2/2 Myotis bat parechovirus 1, 2/2 Myotis bat
  parechovirus 2 - every bat genome in the set - plus 8/202 human
  parechoviruses. The bat genomes route and get their other ten features, so
  this is a profile-coverage gap specific to 3A. No bat 3A sequence exists in
  the collection to build from.
- **One unrouted genome**: 3413225.3 Bos parechovirus 1 NZ1, 6,134 nt, the
  only bovine parechovirus in BV-BRC. A 26th reference would route it; the
  project budget is 25.
