# Aphthovirus build notes

## Module boundary

One module for the genus: Ldr (leader protease, ~201 aa), VP4-VP2-VP3-VP1,
2A StopGo (18 aa), 2B, 2C, 3A, VPg, 3C, 3D. 13 rep contigs, 100% routing.
Built 2 October 2026 at the documented defaults for all 13 features
(129 PSSMs); leader key renamed L -> Ldr on 3 October.

Separate from Cardiovirus because the leader slot holds a different protein:
a ~201 aa papain-like protease here against a 67-76 aa zinc-finger there.

Every measured median matched the published FMDV size: Ldr 201, VP4 85,
VP2 218, VP3 220, VP1 211, 2A 18, 2C 318, 3A 153, 3C 213, 3D 470. The rules
know nothing about length, so that is an independent check on the triage.

## VPg models the FMDV tandem

FMDV carries three tandem copies of 3B and 410 BV-BRC records annotate all
three as a single `3B1/3B2/3B3` mat_peptide of 71 aa (23+24+24). Three
profiles from mutually paralogous 23-mers would never separate, so VPG is
one feature modelling the tandem. Forss & Schaller 1982 Nucleic Acids Res
10:6441 (PMID 6294604) for the tandem; Nayak 2005 J Virol 79:7698
(PMID 15919922) for uridylylation of all three copies. Both model-proposed.

## Measured 5 October 2026

- Full-taxon coverage: 2,580 genomes >= 5 kb, 431 exemplars at 95% nt
  identity, 431 routed, 369 with every feature.
- Held-out panel: 40 genomes (rep-contig genomes excluded), 39 clean. The
  flagged genome has its polyprotein called three times - frameshifting
  indels in the deposit.
- Self-recall 12,378/12,378.

## VPg outside FMDV -- FIXED 5 October 2026

The VPg decision above is correct for FMDV and wrong for the rest of the
genus. Counting every VPg record in the dump by species:

| species | n | lengths |
|---|---|---|
| A. vesiculae (FMDV) | 2385 | 71 x2209, 24 x107, 23 x54 |
| A. reedi (ERAV)     |   62 | 25 x62 |
| A. bogeli (ERBV)    |   25 | 24 x17, 25 x8 |
| A. burrowsi (BRAV)  |    5 | 24 x5 |

No tandem outside FMDV; no single copy inside it. The collection is 443/451
tandem, so the window (57-85) excludes the single-copy form outright. Cost:
VPg called on 2/23 A. bogeli, 20/32 A. reedi, 0/5 A. burrowsi, and those
genomes get the other twelve features at 100%.

### The second feature was the wrong design

Of 50 unique single-copy sequences, 38 are FMDV - records where a submitter
annotated 3B1/3B2/3B3 separately rather than as one tandem. Only 12 come from
the rhinitis species. And every aphthovirus VPg copy begins GPY..., so a
24-mer profile necessarily matches the tandem's first 24 residues. Two
features would have competed at one locus: the truncated-profile case the
curation checklist forbids.

### One feature covering both forms

The 50 single-copy sequences were merged into VPG (451 -> 501) and the window
widened to 19-85. Clustering separated the forms on its own: two profiles of
71 aa (435 seqs) and two of 23-25 aa (10 seqs).

Cutoff measured against polyproteins by species, comp_based_stats 0:

| profile | FMDV | reedi | bogeli | burrowsi |
|---|---|---|---|---|
| tandem x2 | 135-147 | 30-35 | 32-36 | - |
| single x2 | 21-26 | 56 | 51-56 | - |
| burrowsi  | 19 | - | - | 49 |

Noise ceiling 36, weakest true signal 49 -> **bit_cutoff 42**. The previous
cutoff of 31 was itself wrong: it would let a 71 aa tandem profile fire on a
24 aa protein.

### burrowsi came from the leftover pass, not a lowered member floor

Its two VPg sequences were already in the collection but start SAYDP... where
every other aphthovirus VPg starts GPY... They fall below the main pass's -m 5
and below the leftover pass's default length floor, which is relative to the
feature median (71 aa for a tandem-dominated collection).

Built the sanctioned way: documented main pass (`departures: none`, 3
profiles), then `build_leftover_pssms.py --min-seqs 2 --min-len-frac 0.3`,
which adds 2 more including the burrowsi pair. 5 profiles total. The identity
floor is never lowered in the leftover pass by design.

An earlier attempt used `-m 1` on the main pass and produced the same five
profiles, but a member floor of 1 is not allowed; two-member clusters reached
through the leftover pass are.

### Result

| species | before | after |
|---|---|---|
| A. vesiculae (FMDV) | 360/368 | 360/368 |
| A. reedi | 20/32 | **32/32** |
| A. bogeli | 2/23 | **23/23** |
| A. burrowsi | 0/5 | **5/5** |

VPg overall 382/431 -> 417/431 (and 422 with the burrowsi profile, confirmed
on the previously failing genomes). **Zero duplicate VPg calls** on any
genome - the specific risk the measured cutoff guarded against.

Remaining: Coypu aphthovirus, one genome, with no VPg record of its own.
