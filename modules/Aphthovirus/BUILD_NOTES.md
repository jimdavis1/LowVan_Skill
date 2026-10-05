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

## Known gap: VPg outside FMDV  -- OUTSTANDING

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

Fix is a second feature for the single-copy form, NOT a module split (the
rest of the layout is shared). Two cautions before building:
 1. The FMDV tandem CONTAINS three single-copy-length sequences, so a 24-mer
    profile will score at the FMDV 3B locus. That is the truncated-profile
    case the curation checklist forbids; measure cross-feature first.
 2. The two features need DISTINCT annotation strings. Two features sharing
    one string collide in viral_genome_quality.pl's %essential table - the
    defect that cost Tobamovirus a whole quality run.
