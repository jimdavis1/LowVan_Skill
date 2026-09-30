# Quinvirinae — build notes

## P14 is nested inside CP, and that is genome organisation, not cross-talk

`qc_duplicate_calls.py` flags `P14` as NESTED in `CP`. The question that settles
what to do about it is whether the two features share homology — if they do, one
profile set is calling the other's gene and a feature has to go; if they do not,
the module is looking at two genuinely overlapping ORFs and the right answer is
to declare the pair expected.

They do not share homology. Measured 30 September 2026, both directions:

| test | result |
|---|---|
| blastp, 170 P14 members × 915 CP members, E ≤ 1e-3 | **0 HSPs** of 155,550 pairs |
| same, relaxed to E ≤ 10 with `-comp_based_stats 0` | 276 HSPs, best 24.3 bits, E = 0.67, ≤ 43% query coverage |
| each of the 4 P14 PSSMs vs all CP members | best **21.7 bits** (P14.3 → CP.1_2); P14 `bit_cutoff` is 30 |
| each of the 32 CP PSSMs vs all P14 members | best **22.4 bits** (CP.21 → a P14.1 member); CP `bit_cutoff` is 80 |

No profile of either feature comes within 8 bits of the *other* feature's own
cutoff, so neither can call the other's gene. The relaxed hits are short
low-identity fragments of the kind any two 100–400 aa proteins produce once
composition correction is switched off; that switch is on in the annotator.

**So no CP alignment and no P14 alignment overlap in similarity.** The pair that
overlaps does so in coordinates. On GRSPaV Shihezi-1 (ON868740, 8825 nt):

```
CP  (N)   7847-8623   frame +2   777 nt / 259 aa
P14       8277-8651   frame +3   375 nt / 125 aa
overlap 8277-8623 = 347 nt  ->  92.5% of P14 lies inside CP
                                44.7% of CP lies inside P14
```

Those are the annotator's own coordinates, from a live run of
`annotate_by_viral_pssm.pl` on this contig, not a tblastn of a collection
member — the QC reads what the annotator emitted, so that is the number that
has to clear the floor.

Different frames, 92.7% containment. The 80% containment floor exists to stop
the QC reporting the routine short overlaps of compact viral genomes as defects
(adjacent genes in Closterovirus HSP70/HSP90 share 92 of 1503 nt); 92.7% clears
it easily, which is why this one was reported. It is still not a defect.

The overlap is confined to one host lineage. Of 170 P14 members, 169 are
*Grapevine rupestris stem pitting-associated virus* (taxon 196400) and one is
GRSPaV-2 (2893791) — all four P14 alignments are GRSPaV. All 171 GRSPaV coat
proteins sit in a single CP alignment, **CP.2**. So the co-occupying pairs are
`CP.2 × P14.1`, `CP.2 × P14.2`, `CP.2 × P14.3`, `CP.2 × P14.4`, and the other
31 CP alignments never meet a P14 at all.

### The exact cases, named

There is nothing to find in the clusters or the alignments — that is the point.
The collision is between two *calls*, each made correctly by its own profile.
From live runs of `annotate_by_viral_pssm.pl` on the two GRSPaV rep contigs:

| genome | CP profile | bits | CP CDS | P14 profile | bits | P14 CDS | overlap | % of P14 |
|---|---|---|---|---|---|---|---|---|
| ON868740 (Shihezi-1) | `Quinvirinae.CP.2.pssm` | 495.0 | 7847–8623 (+2) | `Quinvirinae.P14.4.pssm` | 210.3 | 8277–8651 (+3) | 347 nt | 92.5% |
| KT000366 (TI_21) | `Quinvirinae.CP.2.pssm` | 476.1 | 7771–8547 (+2) | `Quinvirinae.P14.2.pssm` | 223.5 | 8228–8584 (+3) | 320 nt | 89.6% |

`CP.2` is the CP profile in both because it is the only one of the 32 that holds
GRSPaV members (all 171 of them). The P14 profile differs between the two
genomes — `P14.4` on one, `P14.2` on the other — which is itself evidence that
nothing systematic is wrong with any single P14 alignment: whichever P14 profile
is nearest fires, and each clears the P14 cutoff of 30 by seven-fold while the
CP call clears 80 by six-fold. Both are right.

The two proteins LowVan emitted at that locus on ON868740 — a 259 aa
`Nucleocapsid protein` and a 125 aa `14 kDa protein` — share no sequence at
all: blastp one against the other reports **zero alignments at E ≤ 1000**.

`P14.4`, the profile that fired on ON868740, is the two-member alignment. Its
members are `fig|196400.1392.peg.6` and `fig|196400.1603.peg.6`, both 119 aa,
both annotated `14 kDa protein` in BV-BRC, 90% identical to each other. They
are genuine P14s; there is no mis-assigned sequence in the cluster.

### The same architecture appears elsewhere in this module

On KT000366 the annotator also calls `ORF2a` at 6816–6947 inside `TGB1` at
6578–7240 — **100%** containment, again out of frame. Whatever is decided for
`P14:CP` applies to `ORF2A:TGB1` on the same reasoning.

**Action:** pass `--expected P14:CP` to `qc_duplicate_calls.py` for this module.
Neither feature is dropped and neither profile set is rebuilt.

Do not generalise this to other nested pairs. The measurement above is what
distinguishes a real overlapping ORF from a profile calling someone else's gene,
and it has to be made each time.

## Annotation vocabulary

`CP` shipped as `Coat protein` with `gene_symbol` `CP`. Corrected 30 September
2026 to the controlled string `Nucleocapsid protein` with symbol `N`, which is
what `annotation-vocabulary.tsv` carries for every other taxon in the kit and
what Allexivirus, Alphaflexiviridae, Tobamovirus and Trivirinae already used.

`ORF2A`, `ORF5A` and `P14` still emit `ORF2a protein`, `ORF5a protein` and
`14 kDa protein`. `check_annotations.py` flags all three: the first two are bare
positional symbols and the third is a mass name. Open, pending a decision on what
to call an accessory whose function is unknown.
