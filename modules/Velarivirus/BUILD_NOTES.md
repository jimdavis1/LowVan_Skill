# Velarivirus build notes

## 2026-10-02 — short accessory proteins discarded by the old length floor

`accessory_clusters()` used to drop any homology group whose median length fell
under 45 aa, on the prior that a profile that short cannot be given a cutoff
that is specific. That floor filtered clusters **before printing them**, so the
casualties never appeared in a log. The floor is now 0 and the groups are
reported and measured (`scripts/test_profile_specificity.py`).

Re-running the same clustering over the preserved unassigned pool (142
sequences) with the floor removed recovers three groups that were discarded:

| group | members | median | range | dominant source labels |
|---|---|---|---|---|
| HYPOTHETICAL | 29 | 31 aa | 31–36 | `hypothetical protein` ×20, `p4` ×5, `ORF2` ×2 |
| P4A | 8 | 36 aa | 36 | `p4 protein` ×7, `P4` ×1 |
| P4B | 5 | 35 aa | 35 | `p4 protein` ×3, `hypothetical protein` ×2 |

**They are one gene, not three.** The three groups share no genome — 29, 8 and 5
sequences on 29, 8 and 5 distinct genomes, pairwise overlap zero — so they are
divergent alleles of a single p4 gene that mmseqs could not join at 25% identity,
not three paralogues. They were merged into one feature, `P4`, 42 sequences on
42 genomes.

Source labels were no help here and were not trusted: the largest group is
mostly `hypothetical protein`, and `p4` appears in all three. The co-occurrence
test decided it.

### Caveat on completeness

The groups above are what the *member* floor (5) also admits. Anything that
failed both floors is still invisible, and the earlier logs cannot be mined for
it because the length filter ran before any printing. The pool is preserved at
`work/accessory_candidates/` (Box) and `/private/tmp/lowvan_work/rescue/` — the
complete list comes from re-running `accessories.py` against it, not from logs.

Related short groups that *cleared* the old floor and were already built:
`PUTATIVE` (11 members, median 51 aa), `P7` (5, 62), `KDAC` (7, 81),
`HYPOTHETICALF` (15, 91), `HYPOTHETICALH` (12, 92), `P12` (11, 89).

### Outcome — measured, and shipped

The merged `P4` was measured with `test_profile_specificity.py` rather than
accepted or dropped on the length prior. First pass, 4 cluster profiles:

```
per-subject view   best on target  min 22.5   best off target MAX 24.1
                   -> a single cutoff loses 1 of 42
```

Every off-target hit came from one profile grazing a 517 aa HSP90 at 21–24
bits — a protein its own HSP90 profiles call at several hundred bits, so it was
never at real risk of being mis-assigned, but it sat above the weakest P4.

Building the leftover profiles fixed it outright (`build_leftover_pssms.py`, 6
leftovers → 1 more profile, 5 total):

```
per-subject view   best on target  min 41.3  10th 64.9  median 75.3
                   best off target median 17.1  99th 22.5  MAX 24.1
                   -> CALLABLE
```

Shipped at `bit_cutoff: 30`, `min_len` 28, `max_len` 39. Measured against all
1,310 unique proteins in the module, not only the P4-bearing genomes.

**A 31 aa protein is callable.** The length floor was wrong about this one.
