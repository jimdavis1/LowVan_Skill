# Ampelovirus build notes

## 2026-10-02 — short accessory proteins discarded by the old length floor

`accessory_clusters()` used to drop any homology group whose median length fell
under 45 aa, on the prior that a profile that short cannot be given a specific
cutoff. That floor filtered clusters **before printing them**, so the casualties
never appeared in a log. The floor is now 0 and the groups are reported and
measured (`scripts/test_profile_specificity.py`).

Re-running the same clustering over the preserved unassigned pool (164
sequences) with the floor removed recovers three groups that were discarded:

| group | members | median | range | dominant source labels |
|---|---|---|---|---|
| KDAA | 52 | 36 aa | 28–39 | `4 kDa protein` ×31, `p4 protein` ×10, `hypothetical protein` ×9 |
| KDAB | 6 | 39 aa | 36–42 | `4 kDa protein` ×5, `hypothetical protein` ×1 |
| KDAC | 5 | 36 aa | 36–41 | `4 kDa protein` ×3, `putative small protein` ×1, `p4` ×1 |

**They are one gene, not three.** The three groups share no genome — 52, 6 and 5
sequences on 52, 6 and 5 distinct genomes, pairwise overlap zero — so they are
divergent alleles of a single 4 kDa gene, not three paralogues. Merged into one
feature, `P4`, 63 sequences on 63 genomes. (The group keys above are the
clustering's own and collide with the unrelated `KDAB`/`KDAC` accessory features
already in this module; they are not those features.)

This was the most expensive casualty of the floor in the class: 52 sequences on
52 genomes is larger than most accessories that did ship, and
`UNASSIGNED_TRACKING.tsv` had already flagged it loudly — `4 kDa protein` ×133
occurrences, `5 kDa protein` ×173.

### Caveat on completeness

The groups above are what the *member* floor (5) also admits. Anything that
failed both floors is still invisible, and the earlier logs cannot be mined for
it because the length filter ran before any printing. The pool is preserved at
`work/accessory_candidates/` (Box) and `/private/tmp/lowvan_work/rescue/` — the
complete list comes from re-running `accessories.py` against it, not from logs.

Related short groups that *cleared* the old floor and were already built:
`KDAB` (30 members, median 45 aa — cleared by one residue), `KDAC` (29, 51),
`P5` (22, 46), `HYPOTHETICALB` (21, 60), `P6A` (13, 51), `P6B` (11, 56),
`HYDROPHOBIC` (11, 57), `P6C` (7, 54), `PUTATIVE` (7, 48), `KDAH` (5, 71).

### Outcome — measured, and shipped

The merged `P4` was measured with `test_profile_specificity.py` rather than
accepted or dropped on the length prior. Scored against all 3,459 unique
proteins in the module, the single highest-scoring non-member (33.4 bits) was
itself annotated `4 kDa protein` — a real P4 that the 5-member floor had left
out. It was folded into the collection, which is why the feature holds 64
sequences and not 63.

First pass, 2 cluster profiles, worst on-target 24.9 against a best off-target
of 24.9 — not separable. Building the leftover profiles fixed it (21 leftovers
→ 5 more profiles, 7 total):

```
per-subject view   best on target  min 37.3  10th 57.2  median 80.7
                   best off target median 16.4  99th 21.8  MAX 24.9
                   -> CALLABLE
```

Shipped at `bit_cutoff: 30`, `min_len` 25, `max_len` 45.

**A 36 aa protein on 64 genomes is callable.** The length floor was wrong about
this one, and it was the largest single thing the floor was costing this module.
