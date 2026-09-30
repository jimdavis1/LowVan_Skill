# Carlavirus — build notes

The largest module in the kit after Hepaciviridae: 316 profiles over six
features, from 3,184 BV-BRC records. Potato virus S and M, carnation latent
virus, chrysanthemum virus B and about ninety more. Replicase, a triple gene
block, coat protein, nucleic-acid-binding protein.

## Split out of Quinvirinae on measurement, not organisation

All six Quinvirinae genera share the replicase / TGB1-3 / CP / NABP layout, so
the organisation rule alone would have kept them together. Routing separated
them: lumped, 25 references reach 81.0% overall, and Carlavirus was consuming
a budget the five smaller genera needed. Split, every genus improved —
Carlavirus by eleven points, Robigovirus from 81.8% to 100%.

## Routing: accepted at 76.3%, 30 September 2026

**Curator's decision: ship at 76.3% and document it. Not tabled, not split
further, not given extra references.**

714 genomes at ≥6 kb fall into 102 clusters at 70% nucleotide identity, and
the distribution collapses immediately:

    155, 88, 44, 44, 22, 21, 21, 20, 18, 10, 9, 8, ... 22 singletons

The top two clusters hold a third of the genus. After the eighth, nothing has
more than 18 members, and 77 of the 102 clusters average 2.2 genomes each.

    references   genomes reached   of 714
         5            353          49.4%
        10            443          62.0%
        25            543          76.1%   <- project budget
        50            632          88.5%
        75            682          95.5%
       102            714           100%

Marginal return per reference: 70.6 for the first five, 5.6 at the budget
boundary, 2.0 from 51-75, 1.2 beyond. The budget of 25 lands in the flat part
of the curve.

**What this costs.** 171 genomes — 23.7% of the genus — match no reference and
therefore receive no LowVan annotation at all. They are not mis-annotated;
they are unannotated. Everything the module *does* route, it annotates well:
86.6% of scored genomes are genuinely clean, all six features are essential on
measurement at 93.3-99.3%, and the duplicate-call audit is empty.

**Why not fixed.** Two options were considered and declined for now:

- *Raise the budget for this module.* 50 references would reach 88.5%. It
  works, but it breaks the project-wide rule that makes "25 references" mean
  the same thing across every coverage audit.
- *Split Carlavirus in two.* Two modules of 25 references cover the same
  diversity as one of 50, which is the mechanism that fixed Quinvirinae,
  Tymoviridae and Betaflexiviridae. The obstacle is that there is no
  *organisational* boundary to cut on — every carlavirus has the same gene
  layout — so the split would be purely on clustering and needs a defensible
  line rather than an arbitrary one. That measurement has not been made.

The figure is stated in §04 of the coverage audit with a warning verdict, so
the limit travels with the module rather than living only here.
