# Why the segment windows are tight, and what they are guarding

On a circular module these bounds are not only the quality script's length
check.  `annotate_by_viral_pssm.pl` consults them before it doubles a contig,
and refuses to double anything outside them, because **a partial contig is not
a circle**: doubling a fragment creates a junction the virus never had and
invites a feature across it.  They are the only guard there is, so they are set
from the measured distribution of COMPLETE genomes rather than loosened for
safety -- loosening them is what would be unsafe.

## Orthohepadnavirus: 3130 - 3400

17,629 complete genomes.  Percentiles 5 / 50 / 95 / 99.9 = 3131 / 3215 / 3221 /
3341.  Known complete genome sizes: HBV 3182-3248, woolly monkey 3179, bat
hepadnaviruses 3149-3215, woodchuck 3320, ground squirrel 3311.

**Below the bound.** 992 contigs fall in 2900-3130.  All but two are named
simply "Hepatitis B virus", the commonest lengths being 3055 (129 records) and
3049 (46), which is 130-190 nt short of a genome -- near-complete deposits, not
a shorter species.  They are refused, which costs the wrap recovery on 5.2% of
Orthohepadnavirus contigs: P and L cross the origin, so on those records they
will not be called at all.  That is the intended trade.  A record missing 130 nt
may be missing them from inside P, and a join built across an invented junction
is a worse answer than a missing call.

**Above the bound.** 47 contigs exceed 3400, the smallest at 3419 and running to
10,123.  Hepatitis B is frequently deposited as a head-to-tail 1.3-mer or 2-mer
construct; a 2-mer is ~6,364 nt and a 1.3-mer ~4,137, and both are refused here.
Doubling either would report every gene twice.

## Avihepadnavirus: 2980 - 3120

318 complete genomes, 3006-3048 with a single 3690 outlier.  DHBV is 3027 nt
(NC_001344.1), heron 3018, stork 3024.  The window is deliberately narrower
than the mammalian one because the genus is: a 2-mer would be ~6,050 nt and is
refused.

## What neither window can catch

A contig truncated but still inside its bounds -- a 3,200 nt HBV record missing
20 nt, say -- is doubled and gains a false junction.  Length cannot detect that,
and nothing else in the pipeline tries to.  It is recorded here as a known limit
rather than papered over.
