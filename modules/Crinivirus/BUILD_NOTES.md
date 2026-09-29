# Crinivirus — binning decisions

Dump: 1,749 unique proteins, 3,903 features, 1,051 genome records. Bipartite.
Binned 3,116 of 3,903 = **79.8%**.

## A mass is not a homology group

The accessory proteins are named by mass, and the masses do not correspond to
homology groups. Measured by clustering the mass-named pool at 25% identity /
50% coverage:

    "p22"        FOUR non-homologous groups, at 188, 191, 192 and 193 aa
    "p26/27/28"  one group holding all three names, plus two further
                 separate groups also calling themselves p26/p27/p28
    "p59"/"p60"  ONE group of 78 sequences
    "p6/7/8/9"   at least seven groups, fragmenting by size

So accessory features are **discovered by clustering and only then named**,
with A/B/C suffixes where one name covers several groups — the Ampelovirus
P20A/P20B precedent. Floors: 8 members (enough to build a profile at `-m 2`)
and 45 aa median (enough to give it a usable cutoff). The discovered "p4"
group is 33 aa and was dropped: a 33-residue profile cannot both fire on its
own members and stay silent across a 9 kb genome.

## p59/p60 and "movement protein" are the HSP90 homolog

The 78-sequence p59/p60 group is the closterovirid HSP90 homolog, which
Closterovirus, Ampelovirus and Velarivirus all name HSP90 rather than by mass.
Shipping it as two mass features would have invented two proteins and split
one that exists.

The 13 records labelled `movement protein` are all exactly 518 aa and blastp
at **98.1–99.4% identity, E=0.0** against the HSP90 collection. In
closterovirids the HSP90 homolog is part of the movement machinery, which is
what that label is reporting. Binned to HSP90; a separate MP feature would
have duplicated it on the same locus.

## The ORF1a/1b fusion is not ORF1a

46 features are annotated as an `ORF1a/1b fusion polyprotein` at 2,413–2,499
aa — the +1 frameshift readthrough product. ORF1a alone is ~1,990 and the
RdRp ~475, and both are called as their own features, as in the shipped
Closterovirus module. POLY is therefore capped at 2,200 and the fusion
excluded explicitly; the rescue later recovered 90 of them into
`POLY.outliers.fasta`, median 2,483, where they are tracked and not used.

## Segment assignment — and the validation it buys

Anchoring RNA1 on the records carrying POLY/RDRP and RNA2 on those carrying
HSP70/CPM gives two sets with **zero overlap**, and every one of the 23
features falls entirely on one side:

    RNA1  POLY RDRP RNASE3 SUPPRESSOR P22A P22B P22C P23 P25 P6B P7A
    RNA2  CP CPM HSP70 HSP90 P26 P27 P28 P6A P7B P8A P8B P9

A discovered cluster that straddled both segments would have been an
artefact. None did. RNASE3 on RNA1 agrees with the SPCSV literature.

## Routing: the gate figure was inflated

`repcontig_budget.py --merge-by-name` reported **100.0%** at budget 25.
The shipped set covers **511/553 = 92.4%**. The gate merged records by genome
name, but 18 Crinivirus names carry 3–8 records — isolate collisions, not
segments. Same defect already corrected for Bromoviridae (94.8% → 69.5%).
`merge_segments.py` gets this right and refuses to merge them
(316 records logged `collision_not_merged`), so the two tools disagreed with
each other as well as with the truth. **92.4% is the honest figure** and
`gate_results.tsv` has been corrected.

## One record would have switched off a length check

`Plant associated crinivirus 1 MIR20SW`, 16,522 nt — 65% longer than the next
largest, almost certainly both segments on one contig. It set RNA2's
`max_len` to 17,348, which is no bound at all for the other 276 RNA2 records
(6,575–8,672 nt). Excluded from the window at 1.5x the 95th percentile. It is
not a rep contig, so routing is unaffected.

Segments as shipped: RNA1 6,548–9,877 nt; RNA2 6,246–9,105 nt.

## The evaluation is a sample, not the whole taxon

Annotation costs ~86 seconds per genome, so the evaluation that decides
`copy_num` and reports per-genome quality was run on a **random sample of up
to 150 genomes per module** (seed 29; Tymovirus has only 102 and so is
complete). Every genome already scored before the cap was kept in the sample
rather than discarded.

150 genomes is ample for the decision it feeds: the test for `copy_num` is
whether a feature is called on at least 85% of routed genomes, and at n=150 a
true rate of 85% is estimated to within about +/-6 points. Capillovirus was
already unambiguous at n=101, all three features calling at 100.0%.

**Whole-taxon coverage is measured separately** by `evaluate_coverage.py`,
which clusters the module's entire BV-BRC holding and annotates one exemplar
per cluster. That, not this sample, is the figure the coverage audit reports
for reach across the taxon.

## Correction, 29 September 2026: 23 features became 17

The feature set above is what was *built*. It is not what shipped.

A duplicate-call audit of the scored genomes found seven features firing
between 82% and 98% of the time inside a longer feature, and two pairs called
on identical coordinates. Five of the eleven features that had cleared the
85% essentiality threshold were counting another protein's locus as their own.

    RNASE3  inside SUPPRESSOR  97%   99% id / 100% cov  -> one gene, merged
    P26     inside P28         83%   23% id / 96% cov   -> one family, merged into P27
    P9, P8A, P8B, P6A, P23     82-98%, no similarity    -> dropped
    HSP70 + HSP90  identical coordinates, 7%            -> genuine paralogues, kept

`RNASE3` and `SUPPRESSOR` are one gene: in *Sweet potato chlorotic stunt
virus* the RNase3 **is** the silencing suppressor. `p26`, `p27` and `p28` are
one diverged family at 23-24% identity over 96-99% coverage — the mass trap
again, this time splitting one protein into three rather than merging three
into one. The five dropped features are 47-90 aa profiles with `bit_cutoff`
25-30 that fire on anything and have no sequence similarity to what they land
inside.

HSP70 and HSP90 stay. Their collections share no sequences, yet the proteins
are 64% identical over 99% of their length — genuine crinivirus paralogues,
and no cutoff separates them without losing one.

**Result: 15 features, quality 81.3% → 90.0% clean, seven features essential
instead of eleven.** This defect *raises* a call rate, which is why
self-recall (92-100% on every feature) and collection disjointness (zero
containment at 90% blastp) both missed it. It is now checked by
`qc_duplicate_calls.py`, run before `apply_copy_num.py`.

### Second correction: the first audit over-dropped

The duplicate-call check that produced the list above had **no containment
floor** — it counted any overlap with a longer feature as nesting, and
adjacent genes in a compact viral genome overlap routinely (Closterovirus
HSP70 and HSP90 share 92 bases). On that faulty measure `P8A` and `P8B` were
reported as 95% and 92% "inside HSP90" and were dropped.

With an 80% containment floor applied, neither is flagged at all. The real
relationships are `p9`↔`p8A` (33-34% identity over 91-99% coverage) and
`p6A`↔`p8B` (29-36% over 79-81%): two duplicate *pairs*, each of which should
survive as one merged feature. Both were restored.

**`P8A` is called on 100% of routed genomes and earns `copy_num`** — dropping
it was a real loss, caught only because the containment bug was found while
checking a different module.

Final: **17 features, 107 PSSMs, 90.0% clean, eight essential.**
