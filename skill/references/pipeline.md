# The clustering / alignment / PSSM pipeline

`Other_Scripts/fasta-cluster-pssm-2.pl` in the Viral_Annotation repo. It reads a
protein FASTA on stdin and produces clusters, alignments and PSSMs.

```bash
cat collections/<Module>/<FEAT>.fasta \
  | fasta-cluster-pssm-2.pl -a "<Annotation string>" -p <Module>.<FEAT> \
      -m 5 -f 0.33 -n 0.75 -c 0 -fd 0.20 -e 3 -efo 0.15 -mi 0.8 -mc 0.8 \
      -p_nterm 65 -n_nterm 3 -nterm_eval_length 10 -nterm-mmseq-id 0.7
```

## What it does

1. drop ambiguous-residue sequences, deduplicate, sort by length
2. MMseqs2 cluster at `-mi` identity / `-mc` coverage
3. keep clusters with at least `-m` members; the rest go to `Leftover_Seqs.aa`
4. MAFFT-align each cluster
5. trim: remove high-gap columns, trim low-occupancy ends, log every action
6. re-evaluate N-terminal conservation; recluster the first
   `-nterm_eval_length` residues of any alignment that fails, splitting it
7. `psiblast -in_msa` each final alignment into a PSSM

## Output layout

```
<FEAT>/
  clusters/<n>.fasta            one FASTA per surviving cluster
  alis/<n>.fa                   raw MAFFT alignment
  corrected_alis/<n>.fa         after trimming -- THIS is what you curate
  pssms/<Module>.<FEAT>.<n>.pssm
  alis-for-reclustering/        alignments that failed the N-term check
  reclustered_alis/<n>_<k>.fa   products of the N-term split
  truncated_alis/               dropped by the truncation QC
  Curation_Report               per-alignment record of every trim
  Leftover_Seqs.{aa,fa}         sub-min_seqs sequences
  nterm_evaluation.log, nterm_reclustering.log, run.log
```

`corrected_alis/` is the contract: PSSMs are built from it, and it is the file
you hand-edit. **Never re-copy the pipeline's output over `corrected_alis/`
after curation has started** — you will silently destroy hand work.

## Parameters that actually matter

`-m 5` (minimum cluster size) is a target, not a rule. Go down to `-m 2` rather
than lose a taxon. If you cannot call all of a taxon's proteins at any setting,
you do not have a module for that taxon.

`-mi 0.8` (clustering identity) is usually the real constraint when coverage is
poor, not `-m`. On Betarhabdovirinae, whose core proteins are very divergent
between species, dropping `-mi` transformed coverage:

| Feature | at `-mi 0.8` | after |
|---|---|---|
| P | 20% | 71% at `mi 0.4` |
| M | 30% | 82% at `mi 0.4` |
| G | 64% | 100% at `mi 0.5` |
| MOV | 63% | 74% at `mi 0.4` |

Record every departure from defaults in a `BUILD_PARAMS` file next to the
feature, with the numbers that justified it. The registry artifact marks those
features so a reader knows they were tuned.

`-n_occ` (default 0.40) is the N-terminal occupancy floor. Raising it is the
blunt instrument for jagged N-termini; splitting the cluster is the sharp one.

## The clustering identity floor

`-mi` is MMseqs2's `--min-seq-id`, passed straight through:

```
mmseqs easy-cluster $tmp.fasta $tmp.mmseq tmp --min-seq-id $mmi -c $mmc --cov-mode 0
```

It is the identity two proteins must share **to enter the same cluster**, and one
cluster becomes one alignment becomes one PSSM. So `-mi` is not a tuning dial for
coverage — it is the threshold at which you declare two sequences to be the same
protein. The documented default is 0.8.

**Never go below 0.6.** When a feature will not cluster at the default, the
temptation is to walk `-mi` down until something comes out, and it always will:
a low threshold reaches high coverage by admitting anything. Betarhabdovirinae P
at `-mi 0.3` produced **one PSSM for 315 sequences at 85% "coverage"** — a single
profile smeared across proteins that share less than a third of their residues.
That number looks like success and means nothing; the profile self-recalls its
own collection and generalises to nothing.

Two other things break underneath you before you get there:

- the truncation QC's `qc_pident 75` is calibrated on the assumption that
  within-cluster identity is at least 0.8, so below ~0.75 it is comparing
  consensuses that were never that similar to begin with;
- column conservation falls far enough that the curation step strips the
  alignment, and you end up with a short profile built from a wide clade.

When a collection is thin, **vary `-m` instead** — a smaller minimum cluster size
reaches sparse data without redefining what counts as the same protein. A ladder
of `-m 5/3/2 x -mi 0.8/0.6` covers the useful space. Assert the floor in code so
it cannot drift:

```python
MI_FLOOR = 0.6
assert all(float(mi) >= MI_FLOOR for _m, mi in LADDER), "ladder breaches the -mi floor"
```

A feature that produces no profile at `-mi 0.6` should produce no profile. That
is a real finding about the collection, not a failure to tune.

## Upstream bugs and patches

Patches against the repo's `Other_Scripts/` are in `assets/patches/`. Seven are
worth knowing about even if you do not apply them.

**1. `mafft --reorder` on a single sequence.** MAFFT exits 1, and the guard
`if (@seqs > 0)` lets a one-sequence alignment reach it, so the `die` aborts the
whole run before any alignment is written. The fix is one character: `> 0` →
`> 1` in `FCP_Main_Utils.pm`. Without it, a single-member cluster kills the job.

**2. Hash-order nondeterminism in the dedup step.** The identical-sequence
collapse iterated a Perl hash, so the surviving order varied between runs. That
destroys MAFFT's similarity ordering and makes PSSMs irreproducible — two runs
on identical input give different profiles. Replace with an order-preserving
pass:

```perl
my @ali6; my %seen;
for my $i (0..$#ali5) {
    my $seq = uc $ali5[$i][2];
    next if $seen{$seq}++;
    push @ali6, ([$ali5[$i][0], $ali5[$i][1], $seq]);
}
```

Keep a floor of 2 — a profile cannot be built from one sequence, so if dedup
collapses an alignment to a single row, put the second one back.

**3. Ordering ties.** `readdir` order and length-sort ties both leak into the
output. Add an id tiebreak to every length sort
(`sort { $b->[3] <=> $a->[3] or $a->[0] cmp $b->[0] }`) and sort cluster
numbers numerically, not lexically, or `10.fa` sorts before `2.fa` and cluster
numbering shifts between runs.

**4. End-trim denominator.** The N- and C-terminal trims computed conservation
against the wrong denominator, so they trimmed on degraded columns rather than
bimodal ones. See `references/curation.md` for the corrected rule and, more
importantly, the Met policy: the trim must **flag**, never delete.

**5. Trim caps — the guards invert exactly when they are needed.**
`FCP_Ali_Utils.pm` removes sequences whose ends are all gaps, and aborts the
trim if too few sequences occupy the terminal columns. That guard is backwards
under load: heavy truncation is *what drives occupancy down*, so the alignments
that most need trimming are the ones where the guard fires and nothing happens.
The N-terminal side had the mirror problem — a column kept at 44% occupancy,
then 56% of the sequences deleted for starting with a dash.

Replace "abort if occupancy is low" with "trim, but never more than a fixed
fraction of the alignment":

```perl
my $max_end_frac = defined($options->{max_end_frac}) ? $options->{max_end_frac} : 0.35;
my @ragged = grep { substr($ali6[$_][2], -$end_dash) !~ /\w/ } (0..$#ali6);
```

with the same shape for `max_start_frac`. Exposed as `-max_end_frac` and
`-max_start_frac` in `FCP_Main_Utils.pm`, both defaulting to 0.35. The cap is
the safety property: a bad call costs at most a third of the alignment, and the
report says how much it took. Pair this with the Met policy in
`references/curation.md` — **trim to the common Met start, do not delete.**

**6. PSSMs outliving the alignment they came from.** `FCP_Nterm_Utils.pm`
rebuilds profiles after N-terminal reclustering. If the rebuild produced
nothing, the previous run's `.pssm` stayed on disk, and nothing downstream could
tell it apart from a current one — a profile attributed to an alignment that no
longer exists. The patch stashes the old file as `.superseded` first and
restores it only if no replacement appears, so a stale profile is never silently
presented as fresh. Same principle as *an output directory must describe exactly
one run* in `references/curation.md`.

**7. The annotator dies when a genome matches nothing.** Three defects in
`annotate_by_viral_pssm.pl` compound into one failure, and the symptom names
the wrong culprit.

`write_fasta` requires its first argument to be an array ref and `confess`es
otherwise, so calling it with an empty feature list kills the run:

```
Bad sequence entry passed to print_alignment_as_fasta at gjoseqlib.pm line 956
```

A genome that routes to a module and then matches no profile is a legitimate
empty result, but it produced **no output at all** -- indistinguishable from a
crash, and in a batch run it silently becomes a "coverage gap". Guard both
writers with `if @gene_seqs` / `if @prot_seqs`.

Underneath that, profile selection accepts zero scores:

```perl
my $best_bit = 0;
unless ($hsp_best_bit < $best_bit) { $best_pssm = $name; ... }   # true at 0 == 0
```

A profile matching nothing becomes `$best_pssm`, and because every later zero
also passes, the **last** profile in the directory wins. Use
`if ($hsp_best_bit > $best_bit)` and skip the feature when nothing clears
threshold. Note the diagnostic trap: the error names whichever profile happened
to sort last, so it looks like that profile is broken. It is not.

Third, the push into `@gene_seqs` is unguarded, so a feature with an empty gene
or protein string makes the same malformed triple by another route. Drop it
with a warning naming the feature instead.

Together these are why one genome in 1356 returned nothing. Individually each
is quiet; combined they turn "this virus has no recognisable proteins" into
"this virus failed to annotate", which is a different claim.

After fixing 2 and 3, verify reproducibility directly — run the pipeline twice
on the same input and diff every output file. Zero differing files is the bar.

## The truncation QC

`assets/patches/FCP_Trunc_Utils.pm` adds a step before curation: take each
cluster's consensus, blast them all against each other, and drop any cluster
whose consensus is a substring of a longer survivor. Defaults
`-qc_pident 75 -qc_qcov 0.95 -qc_lenratio 0.92`.

This exists because a cluster is often just a truncated form of a longer
protein, and its profile competes with the full-length one. It caught a real
case here: 22 rabies nucleocapsids at 346 aa, 99.7% identical to the 450 aa
consensus. Dropped clusters land in `truncated_alis/` so you can inspect them.

Note what this does *not* catch: a cluster of the wrong protein entirely, at
full length. That needs `scripts/qc_cross_feature.py`.
