# Curating the alignments

The pipeline produces drafts. Curation is where the module is actually made, and
it is manual work in an alignment editor. The goal, in the LowVan author's own
words: **clean, compact alignments with strong N-terminal conservation and good
C-terminal conservation.**

## The curator's checklist

From the LowVan README, non-negotiable before a PSSM is considered complete:

1. **Every profile starts at the right Met** — or you understand why it does
   not, and it is flagged.
2. **Coordinates are right** — is it starting from the correct Met, not an
   internal one?
3. **No PSSM is a truncated version of a longer PSSM.** Full-length only.
4. **The sequences actually encode the protein the key claims.** Blast against
   NR when in doubt. Source annotations are wrong often enough to matter.
5. **BV-BRC is missing proteins these viruses encode.** They exist in the
   literature and NR. Pull exemplars in rather than pretending the gene is
   absent.

Point 4 is the one that bites hardest, because a mislabelled cluster looks
perfectly healthy. Automate it with `scripts/qc_cross_feature.py`.

## Reading the Curation_Report

Per alignment: columns cut, sequences with end dashes removed, N-terminal stop
column and residue, and the N-terminal flag. Two things to know:

- The `Cut_Cols` and `Seq_w_End_Dashes_Removed` labels were **swapped** in the
  upstream header. Check the patch is applied before trusting the columns.
- `Nterm_Flag` set means the alignment needs a human. It does not mean the
  alignment is broken.

## The Met rule: flag, never delete

A profile should start at the initiating Met. The wrong way to achieve that is
to trim forward until you find one — that deletes real protein when the
truncated form happens to be the majority.

Concretely: an alignment whose correct N-terminus is `MFK/RV` but where the
majority of members are truncated will, under a search-for-Met rule, get trimmed
to a downstream Met and lose the correct start entirely. The rule that shipped
here trims only on **occupancy and conservation**, then reports where it stopped
and whether that column is a Met:

```perl
my $occ  = $col_sum{$k} || 0;
my $cons = $occ ? ($aa_sum->{$k}->{$mc_res} / $occ) : 0;
my $occf = $occ / $n_seqs_internal_dash;
if ($occf < $n_occ || $cons < $frac_n_term) {
    $exclude{$k} = 1; $n_term_trimmed++;
} else {
    $nterm_stop = { col => $k, res => $mc_res, occ => $occf,
                    id => $cons, kept => $occ };
    last;
}
```

Non-Met starts are rare but real. A flagged alignment you check by hand is
correct; a silently truncated protein is not.

## Splitting on N-terminal multimodality

The most common real problem: one cluster contains several distinct N-termini —
`MIL…`, `MVPW…`, `MVL…`, `MVP…` — so no single profile has a conserved start.
Percent identity does not separate them, because they are similar over the body
and differ only at the end that matters.

What works: cluster on the **first ~25 residues only**, which is what the
pipeline's `-nterm-mmseq-id` reclustering does. Curating by hand, the same idea:

- group by the leading motif, largest groups first
- aim for **just enough** splits to get good N-terminal conservation, not a
  profile per variant
- give the remainder a final catch-all bin rather than forcing it into a group
- keep the ends as similar as possible within each bin — that is the objective,
  not maximising within-bin identity

A split that produces two strong profiles and one weak catch-all beats six
narrow ones.

## Identical sequences

Collapse identical rows before building the profile — duplicates bias the
position-specific counts toward whichever strain was sequenced most. Keep a
**floor of 2**: if collapsing leaves one row, put a second back, because a
one-sequence PSSM is not a profile.

## Jalview and the `/start-end` suffix

Jalview writes ids as `fig|38767.13.CDS.1/1-533`. Strip the `/<start>-<end>`
suffix. It does not affect the profile — the PSSM's query id comes from the
pipeline, not the alignment header — but it breaks downstream id parsing and
makes the alignment set inconsistent.

Jalview also leaves `.bak001`, `.bak002` files when it saves over an existing
file. Those are the most reliable evidence of which alignments were hand-edited,
because cloud-sync flattens mtimes.

## Mature peptides and signal peptides

Derived features are coordinate operations on a curated parent alignment, not
binning products. They have no collection of their own, which is why
`qc_cross_feature.py` skips them and the registry reports them as `derived`.

**Mature peptide** (`G_MAT`): the parent minus the signal peptide. This works
well — 505 columns and 1000+ bits for Alpharhabdovirinae G. Set
`upstream_ext: 0` in the JSON: the mature N-terminus is a cleavage site, not a
start codon, and leaving the default sends the annotator upstream to the
precursor's own Met so both features land on the same coordinates. See
`references/json-schema.md`.

**Signal peptide** (`G_SP`): buildable, but **one profile per cluster, never one
for the feature**. A signal peptide is ~19–30 residues of hydrophobic core with
almost no conserved sequence, so a profile transfers badly — the
Alpharhabdovirinae lyssavirus profile scores 47 bits on rabies, 17 on Irkut
virus and 16 on Australian bat lyssavirus, and those are still lyssaviruses. A
single family-wide signal-peptide profile is hopeless; ten per-cluster ones
work, at 39–46 bits median with 100% self-recall.

Do not predict the cleavage site if you already have the mature peptide: the
signal peptide is simply the precursor minus the mature form, per sequence.
`scripts/build_sp_pssms.py` does that subtraction and builds one PSSM per
cluster. Only reach for SignalP where no mature alignment exists.

One measurement trap here cost this project a wrong conclusion. Scoring a
19-residue profile with `psiblast` against proteins applies composition-based
statistics, which roughly halves the score:

```
same PSSM, same subjects:  -comp_based_stats 2  ->  median 22 bits
                           -comp_based_stats 0  ->  median 43 bits
```

That 22 was read as "signal peptides cannot carry 30 bits" and the feature was
abandoned. The annotator's actual search is `tblastn` against the genome, where
the same profile scores 47. Measure with the search you will actually use.

## When to stop

An alignment is done when its ends are clean, it starts at a real Met (or is
flagged), it is full length, and its members are all the protein you think they
are. Then rebuild the PSSM — `scripts/rebuild_pssms.py` finds every alignment
whose profile is now stale.

## An output directory must describe exactly one run

Every rebuild replaces a feature directory, and the tempting way to write that
step is to enumerate what to copy: the clusters, the alignments, the PSSMs, the
reports you know about. That is wrong, and it fails silently.

The pipeline does not emit a fixed set of files. `nterm_reclustering.log` and
`reclustered_alis/` appear only when some cluster has ragged N-termini;
`truncated_alis/` and `truncation_qc/` only when something was truncated;
`Leftover_Seqs.fa` only when sequences failed to cluster. When a rebuild stops
needing one of them — usually *because the collection got better* — the old copy
is not overwritten. It stays, and it is now describing a build that no longer
exists.

`Betarhabdovirinae/P` accumulated three runs in one folder:

```
09:47  nterm_evaluation.log                      run A  (a failed pass-1 attempt)
13:08  clusters/ alis/ pssms/ Curation_Report    run B  (what is actually installed)
21:21  nterm_reclustering.log run.log tmp/       run C  (a week earlier)
       truncation_qc/ *.mmseq_*
```

Read straight, that folder says cluster 3 was split into two N-terminal
subclusters. It was — in run C. Run B's `nterm_evaluation.log` says cluster 3 is
clean and `reclustered_alis/` is empty, because the re-exported collection fixed
the ragged termini that caused the split. Nothing flagged the contradiction; the
mtimes were the only evidence, and mtimes are exactly what a synced drive
flattens.

So: clear the destination first, then copy the whole run directory across. Keep
only files a human wrote — `BUILD_PARAMS`, `PROVENANCE`, `FLAGGED_NTERM`,
`rejected/` — and let the pipeline own everything else.

```python
KEEP = {"BUILD_PARAMS", "FLAGGED_NTERM", "PROVENANCE", "rejected"}

def mirror(src, dest):
    for name in set(os.listdir(dest)) - KEEP:
        p = os.path.join(dest, name)
        shutil.rmtree(p, ignore_errors=True) if os.path.isdir(p) else os.remove(p)
    for name in os.listdir(src):
        s, d = os.path.join(src, name), os.path.join(dest, name)
        shutil.copytree(s, d) if os.path.isdir(s) else shutil.copy2(s, d)
```

This is the same principle as *curation must be a rule, not a file*, pointed the
other way. There, a decision that lives only in a file is lost on the next
rebuild. Here, an **artefact** that survives a rebuild it did not come from is
worse than lost — it is evidence for a claim nobody made.
