# Changes to Viral_Annotation — Rhabdoviridae module work

Draft. Every entry names the file, what changed, and the measurement or
observation that prompted it. Patches against the unmodified originals are in
this directory as `<file>.patch`.

**Nothing here has been pushed.** The repository's clustering code is still
original: `Other_Scripts/FCP_Ali_Utils.pm`, `FCP_Main_Utils.pm` and
`FCP_Nterm_Utils.pm` on disk match `git HEAD`, and the patched copies live only
here. Three tracked files differ from HEAD, listed under "Already applied".

---

## Already applied to the working clone

### `annotate_by_viral_pssm.pl` — four defects, one visible symptom

All four surfaced from a single genome (`11292.18743`) out of 1,356 that
produced no output at all.

1. **Zero-scoring profiles were selected.** The test was
   `unless ($hsp_best_bit < $best_bit)` with `$best_bit` initialised to 0, which
   is true when both are zero. A profile matching nothing became `$best_pssm`,
   and since every later zero also passed, the *last* profile in the directory
   won. Now `if ($hsp_best_bit > $best_bit)`, with the feature skipped when
   nothing clears threshold.
   *Diagnostic note: the resulting error named whichever profile sorted last,
   which made a newly built profile look broken when it was not.*

2. **`write_fasta` was called with an empty list.** It requires an array ref and
   `confess`es otherwise, so a genome that routed but matched no profile — a
   legitimate empty result — killed the run and produced nothing, which is
   indistinguishable from a crash. Both writers are now guarded.

3. **Features with empty sequences were pushed unguarded**, producing the same
   malformed entry by another route. Now dropped with a warning naming the
   feature.

4. **The taxon assignment was lost on zero-feature genomes.** The GTO builder
   reads close-genome fields out of the *feature rows* of this script's stdout,
   so no rows meant no `close_genomes` and no `viral_family` in the GTO. The
   script now emits one `no_features_called` row carrying module, reference
   contig, accession and bit score. The GTO parser picks those up but builds no
   feature, because it only does so for `CDS` or `mat_peptide`.
   *Observed: a genome matching the rabies reference at 8,183 bits was being
   recorded as nothing at all.*

### `Viral_PSSM.json` — Rhabdoviridae module data replaced

Four modules rebuilt (Alpharhabdovirinae, Betarhabdovirinae, Dichorhavirus,
Novirhabdovirus). `Alphagymnorhavirus` removed. See "Module content" below.

### `Other_Scripts/New-annotate-viral-taxon.pl`

Your file, copied in unchanged from the working directory. Not modified here.

---

## Not applied — patched copies in this directory

### `FCP_Ali_Utils.pm`

- **End-trim caps.** The guard aborted the trim when too few sequences occupied
  the terminal columns, which inverts under load: heavy truncation is what
  drives occupancy down, so the alignments most needing a trim were the ones
  where nothing happened. Replaced with "trim, but never more than a fixed
  fraction" — `max_end_frac` and `max_start_frac`, default 0.35. The cap is the
  safety property: a bad call costs at most a third of the alignment, and the
  report says how much it took.
- **Retention floor separate from the formation floor.** `-m` decides which
  clusters are worth *forming*; a new `min_keep` (default 2) decides which
  formed-and-curated alignments are worth *keeping*. Gating both on `-m`
  discarded alignments that cleared the bar at clustering and then lost a
  sequence to curation — Alpharhabdovirinae M cluster 18 and Dichorhavirus M
  cluster 3 both went 6 → 4 against `-m 5` and were written with no profile.

### `FCP_Main_Utils.pm`

Exposes `-max_end_frac`, `-max_start_frac` and `-min_keep`, and documents them
in the usage block.

### `FCP_Nterm_Utils.pm`

- **Stale profiles could outlive their alignment.** If N-terminal reclustering
  produced nothing, the previous run's `.pssm` stayed on disk and nothing
  downstream could tell it from a current one. The old file is now stashed as
  `.superseded` and restored only if no replacement appears.
- **The parent alignment is restored with it.** It is moved to
  `alis-for-reclustering/` on the assumption children will replace it; when none
  do, restoring the profile without its alignment leaves a PSSM whose source is
  filed under a directory named for an input. Five Betarhabdovirinae alignments
  were in that state.
  *Measured alternative, rejected: lowering the N-terminal clustering identity
  rescued only 2 of those 5, and did so by orphaning the other members. A whole
  parent alignment beats a fragment of it.*

### `FCP_Trunc_Utils.pm` — new file, absent from the repository

Truncation QC run before curation: take each cluster's consensus, blast them
against each other, drop any cluster whose consensus is a substring of a longer
survivor. Defaults `-qc_pident 75 -qc_qcov 0.95 -qc_lenratio 0.92`. Caught 22
rabies nucleocapsids at 346 aa that were 99.7% identical to the 450 aa
consensus.

---

## Module content and conventions (data, not code)

- `-mi` never below **0.6**. It is the identity at which two proteins are
  declared the same protein, and one profile per cluster means a looser value
  buys coverage by admitting anything: Betarhabdovirinae P at `-mi 0.3` gave one
  profile for 315 sequences at 85% "coverage".
- Length bounds derived from **p1/p99** of the collection rather than min/max,
  because the collections still contain truncated proteins. Across the essential
  features this took members-outside-their-own-bounds from 384 to 15.
- `PMID` and `PMID_claude_generated` are **disjoint**: a model-proposed citation
  appears only in the latter, so a feature with nothing read carries no `PMID`.
- Uncharacterized proteins are **one feature per homology group** (`UNC1`…),
  each with its own profiles, all sharing the S1 Table annotation, so any one can
  be renamed when its function is established.
- `cleave` is set only on termini that are genuinely protease cut sites.
  `G_MAT` carried `downstream_ext: 0`, which drifted the mature C-terminus in
  **40% of 1,108 genomes**, median 28 codons and out to 1.2 kb. Corrected to
  `cleave="n"` — the signal-peptidase site is a real cut, the C-terminus is the
  precursor's own stop codon. `check_cleaved_ends.py` now measures this.
- Second-pass **leftover profiles** at a member floor of 2. Sequences too
  divergent to join a cluster at `-mi 0.6` were being discarded, and they were
  exactly the ones the module then failed to call: of 115 Alpharhabdovirinae
  genomes missing P, 85 were closest to a leftover and **none** to a clustered
  member. 331 new profiles across 10 core features took the installed set from
  484 to 916, and recovered 60 of those 115 from zero.

---

## Analyses that are now scripts, not prose

Each of these was performed by hand during the Rhabdoviridae build and has been
written up as a script so the next taxon does not repeat the reasoning:

| Script | Replaces |
|---|---|
| `analyze_feature_gap.py` | the region-extraction / ORF-call / leftover-vs-clustered test that found the P gap |
| `check_cleaved_ends.py` | the end-creep measurement that found the `G_MAT` defect |
| `merge_segments.py` | merging Dichorhavirus records, done by hand for 56 records |
| `minimal_refs.py` | the greedy set cover over unrouted genomes |
| `annotation_rarefaction.py` | the vocabulary-growth curve |

What deliberately stayed prose, because it is judgement and does not transfer:
which features flank the gap (`--left`/`--right`), how many reference contigs
are worth carrying, which termini are genuine protease sites, and whether a
taxon has enough sequences to model at all.

## Measured effect of the changes

Over 1,356 distinct genomes (one exemplar per 95%-identity cluster of the 25,237
Rhabdoviridae records BV-BRC holds between 1 and 50 kb):

| | before | after |
|---|---|---|
| installed profiles | 484 | **955** |
| genomes with all five core proteins | 62% | **69%** |
| genomes with at least three | 83% | **87%** |
| genomes returning no call | 315 | **259** |
| phosphoprotein called, in genomes | 743 | **863** |
| matrix protein called, in genomes | 849 | **938** |
| mature G ending at the precursor stop | 60% | **98%** |

The two proteins that moved most, P and M, are the two the leftover profiles
targeted. `G_mature` call count was unchanged at 974 genomes before and after
the `cleave` correction, so that fix cost nothing and bought 369 genomes a
correct C-terminus.

## A conclusion that was corrected

An earlier version of the coverage report concluded that the routing gap could
not be closed with more reference contigs: that 200 of 315 unrouted genomes
matched no other unrouted genome, and full coverage needed 242 references.

That all-vs-all comparison had been run with settings that under-sampled the
matches -- 7,283 HSP rows where a complete run over a smaller set finds 32,076.
Re-run with the router's own parameters (`-evalue 0.5 -reward 2 -penalty -3
-word_size 11 -soft_masking false`, matching the annotator), the answer inverts:
36 of 259 are singletons, and 43 references would close the gap entirely.

Applying the rule that a reference is only admissible when a module has profiles
for that taxon, 140 of the 259 are eligible and 10 references recover just over
half of them. `minimal_refs.py` enforces the admissibility rule and prints the
curve; where to stop on it remains a human decision.

## Still to decide before pushing

- Whether the clustering patches go in at all, and if so whether as commits to
  `Other_Scripts/` or as a separate branch.
- `FCP_Trunc_Utils.pm` is a new file; the repository has no equivalent.
- 10 features carry model-proposed citations that no human has read.
