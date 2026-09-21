# Changes to Viral_Annotation — Rhabdoviridae, Togaviridae, Matonaviridae and Alsuviricetes module work

Draft. Every entry names the file, what changed, and the measurement or
observation that prompted it. Patches against the unmodified originals are in
this directory as `<file>.patch`.

**Nothing here has been pushed.** The repository's clustering code is still
original: `Other_Scripts/FCP_Ali_Utils.pm`, `FCP_Main_Utils.pm` and
`FCP_Nterm_Utils.pm` on disk match `git HEAD`, and the patched copies live only
here. Three tracked files differ from HEAD, listed under "Already applied".

---

## Alsuviricetes module work — Alphaflexiviridae and the Allexivirus split

### `install_module.py` — retire stale alignments, as we already do for profiles

The installer pruned profiles a rebuild no longer produced but had no
equivalent for alignments, and it copies into the destination without clearing
it. Every alignment written by an earlier build therefore survived a rebuild
that renumbered or dropped its cluster.

Alphaflexiviridae shipped 276 profiles against 289 alignments: thirteen files
from two earlier builds, including an NABP cluster whose six members had since
fallen to the leftover pool. Alignments are provenance rather than runtime
data, so nothing was miscalled — but the counts no longer described the
module, and the same gap is what made Hepeviridae read 15 profiles against 13
alignments. Backfilled: Alphaflexiviridae 276/276, Hepeviridae 63/63,
Tobamovirus 98/98.

### `build_leftover_pssms.py` — integers, not `loN`

Second-pass clusters were named `lo1`, `lo2` … so they could not collide with
the main pass. Collision avoidance never needed a prefix, only a refusal to
restart at 1 — which is what the function's own comment already claimed it
did. The cost was that the id is written into every annotated genome as
`family_assignments -> ["LOWVAN", "<Module>.<FEAT>.<id>", …]`, so a second id
shape existed only in modules this kit built. 684 profiles renumbered across
twelve modules; `renumber_leftover_profiles.py` performs the migration.

That script had a bug worth recording because it is silent: it renamed the
repository's profiles and alignments and the working directory's alignments,
but **not the working directory's `pssms/`**. `install_module.py` installs
from there, so the next install copied the `loN` profiles back and retired the
renumbered ones. Alphaflexiviridae returned with 124 `loN` profiles against
integer-named alignments and six mismatched features.

### Rhabdovirus alignments — 223 were never aligned

The kit ships alignments and not profiles, because profiles are derived and an
order of magnitude larger; `pssms_from_alignments.py` regenerates them. That
contract only holds if every shipped alignment is actually aligned.

223 were not — raw unaligned sequence sets with no gap characters at all, 152
in Alpharhabdovirinae and 71 in Betarhabdovirinae. Those profiles could not be
rebuilt from the kit, which made the two largest rhabdovirus modules unusable
by anyone starting from this repository. Realigned with `mafft --auto
--anysymbol`; verified end to end by staging Alpharhabdovirinae into an empty
working directory and regenerating 529 of 529 profiles.

### `upstream_ext` is a per-feature measurement, not a policy

`scan_to_met_start` permits the annotator to scan upstream for a Met when a
PSSM match does not start with one. Whether that helps is a property of the
feature and the taxon, so it was measured rather than assumed: the annotator
was instrumented across 273 annotated genomes and every firing recorded as
finding a Met or running to the previous in-frame stop.

```
                Allexivirus      Potexvirus
  feature     found / ran-to    found / ran-to    ext
  CP             12 / 2            26 / 6          1 / 1
  TGB2            0 / 0            13 / 4          1 / 1
  ALLEXI_40K     19 / 4             - / -          1 / -
  NABP            0 / 0             1 / 0          1 / 1
  REP             1 / 4             3 / 6          0 / 0
  TGB1            0 / 18            7 / 18         0 / 0
  TGB3            2 / 190          38 / 10         0 / 1
```

TGB3 is the row that forced Allexivirus into its own module: the scan is right
for Potexvirus and catastrophic for Allexivirus, because allexiviruses
initiate TGB3 at a **CUG** rather than an AUG (Lezzhov et al. 2015, J Gen
Virol 96:3159-64, PMID 26296665, by site-directed mutagenesis on shallot virus
X). One feature cannot hold both settings.

Two results here would have been missed by setting the field from policy:
**TGB1 and REP fail in both genera**, and are now 0 in both modules.

Measured effect, same genomes and tools before and after the split: genuinely
clean 227/271 (83.8%) to 249/273 (91.2%); `Feature is too long` 30 to 4; TGB3
calls not starting at a Met 100/266 (37.6%) to 5/164 in Alphaflexiviridae.

### `scan_to_met_start` — reported, deliberately not patched

When no upstream Met exists the function returns the codon after the previous
in-frame stop rather than declining, so the emitted CDS runs to the edge of
the upstream ORF and begins at a non-Met. Both strand loops have it; the call
site's own comment says "this must be turned off if there is a non-AUG start".

A patch making the scan decline was written and measured — all 30 over-long
flags removed, genuinely clean 75.0% to 93.2% on a 148-genome subset, and a
near no-op where profiles are adequate (Potexvirus TGB3 median length 76 to
75). It was **reverted at the curator's direction**; the original behaviour
stands, and the per-feature `upstream_ext` settings above are the containment.

Recorded here because it shapes every coverage figure in the Alsuviricetes
audits, not as a pending change.

---

## Togaviridae module work

Two annotator changes, both of which affect every taxon, not only this one.

### `annotate_by_viral_pssm.pl` — per-feature `internal_stop`, and `-mis`

Most alphaviruses carry an in-frame opal (UGA) **six codons before the end of
nsP3**, read through by roughly 10% of ribosomes. tblastn breaks at a stop by
design, so the match was cropped there and nsP3 came out seven residues short on
eight genomes in nine.

A feature may now declare `internal_stop: 1`, and one `$readthrough` decision
gates both the scan-to-stop and the crop. `-mis N` (default 1) caps how many
stops are tolerated: one readthrough codon is expected, a run of stops is a
broken genome and is still cropped. `keep_stop` remains, unchanged, as the
debug-time override.

The flag is necessary but not sufficient. A profile trained on sequences that
all stop at the opal has never seen what follows it, so that feature's
alignments must also be rebuilt spanning the stop — from the end of the
preceding protein to the start of the following one. After that, Sindbis nsP3
is emitted at 556 aa ending `...SRRTEY*LTGVGG`, which is GenBank's own
readthrough translation exactly.

### `get_transcript_edited_features.pl` — screen the protein before emitting it

The four inclusion gates (identity, coverage, gap count, gap runs) all describe
the **nucleotide** match. Nothing looked at what the gap-fill translated to.
Because the fill takes its bases from the *subject*, an N-masked target yields
an X-bearing protein, and the annotator emitted it: **31 of 293** Togaviridae TF
calls carried an `X` and one carried an internal stop. No amount of reference
curation prevents this — the defect is in the genome being annotated, not in
the reference.

The translation is now screened, and a failure demotes to `partial_cds`, which
is what that branch already existed for. Measured on the taxon: 293 calls, 263
of them clean.

### Skill and tooling

- `install_module.py` now prunes stale profiles **inside** features that still
  exist. `_retire()` only removed whole retired feature directories, so an
  orphaned `.pssm` from an abandoned build stayed installed and callable —
  `Togaviridae.P123.6.pssm` did exactly that.
- `validate_calls.py` names the special features it does not evaluate. It runs
  `annotate_by_viral_pssm.pl` only, so its count is lower than the pipeline's.
- `evaluate_coverage.py` warns when the JSON declares a special feature: its
  flat-table path cannot represent one, so that coverage is absent, not low.
- `skill/references/special-features.md` is new and covers all of the above.

### Known, not fixed

On **NC_043402** the nsP1234 match contains both the opal and the ORF's genuine
terminator, so `n_stops` exceeds `-mis`, readthrough is refused, and the call is
cropped back onto P123's coordinates as a duplicate. `-mis 2` does not fix it —
it reads through the real terminator instead. Distinguishing an opal from a
terminator inside one match is a design decision, not a patch. One genome in 39.

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

---

## Matonaviridae module work

Two skill-script fixes, both applied here, plus one documentation contradiction
that needs a decision. Nothing in `build/` or `annotate/` changed — this taxon
needed no annotator work at all.

### `skill/scripts/collect_registry.py` — `nondefault` tested the wrong thing

```python
("nondefault", os.path.exists(bp)),          # before
("nondefault", _has_departures(bp)),         # after
```

`SKILL.md` instructs writing a `BUILD_PARAMS` file next to **every** feature, and
`run_pipeline.sh` does, containing `departures   none`. Testing for the file's
existence therefore marked every feature in the module as built with non-default
parameters. The published Matonaviridae registry claimed *"non-default build
parameters (7)"* for a module built entirely at the documented defaults, with a
marker beside all seven rows.

`_has_departures()` parses the `departures` line and returns true only for a real
departure. Diff: `patches/collect_registry.nondefault.patch`.

### `skill/scripts/gen_registry.py` — emitted class names its own CSS does not define

The generated HTML wrote `class="th"`, `class="tag"` and `class="tnote"`; the
stylesheet in the same file defines `.taxhd`, `.taxmeta` and `.taxnote`. Every
per-module header, metadata line and prose note rendered unstyled — no serif
heading, no rule, no muted colour. This affected the Rhabdoviridae and
Togaviridae registry pages already in `reports/` (6 module blocks and 1
respectively); the class names in both were corrected in place, since the
`registry.json` they were generated from is not in the repo.

Diff: `patches/gen_registry.classnames.patch`.

### `skill/assets/annotation-vocabulary.tsv` — 7 Matonaviridae rows

`Nonstructural polyprotein` / p200, `Protease p150 protein` / p150,
`RNA-dependent RNA polymerase` / p90, `Structural polyprotein` / p110,
`Nucleocapsid protein` / C, `Mature envelope glycoprotein E2` / E2,
`Mature envelope glycoprotein E1` / E1.

Six of the seven annotation strings already existed; only `Protease p150
protein` is new, patterned on Togaviridae's `Protease nsP2 protein`. Naming
only the demonstrated activity is deliberate — p150's protease is proven while
its methyltransferase is homology-assigned, and p90's polymerase is proven while
its helicase is homology-assigned, so p90 reuses the existing string rather than
spending new vocabulary on an inferred function.

**The `PubMed IDs` column is deliberately empty on all seven rows.** Every
citation in this module is model-proposed, and the vocabulary table has no
provenance column, so a flagged citation entering it would silently become a
curated one.

### Still open: `references/json-schema.md` contradicts `check_dlits.py`

The document says a model-proposed citation is listed in **both** `PMID` and
`PMID_claude_generated`:

> So any citation a model supplied is listed **both** in `PMID` and in
> `PMID_claude_generated`

`check_dlits.py` reports exactly that as a problem — *"CITATION IN BOTH FIELDS
… they must be disjoint"* — and the shipped Togaviridae module **is** disjoint:
`PMID` holds the one human-read DLIT and `PMID_claude_generated` the proposals,
zero overlap. Code and shipped data agree, so the document is the outlier.
Matonaviridae follows the code. The doc is unchanged pending a decision on which
convention is intended.

### Not applied, but measured and written up

`example-matonaviridae/README.md` records six further kit issues this build hit,
three of which fail silently while reporting success: `-download-only` shells
out to `BVBRC_clean_contigs.pl`, which the kit does not ship, and downloads
nothing with exit code 0; `make_gto.py` accepts only `.fna`/`.fasta`/`.fa` and
so silently skipped all 217 genomes written as `.contigs` by the kit's own
downloader; and `run_gto_eval.py` runs with `cwd=wd` while passing `--gto-dir`
through unchanged, so a relative path fails for every genome. None is patched
here — each is a small change with a design choice attached.


---

## Alsuviricetes class build — Hepeviridae

First module of a class-wide build (26 modules planned, partitioned by genome
organisation). Two annotator/pipeline changes here affect **every** taxon.

### `Other_Scripts/FCP_Ali_Utils.pm` — `-fd` measured against the wrong denominator

Patch: `patches/FCP_Ali_Utils.frac_dash_denominator.patch`

`process_alignment` removed sequences whose dash fraction exceeded `-fd`
(default 0.20) **before** excluding columns below the `-f` occupancy floor, and
measured that fraction over every column MAFFT emitted. Since the `-f` filter
exists precisely to discard columns only a few sequences occupy, an ordinary
sequence was scored as mostly gaps because *somebody else* had an insertion.

Hepeviridae ORF1 found it:

```
alignment                   967 sequences x 2363 columns
ungapped length             min 1483  median 1704  max 1790
dash fraction per sequence  min 0.242 median 0.279 max 0.372   (-fd cut is 0.20)
sequences under the cut     0 of 967
columns below 5% occupancy  615, longest runs 93, 77, 66, 57, 54, 51, 45, 41
```

All 967 were deleted, the alignment went to zero rows, and
`print_alignment_as_fasta` died — **0 PSSMs from 1,428 sequences, and the run
still exited 0**. The other three features had already been written, so it read
as a partial success. Fix: score each sequence over the columns that will
survive the `-f` exclusion. Only the denominator changes. After it:
`1428 seqs -> 8 alignments, 8 PSSMs`, cluster 1 keeping 908 of 967 with 654
low-occupancy columns cut.

This will bite any taxon with a variable-length region inside a conserved
protein — in Alsuviricetes, most of them.

### `annotate_by_viral_pssm.pl` — a genome that routes but matches nothing

Patch: `patches/annotate_by_viral_pssm.empty_result.patch`

Three defects compounding, all documented in `references/pipeline.md` as bug 7
but not previously shipped as a patch. `write_fasta` confesses on an empty
list, so a genome that routes to a module and clears no profile's threshold
produced **no output at all and exit 255** — indistinguishable from a crash, and
in a batch run it silently becomes a coverage gap. Underneath,
`unless ($hsp_best_bit < $best_bit)` is true at `0 == 0`, so a profile matching
nothing became `$best_pssm` and the **last** profile in the directory won.

Seen on Fish-associated hepevirus (OP933684), the only near-complete genome of
its species in BV-BRC. Fix: guard both writers with `if @gene_seqs` /
`if @prot_seqs`, use `>` for profile selection, and `next` past a feature with
no result rather than dereferencing undef.

### `skill/scripts/` fixes

| script | defect |
|---|---|
| `build_leftover_pssms.py` | wrote **unaligned** sequences into `corrected_alis/` — the directory SKILL.md calls the curation contract. 12 of 38 were then unrebuildable (`psiblast: Repeated Seq-IDs`, which is misleading — the ids are unique, the rows are ragged); the other 26 had coincidentally equal-length rows and would have been silently rebuilt **wrong**. Also `--mi` defaulted to 0.6, silently lowering the identity floor below the feature's own. |
| `run_gto_eval.py` | neither `--gto-dir` nor `--out`/`--report` was absolutised while the subprocesses run with `cwd=<out>/_wd`. **669 of 669 genomes failed**, the run exited 0, and `quality.tsv` was written containing only a header. |
| `collect_synmap.py` | `Rhabdoviridae_Viral_PSSM.json` was hardcoded — it could only ever work for the family it was written for. |
| `collect_registry.py` | counted `.outliers` in the self-recall denominator. On a partial-submission-heavy family that is fatal to the number: Hepeviridae ORF1 read 27.1% against 5,350 sequences instead of 97.5% against its 1,429 in-range ones, because 3,921 fragments — some 8 aa — cannot clear an 800-bit cutoff and were never meant to. Now reports both. |
| `gen_registry.py` | headlined the psiblast proxy and pooled all features. Now prefers the annotator's own genome-level rate and headlines the **median** per-feature rate: pooling averaged three core proteins with a Rocahepevirus-only accessory and produced a meaningless "76%". |

### `PMID_claude_generated` — the documentation was wrong

`SKILL.md` step 7 said a proposed citation "goes in `PMID` and in
`PMID_claude_generated` alongside it", and `references/json-schema.md` described
the same nested form with a worked example. Both contradict `check_dlits.py`,
which reports an id in both fields as an error, and contradict all 38 installed
modules, which use the disjoint form.

Worse, the phrase "a citation you did not read is not a DLIT" reads as implying
that one you *did* read is — and that is how this build initially put four ids
in `PMID`. Rewritten in both files: every model-supplied citation goes in
`PMID_claude_generated` and nowhere else, the lists are disjoint, the test is
**who supplied the id** rather than how much work went into it, and reading a
paper does not clear the flag — only a curator deleting the registry entry does.
Also added: never put a flagged id in `annotation-vocabulary.tsv`, which has no
provenance column.

### Not fixed — `coverage_cutoff` does not measure coverage

Reported, deliberately **not** patched, because it changes calling behaviour for
every installed module and is the project owner's decision.

`matching_tblastn_hsps_json` computes
`length(hseq) / (q_to - q_from + 1)` — the aligned hit length over the aligned
**query span**, which is the alignment measured against itself and is 1.0 for
any gapless local alignment however little of the profile took part. An HSP
covering 20% of a profile passes `coverage_cutoff 0.65` trivially; the only
thing the test can reject is a gappy alignment. `references/json-schema.md`
describes it as a "subject coverage floor".

Found on Endornaviridae: one profile produced three HSPs against a single
19,406 nt genome (666, 496, 215 bits), each extended to a start and stop, so one
polyprotein became three nested CDS at 521-19342, 7316-19342 and 15617-19342.
`evaluate_module.py` reports **0 duplicates** for this, because it compares
exact coordinates. A real test would be `(length(hseq) - gaps) / profile_length`.

### New in `skill/scripts/`

`gen_rarefaction.py` promoted from `example-matonaviridae/` and parameterised
(`--taxon`, `--rarefaction`, `--out`); it previously resolved its input relative
to its own location. Plus eight programs this build needed that are not
taxon-specific: `fetch_features.py`, `fetch_seqs.py`, `fetch_contigs.py`
(batched BV-BRC dumps — `query_PATRIC_bob.pl` issues one HTTP round trip per
input line and rebuilt `P3DataAPI` each time, which projected to ~4 h for 87,179
genomes against ~3 min batched), `subset_dump.py`, `build_features.py`,
`rescue_unassigned.py`, `reroute_outliers.py`, `qc_truncation_symmetric.py`.

The last of those exists because the pipeline's own truncation QC tests
containment in one direction only: Hepeviridae ORF3 cluster 8 is 75 aa and
**100% identical to residues 25-99** of a 122-aa ORF3, and it passed because the
row with cluster 8 as BLAST query aligned 69 of 75 residues (qcov 0.92, under
the 0.95 cut) while the reverse row aligned all 75. It also covers
`reclustered_alis/`, which the pipeline's QC never sees.
