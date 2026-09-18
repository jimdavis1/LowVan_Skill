# LowVan kit feedback -- Alsuviricetes build

## check_dump.py counts a blank sequence as present

`fetch_seqs`/`query_PATRIC_bob.pl` emit `md5<TAB>` with an empty second field when
`feature_sequence` has no record for an md5. `check_dump.py` joins on the md5 and
reports "every md5 named by id_md5 has a sequence", so the hole is invisible:

```
--- completeness ---
  distinct md5 named by id_md5   : 69231
  distinct md5 with a sequence   : 69231
  ok    every md5 named by id_md5 has a sequence
```

71 of those 69231 rows are empty. Suggest counting rows whose sequence field is
zero length and reporting them separately, the way the blank-md5 note already
works. Found on Alsuviricetes, 17 September 2026.

## query_PATRIC_bob.pl is one HTTP round trip per input line

It also constructs a new `P3DataAPI` per line. On 87,179 genomes that projected to
about four hours for `id_md5` alone. Replacing it with batched `in(genome_id,(...))`
POSTs to `/api/genome_feature/` over a 12-way thread pool took 3 minutes for the
same query -- roughly 80x. Scripts in `work/fetch_features.py` and
`work/fetch_seqs.py`; they also return coordinates, which the gene-order analysis
that drives module partitioning needs anyway.

## references/pipeline.md documents a fix that assets/patches does not contain

Bug 5 in `references/pipeline.md` ("Trim caps -- the guards invert exactly when
they are needed") says the fix is:

> Replace "abort if occupancy is low" with "trim, but never more than a fixed
> fraction of the alignment" ... Exposed as `-max_end_frac` and
> `-max_start_frac` in `FCP_Main_Utils.pm`, both defaulting to 0.35.

Neither option exists in `assets/patches/FCP_Ali_Utils.patch` or
`FCP_Main_Utils.patch`, and after applying all three patches cleanly
`grep -c max_end_frac Other_Scripts/*.pm` is 0. The `n_occ` option from the same
section *is* in the patch, so the omission looks accidental rather than
deliberate. Either ship the cap or drop the paragraph -- as it stands a reader
follows the reference, greps for the option, and cannot tell whether their
patch failed.

Also worth noting for whoever applies these: the three `.patch` files are
unified diffs with `a/Other_Scripts/` prefixes, so they need `patch -p1` from
the **repo root**, while `FCP_Nterm_Utils.patch` is an ed-style diff that has to
be applied to the file directly. A one-line note in the reference would save
the guesswork.

Found on Alsuviricetes, 17 September 2026.

## An eighth pipeline bug: -fd is measured against columns that are about to be deleted

**This one silently destroyed the largest feature of the first module built.**
Patch in `work/patches/FCP_Ali_Utils.frac_dash_denominator.patch`.

`process_alignment` in `FCP_Ali_Utils.pm` does two things in this order:

1. remove sequences whose dash fraction exceeds `-fd` (default 0.20),
   measured over **every column MAFFT emitted**;
2. exclude columns whose occupancy is below `-f` (default 0.33).

Step 2 exists precisely because an insertion carried by a few sequences gets its
own alignment columns and everyone else gets a dash. Measuring step 1 against
those columns scores an ordinary sequence as mostly gaps because *somebody else*
has an insertion.

Hepeviridae ORF1 found it. 967 sequences, 1483-1790 aa, tightly homogeneous:

```
alignment                   967 sequences x 2363 columns
ungapped length             min 1483  median 1704  max 1790
dash fraction per sequence  min 0.242 median 0.279 max 0.372   (-fd cut is 0.20)
sequences under the cut     0 of 967
columns below 5% occupancy  615, longest runs 93, 77, 66, 57, 54, 51, 45, 41
```

So all 967 were deleted, `@ali3` went empty, `$ali_len` came out undefined, six
"Use of uninitialized value" warnings scrolled past, and then:

```
Bad sequence entry passed to print_alignment_as_fasta
 at gjoseqlib.pm line 956
```

which aborted the whole feature. **ORF1 produced 0 PSSMs from 1428 sequences and
the run still exited 0.** The other three features had already been written, so
the summary looked like a partial success rather than a failure. Nothing in the
output says "the most important protein in this module cannot be called".

Note this is the same empty-alignment crash as documented bug 7, but reached
from the build side rather than the annotator side, and guarding the writer
would only have converted a silent catastrophe into a silent empty feature. The
denominator is the actual defect.

**Fix**: compute each sequence's dash fraction over the columns that will
survive the `-f` exclusion. Only the denominator changes -- a sequence that is
genuinely fragmentary relative to the conserved core is still removed, which is
what `-fd` is for. After the fix:

```
ORF1   1428 seqs -> 8 alignments, 8 PSSMs     (was 1 alignment, 0 PSSMs)
  cluster 1: 967 -> 908 kept, 654 low-occupancy columns cut, Ali_Len 1709, starts at M
```

This will bite any taxon with a variable-length region inside a conserved
protein, which in Alsuviricetes is most of them -- the HEV ORF1 hypervariable
and poly-proline regions here, and the same shape is visible in the
Closteroviridae and Betaflexiviridae collections.

**Suggestion beyond the patch**: make `fasta-cluster-pssm-2.pl` exit non-zero,
or at least print a clear final-line warning, when a feature finishes with zero
PSSMs. A silent zero is indistinguishable from a feature that was not requested.

## build_leftover_pssms.py --mi defaults to 0.6, contradicting both docs

SKILL.md, step 3: "It re-clusters the leftovers on their own at the **same**
`-mi` -- the identity floor is never lowered". The script's own docstring agrees:
"--mi the feature was built at". But `ap.add_argument("--mi", type=float,
default=0.6)`, so the default silently lowers the identity floor for any feature
built at the documented default of 0.8.

That is exactly the failure mode the `-mi` floor discussion in
`references/pipeline.md` warns about, arriving through the back door: a profile
built at 0.6 from sequences that failed to cluster at 0.8 is a wider clade than
the feature's other profiles, and nothing in the output says the threshold
differed.

Either default it to the feature's build `-mi` by reading `BUILD_PARAMS`, which
is sitting right next to the alignments, or make `--mi` required.

Found on Alsuviricetes Hepeviridae, 17 September 2026.

## build_leftover_pssms.py writes UNALIGNED sequences into corrected_alis/

Patched in both copies; `.bak` alongside. This is the one I would fix first.

`build_pssm()` aligns its group internally with MAFFT and builds the profile
from that alignment -- correctly. But the caller then wrote its **own
unaligned input** to `corrected_alis/<cid>.fa`:

```python
with open(os.path.join(fd, "corrected_alis", cid + ".fa"), "w") as fh:
    for h, s in sub:
        fh.write(">%s\n%s\n" % (h, s.replace("-", "")))
```

`sub` comes from `Leftover_Seqs.aa`, which is the unaligned file. So the
profiles are right and the files on disk are not the alignments they came from
-- in the directory SKILL.md calls "the contract: PSSMs are built from it, and
it is the file you hand-edit".

Two consequences, and the quiet one is worse:

- **Loud:** `rebuild_pssms.py` fails on any such file whose rows differ in
  length. On Hepeviridae that was 12 of 38, reported as
  `psiblast produced nothing`. The real message from psiblast is
  `Repeated Seq-IDs detected in multiple sequence alignment file`, which is
  misleading -- the ids are unique; the rows are ragged. Renaming every row to
  `s0..s3` reproduces it exactly.
- **Quiet:** the other 26 files have rows that happen to be equal length after
  gap-stripping, so they look like valid alignments. `rebuild_pssms.py --write`
  will happily rebuild a profile from unaligned sequences and report success.
  Any leftover group whose members genuinely need gaps gets a silently wrong
  profile on the next rebuild.

**Fix**: have `build_pssm` expose the MSA it used and write that, with the
original headers restored. After the patch, Hepeviridae goes from
`51 reproduce / 0 stale / 12 failed` to `63 reproduce / 0 stale / 0 failed`.

## rebuild_pssms.py cannot reproduce a freshly built module

"`rebuild_pssms.py` reports 0 stale" is in SKILL.md's *What "done" looks like*,
but it cannot be true immediately after a pipeline run. The pipeline builds
profiles through `BlastInterface::alignment_to_pssm` plus
`FCP_PSSM_Utils::process_pssm_file`; `rebuild_pssms.py` calls `psiblast`
directly and pipes through `norm_pssm.py`. The score matrices are identical --
the two differ only in the query-id block:

```
pipeline   local str "2"   +   descr { title "2 Untitled PSSM 2" }
rebuild    local str "1"       (no descr block)
```

57 bytes in a 65 kB file. Harmless at runtime, because
`annotate_by_viral_pssm.pl` selects profiles by **filename** and never reads the
embedded title. But a byte comparison calls all 25 of them STALE on a module
nobody has touched, which trains the reader to ignore the check -- exactly when
it is supposed to be catching real drift.

Suggest either normalising the pipeline's output through `norm_pssm.py` too, or
having `rebuild_pssms.py` compare the score matrix rather than the whole file
and report query-block differences separately.

Found on Alsuviricetes Hepeviridae, 17 September 2026.

## run_gto_eval.py: --gto-dir is not absolutised, so every genome fails

Confirmed the relative-path bug already listed in the Matonaviridae feedback,
with a diagnosis and a one-line patch. Applied to both copies, `.bak` alongside.

`main()` absolutises `--repo` (line 42) but not `--gto-dir`, and the annotator
is then invoked with `cwd=<out>/_wd` (lines 69, 76). So a relative `--gto-dir`
resolves against the wrong directory:

```
Could not open input file gto/1016879.45.gto: No such file or directory
  at /Applications/BV-BRC.app/deployment/lib/SeedUtils.pm line 2868.
```

On Hepeviridae that was **669 of 669 genomes**, and the run still exited 0 with:

```
  scored 0 genome(s)
  669 annotation failure(s), 0 quality failure(s)
  wrote quality.tsv
```

`quality.tsv` was written containing only a header. A reader who trusts the exit
code and the "wrote quality.tsv" line gets an empty quality section and no
indication that nothing was scored. Suggest also exiting non-zero, or at least
printing a warning, when the scored count is zero but GTOs were supplied.

Fix: `glob.glob(os.path.join(os.path.abspath(args.gto_dir), "*.gto"))`.

## coverage_cutoff does not measure profile coverage, and is effectively inert

`annotate_by_viral_pssm.pl`, `matching_tblastn_hsps_json`:

```perl
my $cov = ((length $results->{hseq})/(abs($results->{q_to} - $results->{q_from})+1));
if (($results->{bit} > $bit_cutoff) && ($cov > $cov_cutoff)) { push @data, $results; }
```

The denominator is the **aligned query span**, not the profile length. So the
ratio is the alignment measured against itself: for any gapless local alignment
it is 1.0 no matter how small a fraction of the profile took part. An HSP
covering 20% of a profile passes `coverage_cutoff 0.65` trivially. The only
thing the test can actually reject is a gappy alignment.

`references/json-schema.md` describes the field as "subject coverage floor,
typically 0.65", which is not what it does. Every module in the repo sets it to
0.65 and none of them is getting the filter they think they are.

Found on Endornaviridae, whose polyproteins are 3,070-7,156 aa. One profile
produced three separate HSPs against a single 19,406 nt genome:

```
Choosing Endornaviridae.POLY.lo8.pssm   665.6
  KX355144   16040-19339   frame 2   665.6
  KX355144     731- 4969   frame 2   495.8
  KX355144    7463- 8848   frame 2   214.6
```

Each was then extended to a start and a stop, so one polyprotein became three
CDS features at 521-19342, 7316-19342 and 15617-19342 -- nested calls sharing a
C-terminus. On a 24-genome panel this turned 26 reference CDS into 35 calls.

Note `evaluate_module.py` reports **0 duplicates** for this, because it compares
exact coordinates and these differ. So the current QC cannot see it either.

A real coverage test would be
`(length(hseq) - gaps) / <profile length>`, with the profile length read from
the PSSM. I have **not** changed this -- it alters calling behaviour for all 38
installed modules and is a decision for the project owner, not a patch to slip
into one build. Raising `bit_cutoff` is not a substitute: the spurious HSPs here
scored 215 and 496 against true signals that go as low as 84.

Suggest also teaching `evaluate_module.py` to report overlapping/nested calls of
the same feature, not only coordinate-identical ones.

---

## CORRECTIONS to earlier entries in this file

Three things I reported above were wrong, and the cause was the same in each
case: I cloned `CEPI-dxkb/Viral_Annotation` fresh and treated it as the
baseline, when the kit's own `build/` and `annotate/` directories are further
ahead. The upstream repo has not received the Togaviridae or Matonaviridae
annotator work; the kit has.

**1. "references/pipeline.md documents a fix that assets/patches does not
contain" — withdrawn.** `-max_end_frac` and `-max_start_frac` are real and
shipped: `build/FCP_Ali_Utils.pm` and `build/FCP_Main_Utils.pm` both carry
them. They are absent only from `patches/`, which is a diff against the
unmodified upstream and is therefore always behind `build/`. The reference is
correct.

**2. "annotate_by_viral_pssm.pl dies when a genome matches nothing" — already
fixed.** `annotate/annotate_by_viral_pssm.pl` already had both writer guards
and `if ($hsp_best_bit > $best_bit)`. I rediscovered a solved problem because I
was reading the upstream clone. Only one line of that entry is new: a guard so
a feature whose profiles all score zero does `next` instead of dereferencing an
undef result set. That is now in the kit's copy.

**3. The `internal_stop` claim needs splitting in two.** The *modules* are fine
— Togaviridae annotates CHIKV S27 completely and correctly, all 13 features,
mature peptides tiling their precursors, `nsP123` and `nsP3` both ending at
5665 as intended. What is true is narrower: **`CEPI-dxkb/Viral_Annotation` on
GitHub has no `internal_stop` support at all**, so anyone cloning that repo gets
alphavirus nsP3 cropped at the opal on isolates that carry one. That is a
distribution gap, not a defect in Togaviridae.

**The general lesson, which belongs in the kit:** `patches/` is a diff against
pristine upstream and drifts behind `build/` and `annotate/` as work
accumulates. **The kit's own directories are the source of truth, not the
patches and not the upstream repo.** A new build should copy `build/*` and
`annotate/*` over its `Viral_Annotation` checkout before starting, rather than
applying `patches/` to a fresh clone — which is what I did, and it cost a
rediscovered bug, a wrong feedback entry, and a Tobamovirus build run against
an older pipeline than the kit ships.
