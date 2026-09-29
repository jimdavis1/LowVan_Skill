# Silent failures

Every entry here produced **no error message**. The build completed, the
numbers looked plausible, and something was wrong. They are collected in one
place because the pattern matters more than any individual bug: in this kit,
the dangerous failures do not crash — they return a confident wrong number.

The general defence is the same each time. **Check the thing you actually
care about, not a proxy for it.** Zero quality failures does not mean the
features were called. A validated dump does not mean the annotations line up
with the sequences. A profile that recalls its own training data does not mean
it calls anything in the wild.

---

## The JSON

### `copy_num` is what makes a feature essential

Nothing in the JSON is named "essential". `viral_genome_quality.pl` treats any
CDS or `mat_peptide` carrying `copy_num` as essential and flags it missing on
every genome that lacks it. A feature without `copy_num` is still called and
still reported — it just never counts against the verdict.

Quinvirinae declares three accessory proteins. ORF5a and ORF2a exist only in
*Robigovirus*, p14 only in *Foveavirus*. Declared with `copy_num: 1` they
produced **258 "missing essential feature" flags** and drove genuinely-clean to
**0.0%** on a module whose core features call at 98–100%. Removing `copy_num`
from those three: **79.6%**.

**The first version of this rule was wrong, and the bug recurred because of
it.** After Quinvirinae I wrote "count the taxa in the feature's collection"
and automated that as *single genus → everything essential*. It promptly
happened again in two more modules, because **accessory ORFs vary by species
WITHIN a genus**: Closterovirus p13 is carried by 9.5% of genomes and
Ampelovirus p20A/p20B by 6.1%, both in single-genus modules. Both scored
**0.0% genuinely clean** while their core features called at 91–100%.

**Set `copy_num` on the fraction of GENOMES that carry the feature, not on how
many genera it spans.** Measure it on the panel:

```python
essential = called[k] / routed >= 0.85
```

That de-escalated 7 of Closterovirus's 13 features (0.0% → 33.8%) and 6 of
Ampelovirus's 11 (0.0% → 52.5%).

Readthrough products and anything variable in presence get none either.

Note the distinction that matters. A feature that is *absent* in some genera
should be non-essential. A feature that is *present but uncallable* — like
Pestiviridae NS5A in the atypical lineages — should stay essential, because
making it optional hides a real gap.

### `segments` takes an object, not a scalar

```json
"segments": {"Single RNA Segment": {"min_len": 8000, "max_len": 11000,
                                    "replicon_geometry": "linear"}}
```

Written as `{"Single RNA Segment": 1}` the module annotates perfectly and then
dies in quality scoring on **every genome** with

```
Can't use string ("1") as a HASH ref ... viral_genome_quality.pl line 302
```

so the module looks built and scores nothing. `install_module.py --check` now
validates the shape; it did not until this bit.

### `min_len`/`max_len` on a derived feature are NUCLEOTIDES

The annotator gates a `non_pssm_partner` feature on `abs($begin - $end)`, which
is a nucleotide span. Pestiviridae NS2 is 457 aa; declared as `380–560` it was
silently rejected on every genome, because the real span is ~1,371 nt. The
Arenaviridae precedent is the one to copy — a ~230 aa Gn with a 450–900 window.

`end_offset` is **subtracted**: to land one base before the anchor's start, use
`+1`, not `-1`.

---

## The dump

### The three `uniq.*` files arrive misaligned

They are line-index aligned and that alignment is the join key. BV-BRC has
delivered them misaligned on **every taxon this project has dumped**:

| taxon | md5 / id_ann / seq | sequence-less md5 |
|---|---|---|
| Hepaciviridae | 100,855 / 103,166 / 100,290 | 565 |
| Pestiviridae | 11,318 / 11,318 / 11,311 | 7 |
| Pegivirus | 2,767 / 2,767 / 2,764 | 3 |

Unrepaired, every index is off by one from the first gap onward and **the wrong
annotation attaches to the wrong sequence**, silently. Run `check_dump.py`
first and `realign_dump.py --write` if it fails.

### Duplicate contig ids lose whole genomes

BV-BRC sometimes returns the same sequence record twice — 3 of 2,006
Pestivirus genomes, 3 of 670 Pegivirus. `rast-create-genome` then dies with
`Attempt to add duplicate contig id` and **loses the entire genome**, not the
duplicate. `fetch_contigs.py` now deduplicates on accession.

---

## The measurement

### A staleness check that cries wolf gets ignored

`rebuild_pssms.py` reported **399 of 808** Hepaciviridae PSSMs as stale when
every score was identical. Three artefacts of the comparison, none of them a
score: psiblast does not preserve the query id, the comparison was asymmetric
(normalised rebuild against raw stored file), and some pipeline runs emit a
`descr` block and some do not.

Worth recording *how the diagnosis went wrong*, because the failure mode is
general. Diffing **one** flagged profile showed only the `descr` block, and all
399 flagged carried it while all 409 clean ones did not — a perfect
correlation. It was coincidental: the block tracks which pipeline run produced
a profile, not why the comparison failed. Fixing it moved 399 to 389.
**Classifying the diffs across a 30-profile sample gave the real answer in one
line: 16 identical, 14 differing in the query id and nothing else.**

### `quality.tsv` could not reproduce its own headline

"Genuinely clean" is no genome, contig **or feature** flag. The table shipped
only the first two, so reading it gave the "good" figure instead — 65.3%
against a genuinely-clean 53.5% on the same panel. `run_gto_eval.py` now writes
`n_feature_flags`.

### A filter that encodes an assumption will appear to confirm it

Hepaciviridae's collections were filtered on `CYS = set("CT")` and
`SERALA = set("SA")` before the profiles were built, dropping 8,807 sequences.
Measuring the surviving termini then "confirmed" the cleavage chemistry — of
course it did. That is the filter reporting its own admission criterion.

Only two results from that check were independent, and both were strong: the
filter allowed **either** C or T at all four NS3-4A sites and never said which,
and the data chose T at the cis site and C at the three trans sites; and NS2,
whose junction entry was `(None, None)`, came out 98% Leu unprompted.

**When you build a chemistry table, measure the termini BEFORE you filter on
them, and say in the audit which figures are independent.**

---

### A script that has never run once looks exactly like a script with nothing to do

`rescue_unassigned.py` had never executed successfully in this project. Every
Alsuviricetes dump crashed it on line 43:

    ValueError: dictionary update sequence element #0 has length 3; 2 is required

BV-BRC writes `uniq.seq` and `uniq.id_ann` with a TRAILING TAB, so every line
splits into three fields and `dict(...)` refuses it. Six modules -- Carlavirus,
Quinvirinae, Closterovirus, Ampelovirus, Velarivirus, Bromoviridae -- were
built and shipped without a homology rescue, and nothing anywhere said so,
because a rescue that never ran leaves the same evidence as a rescue that
found nothing: no `RESCUED` artifact.

Behind the crash sat a worse bug that the crash was hiding. `assigned` -- the
set of sequences already binned, which the rescue must skip -- was built by
taking the last pipe-delimited field of each collection FASTA header and
treating it as an md5:

    assigned.add(line.rstrip("\n").rsplit("|", 1)[-1])

but `build_collections.py` writes the FEATURE ID as the header
(`>fig|28347.83.peg.1`), so that expression yields `28347.83.peg.1` and never
matches an md5. The set was permanently empty. **Had the trailing-tab crash
been fixed on its own -- the obvious one-line fix -- the rescue would have
re-adopted every already-binned sequence and silently duplicated the entire
collection set.** The crash was the only thing preventing a much quieter
corruption.

Three lessons, in order of how much they cost:

- A tool that has never produced output is not the same as a tool that found
  nothing. **Check that each step in a pipeline has actually run at least
  once**, by looking for its artifact, not by its silence.
- When a script dies on input parsing, fix the parse *and then read the rest
  of it*. The first error can be the only thing standing between you and a
  defect that produces plausible output.
- Assertions about identity -- "the md5 is the last field" -- decay when the
  producer changes format. The comment asserting it was still there, still
  confidently wrong.

### A gate figure and the shipped reference set must agree, or one of them is wrong

`repcontig_budget.py --merge-by-name` reported Crinivirus routing at **100.0%**
at budget 25. `build_rep_contigs.py`, run with the same budget on the same
data, covered 511 of 553 records = **92.4%**. The skill already says these two
must reproduce each other, so the disagreement was the finding.

Cause, and it is the second instance of it: `--merge-by-name` groups records by
genome name to collapse the segments of one multipartite genome. But 18
Crinivirus names carry between 3 and 8 records -- those are ISOLATE
COLLISIONS, several independent isolates deposited under one name, not
segments. Merging them counts several genomes as one and inflates coverage.
Bromoviridae had already been corrected for exactly this (94.8% -> 69.5%);
Crinivirus had not.

**For a multipartite taxon, never accept a merged coverage figure without
first checking the distribution of records per name.** Exactly `n` records for
an `n`-segment virus is a segment set; more than `n` is a collision.
`merge_segments.py` gets this right -- it refuses to merge a name whose record
count exceeds the segment count, and logs it as `collision_not_merged` -- so
the two tools disagreed with each other as well as with the truth.

### One record can switch off a length check for the whole taxon

Crinivirus segment windows are computed from observed contig lengths. One
record -- "Plant associated crinivirus 1 MIR20SW", 16,522 nt -- is 65% longer
than the next largest and almost certainly both segments on a single contig.
Taken at face value it set RNA2's `max_len` to 17,348, which is not a loose
bound but **no bound at all** for the other 276 RNA2 records, all of which sit
between 6,575 and 8,672 nt.

A window derived from extremes inherits every outlier's opinion. Bound it
robustly -- discard anything beyond 1.5x the 95th percentile -- and say in the
build notes which records that excluded. The cost of excluding a real genome
from a *window* is that it may fail a length check; the cost of including a
bad one is that nothing ever fails again.

### A feature can fail to build and leave a directory named after a hash

Crinivirus P6A -- 14 sequences of 53-55 aa -- produced zero alignments at the
default `-m 5`, because no cluster reached five members. The pipeline still
created an output directory, named `85c1dab168` after a content hash rather
than after the feature, containing empty `alis/` and `pssms/`.

`build_features.py` caught it and said so:

    P6A            FAILED -- no output directory carries P6A profiles (1 new dirs)

which is the behaviour you want, but note what the failure looks like from one
step further away: 22 of 23 features built, a full run, no non-zero exit, and
a stray hash-named directory that the next glob would happily walk into.
**Count the features you declared against the features that built**, every
time. `-m 2` fixed it.

### `--jobs N` is not N processes, and an over-subscribed box looks like a hang

`run_gto_eval.py --jobs 8` on a 10-core machine drove the load average to
**138** and completed **zero** genomes in two minutes. It was not deadlocked
and nothing had crashed — each annotator spawns roughly ten short-lived BLAST
children, so eight jobs is about eighty runnable processes, and every one of
them was getting a tenth of a core.

Three things made this hard to read, and all three are worth knowing:

- **Load average lags.** It is a decaying one-minute average, so it keeps
  climbing for a minute after you fix the cause and keeps falling for a
  minute after you make things worse. Do not tune against it; count finished
  genomes over a fixed interval instead.
- **`pkill -f run_eval.sh` does not kill the work.** The driver dies, the
  perl annotators and their BLAST children are reparented and keep running.
  Two restarts stacked three generations of the same run against each other,
  which is what actually produced the 138. Kill the whole family explicitly —
  driver, `run_gto_eval.py`, `annotate_by_viral_pssm`,
  `viral_genome_quality`, then `blastn`/`psiblast`/`tblastn` by name — and
  confirm the count is zero before restarting.
- **A restart re-scores before it advances.** Resume skips a genome only when
  its `.qual.gto` exists, so the first minutes after a restart rebuild
  quality for genomes that already had an `.ann.gto`. Polling in that window
  shows a flat count and reads like a hang.

Rule of thumb: **`--jobs` at about a third of the core count** for this
pipeline, and measure the rate over two minutes before deciding it is wrong.

**Never run two evaluations at once.** This was written down as a lesson and
then re-tripped three times in the same session — once by restarting without
killing the children, once by starting a re-check beside a backfill, once by
doing it again an hour later. Documentation did not fix it, so
`run_gto_eval.py` now refuses to start when `annotate_by_viral_pssm.pl` is
already running (`LOWVAN_ALLOW_CONCURRENT_EVAL=1` overrides). Serialise and
chain instead:

```bash
while ! grep -q "DONE" first.log; do sleep 30; done; ./second.sh
```

A lesson that only exists in prose gets re-learned. If the mistake is cheap to
detect, make the tool detect it.
Reading PSSMs from a cloud-synced directory costs real time too — copying the
runtime to local disk first takes nine seconds and is always worth it.

### The budget estimate and the run disagree in both directions

Two modules in the same batch showed the reference-budget figure missing the
annotation run, in opposite directions, and both are worth knowing because
only one of them is a defect.

**Tymovirus: budget 92.2%, run 95.1%.** The budget picks references by
clustering and counts a genome covered if it falls in a covered cluster.
BLASTn at annotation time is more permissive than that clustering, so three
genomes routed that the budget had written off. Not a defect — the budget is
simply conservative, and **the conservative number is the one to plan with**
and to put in the coverage audit.

**Crinivirus: gate 100.0%, shipped set 92.4%.** Here the estimate was
*optimistic*, and that always means something is wrong. `--merge-by-name`
had collapsed 18 names carrying 3-8 records each — isolate collisions, not
segments — so several genomes were being counted as one.

The rule that separates them: **an estimate below the run is tolerable, an
estimate above it is a bug.** If the budget claims more coverage than the
annotator delivers, find out why before shipping; if it claims less, record
the budget figure and move on.

### A call rate of 0% can mean the feature was never offered a genome

Crinivirus P22A: 26 sequences, a profile that recovers 26 of 26 of them, and
**called on 0.0% of routed genomes**. Nothing was broken. Every one of its 26
source records sits on a genome record shorter than 6 kb, and the evaluation
population is genomes at or above 6 kb, so the feature was never put in front
of a genome that could carry it.

The same effect, in weaker form, across the module:

    feature   source records <6 kb   called on
    P22A          26 of 26             0.0%
    RNASE3        41 of 57            48.7%
    P7A           35 of 45            55.3%
    P22C           0 of 11            13.3%      <- genuinely rare
    P25            0 of  8             5.3%      <- genuinely rare

The two populations are drawn differently and always have been:

  collections    every protein record in the dump, whatever its genome length
  evaluation     only genomes long enough to annotate, >= 6 kb

So a low call rate has two possible causes that look identical in the report —
the protein is rare, or the protein's genomes are partial. **Before reading a
low call rate as a weak profile, check the length distribution of the source
records.** One query separates them.

This does not endanger the `copy_num` decision, and that is worth stating:
both causes argue for the same answer. A feature that cannot be called on
whole genomes must not be marked essential either way, so the safe direction
is preserved. What it endangers is the *explanation* written next to the
number in the coverage audit.

### Two profiles on one ORF: the defect that raises a feature's call rate

Crinivirus shipped 23 features. Five of them were called on 96-100% of routed
genomes and duly earned `copy_num`. They had earned nothing: they were firing
on a locus another feature already occupied.

    RNASE3  inside SUPPRESSOR   97% of its calls
    P9      inside HSP90        98%
    P8A     inside HSP90        95%
    P26     inside P28          83%
    P23     inside P22B         82%
    HSP70 + HSP90 on IDENTICAL coordinates in 7% of genomes

**This raises the call rate, which is why it survives every other check.** A
feature that fires inside a longer one is present on nearly every genome that
has the longer one, so it looks essential. Self-recall does not see it — each
profile recovers its own collection perfectly. Collection disjointness does
not see it either: blastp at 90% found *zero* containment between these
collections and the ones they were firing inside.

`scripts/qc_duplicate_calls.py` looks for it in the scored genomes, and
separates three cases that need different answers:

  EXACT      identical coordinates. RNASE3 and SUPPRESSOR turned out to be
             99% identical over 100% coverage -- in Sweet potato chlorotic
             stunt virus the RNase3 IS the silencing suppressor, so this was
             one gene declared twice. Merge.
  NESTED     shorter inside longer. If the collections are homologous over
             most of their length it is one family split by the mass trap --
             P26/P27/P28 were 23-24% identical over 96-99% coverage, one
             diverged protein under three masses. If they are NOT homologous,
             the short profile is promiscuous: P8A, P9, P8B, P6A and P23 are
             47-90 aa with bit_cutoffs of 25-30 and fire on anything.
  PARALOGUE  both long, collections disjoint but similar. Crinivirus HSP70
             and HSP90 share no sequences yet are 64% identical over 99%
             coverage -- genuine paralogues the profiles cannot separate.
             Document it; do not "fix" it.

**Declare the overlaps that are real.** Capillovirus CP is fused inside its
ORF1 polyprotein and its MP overlaps ORF1; Tymovirus MP overlaps the
replicase by definition. Those are architecture, not defects, and the check
takes `--expected CP:POLY,MP:POLY` so it does not cry wolf. A check that
reports known biology as a finding gets ignored, and then the real finding
gets ignored with it.

**Run it before `apply_copy_num.py`,** because its whole failure mode is
manufacturing call rate.

## The last mile

A module is not built when its PSSMs exist in a temp directory. Twice in this
project a module was measured, reported as finished, and found later to exist
only under `/tmp` — not in the kit, not in Box, not in the runtime, not
committed. Use `ship_module.py --check` before calling anything done.
