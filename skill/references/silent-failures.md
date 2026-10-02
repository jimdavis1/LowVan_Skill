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
Bromoviridae had already been corrected for exactly this (94.8% -> 74.2%);
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

**In the closterovirids this is most of what the check finds, and it is
biology.** Crinivirus HSP70/HSP90 are 64% identical over 99% coverage with
disjoint collections; Closterovirus CP/CPM are 32% over 74%, because the
minor coat protein IS a duplicated coat protein gene. Both pairs are real,
both genes present, both `copy_num` assignments correct. Raising a cutoff
until the pair separates simply loses one of them.

**A containment floor is not optional.** The first version of this check
counted ANY overlap with a longer feature as nested, and adjacent genes in a
compact viral genome overlap as a matter of course -- Closterovirus HSP70
(1,794 nt) and HSP90 (1,503 nt) share 92 bases. Without the floor it reported
an entire healthy gene block as defective: Velarivirus was wrongly called at
62% clean, and two real Crinivirus features were dropped before the error was
caught. Require the shorter feature to be **80% contained** in the longer one.


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

---

### Derived numbers go stale when their source changes, and look authoritative

`min_len`/`max_len` in a module JSON are computed from the collection. Nothing
says so in the file, and they are indistinguishable from curated values. After
the homology rescue added sequences to six modules, 13 features carried bounds
computed from the *previous* collection, and `viral_genome_quality.pl` then
flagged "Feature is too short" on proteins the module had just learned to
model. Velarivirus `RDRP` reached 455 aa against a `min_len` of 466.

The derivation was confirmed before it was trusted: across 36 features, every
bound that had *not* been touched matched its collection's min and max exactly.

**After any step that changes a collection, refresh everything derived from it
in the same pass.** Rebuilding the profiles is the obvious half; the JSON is
the half that is easy to miss, because it lives in a different file and looks
hand-written.

---

### A control is not optional when the baseline predates an unrelated change

The six rescued modules were first scored against `measurements/quality.tsv`,
the only surviving per-genome baseline. It predates the 30 September `copy_num`
and vocabulary edits. On Closterovirus those edits alone retire 184 "missing
essential feature" flags, because six accessories stopped being essential — so
the rescue appeared to take the module from **0% to 36.5%** when its real
contribution was **27.0% to 36.5%**. Ampelovirus looked like 0% to 50.5% and is
actually 52.5% to **50.5%**, a net loss. Two of the three headline improvements
were somebody else's work and one was a regression wearing its clothes.

The fix is cheap and mechanical: score the *shipped* module on the *same*
panels with the *current* JSON, and compare against that.

**Before quoting a delta, date the baseline and ask what else changed in
between.** A stored table carries no record of the configuration that produced
it.

---

### The build log is a record of an intention; the collection is the fact

Bromoviridae `REP2.fasta` held 587 sequences where its build log said 590. The
three missing sequences were restored on the strength of the log. They had been
removed deliberately: blastp puts all three at 97-99.6% identity to `REP1`
across all 993 residues, so they are `1a` proteins that BV-BRC labels
`RNA-dependent RNA polymerase`. Restoring them raised `REP2`'s `max_len` from
886 to 993, produced a profile that called the 1a ORF a second time, and cost
8 genomes to *Too many contigs for RNA2*.

The signal that would have caught it immediately: **`REP2` held no sequence
above 900 aa before the change and four after.** A collection whose length
extremes move after a rescue has had its membership changed, not merely
extended. Count, min and max — not count alone.

---

### The answer was already published, on the page describing the module

Bromoviridae `REP2.fasta` held 587 sequences against a build log saying 590,
and the three missing ones were restored on the strength of the log. Section 08
of that module's own coverage audit — written a day earlier, published, and
sitting in `reports/` — says exactly what they are:

> Three 1a sequences had been binned into the REP2 collection. They formed
> their own cluster, produced their own profile, and that profile fired on 1a
> […] Removing the three sequences, rebuilding REP2 and tightening `max_len`
> to 900 took the module from **43.9% to 89.2% genuinely clean**.

So a documented fix worth 45 points was undone, cost 8 genomes, and was then
re-derived from scratch by blastp. The measurement agreed with the page, which
is some comfort, but the page was free.

It also cost a second mistake. `max_len` 900 sits deliberately *above* the
886-aa collection maximum — a margin chosen to keep REP1 and REP2 separable —
and the bounds-refresh script, which derives bounds from the collection,
narrowed it to 886. **A bound tighter than its collection is a decision; a
bound looser than its collection is stale.** The script now only ever widens.

**Before changing a module, read its coverage audit.** The audits exist to
record why a module is the way it is, and a discrepancy between two files is
exactly the question they were written to answer.

---

### A coverage floor on the query does not screen a short protein

`rescue_unassigned.py` adopted at ">= 80% identity over >= 60% of the query",
and the second clause reads like a guard against fragments. It is not one.
`ln / qlen` is coverage of the **query**, and a truncated protein aligns over
nearly all of *itself*, so it scores a high query coverage by construction. The
floor rejects a partial *alignment*; it says nothing about a partial *protein*.

What bounds truncation is the feature's length window, applied after adoption.
Measured over the six modules rescued at 0.60: of 2,196 sequences adopted into
collections, **20 (0.9%) fell below their collection's previous minimum
length** — worst Bromoviridae `REP1` at 727 aa against a prior floor of 908,
and Closterovirus `POLY` at 2137 against 2421. Those 20 are exactly what
invalidated 13 features' `min_len`/`max_len`.

The floor is now 0.85 on the curator's instruction, which tightens adoption
because most borderline hits really are partial alignments. **If short proteins
must be excluded rather than merely bounded, the instrument is a floor at the
collection's existing minimum.** Say which of the two a threshold is doing.

---

### A median is the wrong yardstick for a bimodal collection

Asked to size how many truncated proteins the rescue had admitted, the first
measurement compared each adopted sequence to its feature's median length and
reported **116 of 2,196 at least 15% short**, with Quinvirinae `CP` the worst
at 33 of 38 adoptions.

That figure was an artefact. Quinvirinae `CP` is bimodal — *Foveavirus* coat
protein is ~259 aa, the other genera ~393+ — and the collection **already held
190 sequences at exactly 259 aa**, 435 of 987 below 275. The 27 adoptions at
259 aa joined an existing dominant mode. They were not truncations, and the
feature's `LEN_OVERRIDE` of (200, 460) exists precisely because its length
varies by genus.

A tight mode at one length with a matching mass name (`28 kDa coat protein`)
is a real protein form; fragments scatter. **Before calling a sequence short,
look at the distribution it is short against** — and prefer the collection's
extremes to its median, which is what the REP2 episode had already shown.

---

### A tracking file was mistaken for a resolution

`SKILL.md` described step 2's outputs as "`UNASSIGNED_TRACKING.tsv` and
`NOT_MODELLED.tsv` **so nothing vanishes silently**." That sentence is why
three closterovirid modules shipped with their accessory proteins unmodelled.

The file does not stop anything vanishing. Everything listed in it *has*
vanished from the module; the file is the receipt. Read as a safeguard, it
turns the unassigned pool into something already handled, and the mandatory
accessory-discovery pass — cluster the pool by identity, declare each group
that reaches the member floor, leave out only true singletons — looks optional.

It also had the failure shape this document already warns about twice: **an
unclustered pool leaves exactly the evidence an empty one does.** No
`ACCESSORY_CLUSTERS.tsv`, no `UNC*` feature, nothing in any log. Velarivirus,
Ampelovirus and Closterovirus were each measured, audited, published and
declared finished in that state.

Measured once a gate existed to ask, across the 16 modules with a tracking
file: **8,757 protein occurrences in 398 distinct strings sit above the member
floor, never clustered.**

| module | strings ≥5 | occurrences | largest |
|---|---|---|---|
| Flaviviridae | 91 | 2,637 | 454 × `alternative reading frame 4 polyprotein` |
| Bromoviridae | 37 | 1,173 | 306 × `replicase` |
| Ampelovirus | 53 | 1,128 | 173 × `5 kDa protein` |
| Quinvirinae | 17 | 680 | 361 × `hypothetical protein` |
| Trivirinae | 7 | 493 | 433 × `hypothetical protein` |
| Hepeviridae | 9 | 491 | 268 × `hypothetical protein` |
| Velarivirus | 34 | 467 | 43 × `putative transmembrane protein` |
| Crinivirus | 33 | 434 | 65 × `p5.2` |
| Closterovirus | 39 | 397 | 34 × `helicase` |

Not all of those are missing features — Bromoviridae's `replicase` is the
ambiguous string its build notes deliberately leave unbinned, and some will
cluster into nothing. **That is the point: the pass that decides was never
run**, so nobody knows which.

Fixed in three places, because documentation alone is what failed:

- `SKILL.md` now says the file is a list of unfinished work and states the test.
- `collection_engine.py` prints a block naming the strings and the command,
  where the file is written.
- `install_module.py` reports the pool at install time — the last gate before
  a module ships.

**A file that records a loss is not a file that prevents one.** When a step
writes a manifest of what it could not do, something downstream has to read it.

---

### A rebuild that shrinks a profile set validates exactly like one that grows it

`build_features.py` starts a rebuild with `shutil.rmtree(dest)`. That removes
the feature's **leftover** profiles as well as its clustered ones, and the
leftover pass is a separate script that nobody re-runs unless they remember.

Rebuilding three merged Carlavirus features took the module from **338 profiles
to 234** — 104 gone, a third of the module — and `install_module.py` validated
it without complaint, because validation asks whether every declared feature
has *at least one* PSSM, not whether it still has as many as it had.
Re-running the leftover pass recovered it to 346.

This is the same shape as the stale-alignment bug the installer was already
fixed for: the check was about presence, and the defect was about quantity.

`build_features.py` now prints how many profiles a rebuild is about to discard,
so the shrink appears in the build log at the moment it happens. **When a step
deletes before it writes, log what it deleted.**

---

### "No PSSM expected" is a statement about the declaration, not about the disk

`install_module.py` printed `1 special (no PSSM expected): RDRP` for
Velarivirus, Ampelovirus and Closterovirus while shipping 10, 18 and 9 RDRP
profiles left over from before those features were redeclared
`special: transcript_edit`. The annotator called the RdRp twice -- once from
the stale profiles, once from the transcript-edit handler -- and the only
symptom was `Genome has too many HSPs for: RNA-dependent RNA polymerase`.
Clean genomes fell to **4.7%** (Velarivirus, from 64.1%) and **2.0%**
(Ampelovirus, from 51.5%), with nothing naming the cause.

The line was true of the JSON and false of the filesystem, and it was printed
in the same run that copied the profiles into place.

**When a validator prints an expectation, make it count the thing instead.**
The installer now counts profiles under every `special` feature and refuses.

---

### Requiring a value you are not going to use

The ORF1b junction reference is `g[stop-60:stop] + ORF1b` -- it needs ORF1a's
**stop** codon and nothing else. The extractor nonetheless computed both bounds
and bailed on either being absent:

```python
atg, stop = orf_bounds(g, best[0], None)
if not atg or not stop: fails["no ORF bounds"] += 1; continue
```

Every 5'-truncated assembly -- ORF1a's start codon lying off the front of the
contig -- was therefore discarded along with its perfectly good junction and
its complete RdRp. **32 genomes across three modules**, each of which then
reported `missing essential feature: RNA-dependent RNA polymerase`. Splitting
the two outputs so the RdRp and its reference need only the stop recovered
Velarivirus 52 -> 59 references, Ampelovirus 89 -> 95, Closterovirus 60 -> 73.
The ORF1a counts did not move, which is correct: a 5'-truncated assembly really
has no ORF1a.

The failure counter said `no ORF bounds` 7/11/11, which read as a data problem
and was a control-flow one. **Name the failure after the thing that was
missing, not after the function that returned it** -- `no ORF1a start codon`
would have been read in a minute.

---

### A reference set built from the genomes you evaluate on

The transcript-edit reference sets were extracted from the same genome panels
the quality run scores. That is self-recall wearing different clothes: the
measured miss rate is a floor, and the real one on unseen genomes is higher.
It is the reason the 32 recovered genomes matter more than their share of the
panel suggests -- they were the members the reference set could not reach even
with the answer in front of it.

---

### The build tree lives in scratch and nothing says to save it

Builds run in `/private/tmp` because Box sync turns many-small-file I/O into
load average with the job itself at single-digit CPU. That instruction lived in
the project's `env.sh` — *"Build in scratch; only the finished module goes to
Box"* — and nowhere in this skill. So the Enterovirus module was built,
evaluated, installed and committed, and the entire build tree stayed in
`/private/tmp`, one reboot from gone: `clusters/`, `alis/`, `corrected_alis/`,
`Curation_Report`, `Truncation_Report`, `FLAGGED_NTERM`, `BUILD_PARAMS`,
`Leftover_Seqs.aa`. None of it is reconstructible from the installed module,
and a curator cannot check the work without it.

Every other taxon already stored this under `<Taxon>/work/<Module>/`. Nothing
in the workflow said so, nothing checked, and the module looked finished.

**An instruction that lives only in the environment is not part of the
workflow.** `archive_workdir.py` now does the preservation with a required
`--dest`, `install_module.py` says so after any install run out of scratch, and
step 9 of SKILL.md documents it.

---

### `corrected_alis/` is not the whole provenance

Every module in this project preserved `PSSM-Alignments/<FEAT>/corrected_alis/`
and nothing else, which reads as complete and is not. The leftover pass and the
N-terminal re-clustering write their alignments to `reclustered_alis/`, so on
Enterovirus `corrected_alis/` held 634 alignments against **701 installed
profiles** — 67 shipped profiles whose alignment was not preserved, 32 of them
in 3D and 25 in VP1.

Nothing reported it. The counts only disagree if you compare them, and the two
live in different trees.

`archive_workdir.py --repo` now compares the preserved alignments against the
installed PSSMs per feature and names the profiles that have none. Note the id
shapes differ — a profile is `<Module>.<FEAT>.<n>.pssm` and its alignment is
`<n>.fa` — and comparing filenames instead of cluster ids reports *every*
profile as missing, which is how the first version of that check behaved.

---

### A parameter that is recorded as applied, and is not

`build_features.py` keeps its defaults in a dict keyed WITHOUT the leading
dash — `DEFAULTS = dict(m="5", mi="0.8", ...)` — and built the command from
`params["m"]`. A `features.json` written with dashes, which is how the flags
appear everywhere else including in this script's own `BUILD_PARAMS` output,
therefore did this:

```python
params = dict(DEFAULTS); params.update({"-m": 2})   # adds a NEW key
cmd = [... "-m", params["m"] ...]                   # still "5"
```

Nothing failed. The run proceeded at `-m 5` and `-mi 0.8` while `BUILD_PARAMS`
recorded `departures: {"-m": 2}` — the file asserted a departure that had not
happened, which is worse than recording nothing.

Avihepatovirus has 277 unique proteins over 264 genomes, so features hold 2 to
5 sequences. Asked for `-m 2` and silently given `-m 5`, **every one of its
twelve features built zero profiles** and the evidence on disk said they had
been built at `-m 2`. Two rounds of lowering `-mi` changed nothing, for the
same reason. With the parameters actually applied, eleven of twelve built at
once.

`_normalise()` now accepts either spelling and is **fatal on an unknown key**,
because the failure mode of the alternative is a parameter that looks set and
is not.

**A departure recorded but not applied is a false record.** Where a script
writes down what it did, build that record from the same values it passed.

---

### `-m 2` does not mean two

The clustering tool treats a cluster of exactly `min_seqs` as small and sends
it to the leftovers: `Processing MMSeqs clusters (min_seqs=2) ... Found 1
clusters, Created 0 cluster files, Put 2 sequences from 1 small clusters into
Leftover_Seqs.aa`. So `-m 2` requires three members, and a feature with two
unique sequences builds nothing from the main pass however low `-mi` goes.

Four features across Avihepatovirus and Senecavirus were in exactly that
state. `-m 1` built all four. The leftover pass would also have reached them,
but only if the main pass had created the feature directory first — and it had
not, because it "failed".
