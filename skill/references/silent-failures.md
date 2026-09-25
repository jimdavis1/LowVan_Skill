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

**Before setting `copy_num`, count the taxa in the feature's collection.**
Genus-specific accessories, readthrough products, and anything variable in
presence across the module get none.

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

## The last mile

A module is not built when its PSSMs exist in a temp directory. Twice in this
project a module was measured, reported as finished, and found later to exist
only under `/tmp` — not in the kit, not in Box, not in the runtime, not
committed. Use `ship_module.py --check` before calling anything done.
