# The four reports

All four are published as Artifacts. They answer the questions anyone reviewing
a module actually asks: *what did the annotation cleanup achieve*, *does the
module work*, *how fast does the vocabulary saturate*, and *what can it not do*.

The **coverage audit is the one that gets skipped**, because unlike the others it
has no single input file — it reads results from steps 10, 11 and 12 together.
Build it last, and build it every time.

Load the `artifact-design` skill before writing any of these pages.

## 1. String collapse — what the cleanup achieved

How many distinct BV-BRC product strings each annotation string absorbs. This is
the number the LowVan manuscript reports as its headline result.

```bash
python3 scripts/collect_synmap.py <workdir>     # -> agg.json, unbinned.json, typos.json
```

Reads `synonyms.tsv` (from the triage step) joined to the module JSON. Produces:

- `agg.json` — per feature key: the annotation string, total feature count, and
  every source string that collapsed into it with its own count and genera
- `unbinned.json` — every string that matched no rule, with counts and genera
- `typos.json` — source strings within two edits of the most common string for
  the same key, which are usually submission typos

Render with the `assets/collapse.css` stylesheet.

**What to say on the page.** The collapse ratio, then the interesting parts:

- the biggest absorbers — one annotation swallowing dozens of spellings is the
  whole argument for the project
- the typo cluster — concrete evidence of what the source data looks like
- the unbinned list, honestly. This is what the module does *not* cover, and
  hiding it makes the page untrustworthy.

## 2. PSSM registry — does the module work

Every declared feature, the profiles built for it, and what fraction of its own
collection those profiles recover.

```bash
python3 scripts/collect_registry.py --workdir . --out registry.json
python3 scripts/gen_registry.py --registry registry.json --taxon <Family> \
    --notes notes.json --date "3 Sep 2026" --out registry.html
```

`collect_registry.py` does the measurement: for each feature it runs that
feature's own PSSMs against that feature's own collection with
`psiblast -in_pssm` and counts sequences scoring at or above the JSON
`bit_cutoff`.

**Be honest about what self-recall means.** It is a floor, not a validation. A
profile that cannot recover the sequences it was built from will certainly not
recover anything new — but recovering them proves very little. Say so on the
page. The held-out evaluation from `evaluate_module.py` is the real test, and
belongs alongside it.

Statuses assigned per feature:

| Status | Meaning |
|---|---|
| `built` | has PSSMs and a collection |
| `derived` | PSSMs but no collection of its own (mat_peptides) |
| `special` | called by an external program (transcript_edit / splice) |
| `too-few` | collection smaller than a profile can support |
| `no-collection` | declared with no sequences behind it at all |

The last two are the gaps, and the page must separate them because they are
different problems. `no-collection` means the protein is not in BV-BRC under any
string the rules match. `too-few` means it is there but too thin — often because
the key is a positional label covering non-homologous proteins, in which case
the right answer is not to build it.

### The notes file

Auto-generated counts explain nothing. Write the prose:

```json
{
  "standfirst": "one paragraph under the headline, HTML allowed",
  "modules": {"<Module>": "why this module looks the way it does"},
  "gaps":    "prose for the 'what was missed' section",
  "derived": "prose about mat_peptide / derived features",
  "why":     {"<FEATURE>": "why this one could not be built"}
}
```

The `why` entries matter most. `n=1` on its own reads as laziness; *"1 sequence;
positional label, not a homology group"* reads as a decision. Every `too-few`
feature deserves one line saying whether it was checked and what was found.

## 3. Vocabulary saturation — the same claim as a curve

```bash
python3 scripts/annotation_rarefaction.py --new coverage_eval/ann \
        --old bvbrc_products.tsv --out rarefaction.json --replicates 100
python3 scripts/gen_rarefaction.py --rarefaction rarefaction.json \
        --taxon <Family> --out rarefaction.html
```

Shuffle the genomes, accumulate distinct annotation strings, average over
replicates. A controlled vocabulary **saturates**, because the *n*th genome
reuses names the first *n*−1 established; free text climbs roughly linearly,
because every submitter spells the same protein a new way.

`--old` is a headerless `genome_id<TAB>product` TSV of the source strings; cut it
from the dump. Report the **ratio**, not just the picture: strings per 100
genomes at the end of each curve says how much collapse the vocabulary bought.

Hepeviridae: free text reaches 95% of its final vocabulary at 533 genomes and is
still climbing, the controlled set at **13**; 10.2 strings per 100 genomes
against 0.6, a 15.8-fold difference. Matonaviridae, a conserved taxon whose
submitters were more consistent, gives 13.9 against 3.2 — only 4.3-fold. Both
numbers are honest; the shape of the curve is the claim, not the ratio alone.

## 4. Coverage audit — the whole module in nine questions

```bash
python3 scripts/gen_coverage_audit.py --facts <taxon>_audit.json \
        --out <taxon>-coverage-audit.html
```

Nine sections, in this order, because it is the order a reviewer asks them in:

| § | question |
|---|---|
| 01 | what the module covers |
| 02 | how it was built |
| 03 | where the cutoffs sit |
| 04 | does it route? |
| 05 | does it call the proteins? |
| 06 | are the calls right? |
| 07 | **what it cannot do** |
| 08 | what the run found |
| 09 | what the vocabulary bought |

The generator takes a **facts file you write by hand**. That is deliberate: the
numbers come from `evaluate_module.py`, `check_rep_contigs.py`,
`evaluate_coverage.py`, `run_gto_eval.py`, `collect_synmap.py` and
`annotation_rarefaction.py`, and a page that scrapes a working directory for
whichever of those happen to have run is worse than one whose inputs are
explicit. Every figure in the output appears in the JSON, so the page can be
checked against `BUILD_NOTES.md` line by line.

Block types available: `tiles`, `bars`, `table`, `pre`, `prose`, `verdict`
(`kind` is `ok`, `warn` or `stop`). Copy a shipped facts file and replace the
numbers rather than starting from the schema.

**Write §07 from the measurements.** "What it cannot do" is the section that
makes the rest credible. Name the gap, give the number behind it, and say why
tuning will not close it — a feature with two training sequences is not a
threshold problem, and saying so is more useful than a percentage. If §07 is
thin, the audit has not been done.

## All four pages

- put the counts in the page, not in the chat
- name what is missing as prominently as what works
- the colours in `assets/registry.css` are a validated categorical pair
  (`#2a78d6` / `#eb6834` light, `#3987e5` / `#d95926` dark) — reuse rather than
  re-pick, and re-run the `dataviz` validator if you change them
- save the generator, and the coverage audit's facts file, alongside the HTML in
  `Reports/generators/` so every page can be regenerated after the next rebuild
