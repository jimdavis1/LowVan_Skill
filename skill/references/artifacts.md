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
page. The held-out evaluation is the real test and belongs alongside it —
`run_gto_eval.py` over a panel drawn from the build's own BV-BRC contigs with
the reference genomes excluded, and the audit should say how the panel was
drawn. A panel that includes the genomes the references came from is self-recall
again, one level up.

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

Rarefy over **proteins**: pool every annotated protein (one row per called
feature, one per source product), shuffle, draw them one at a time and count
distinct strings, averaged over replicates. The x-axis is proteins sampled and
both curves start at (0, 0). A controlled vocabulary **saturates**, because the
*n*th protein reuses a name already seen; free text keeps climbing, because
every submitter spells the same protein a new way.

The kit shuffled whole genomes until 4 October 2026. That starts each curve at
one genome's vocabulary instead of 0 and puts genomes on the x-axis; every
saturation page built before then needs regenerating.

`--old` is a headerless `genome_id<TAB>product` TSV of the source strings; cut it
from the dump, one row per protein. Report the **ratio** at equal depth: the
two pools differ in size, so compare distinct strings at the protein count both
curves reach (`common_n` in the JSON). The Hepeviridae and Matonaviridae figures
once quoted here were genome-level and are withdrawn until regenerated.

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
numbers come from `run_gto_eval.py`, `check_rep_contigs.py`,
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


## Two things that block these pages after the fact

**The collapse page needs `synonyms.tsv`, and nothing warns you if it is
missing until you try to build the page — long after the module shipped.**
No `build_collections.py` emitted it before 29 September 2026, so ten
Alsuviricetes modules reached the repo with no way to produce their headline
artifact. `collection_engine.py` writes it now.

For a module already built, `scripts/synmap_from_collections.py`
reconstructs it from the dump and the collections on disk **without re-running
the binning** — which would truncate the collections and discard the rescue.
Be aware it answers a slightly different question:

    engine's file        a RULE matched this string        vocabulary
    reconstructed file   this string reached a collection  post-QC membership

The reconstructed version counts sequences dropped on length, ambiguity or a
cleavage check as unbound, so its bound fraction is lower. Say which one a
page used, and do not compare the figure across modules built the two ways.

**The saturation page needs per-genome feature tables, which only exist if
the evaluation was run with `--tbl-dir`.** They are written into `_wd` and
deleted with it otherwise. If the module has already been evaluated without
it, the page cannot be built without re-annotating every genome — at roughly
86 seconds per genome, that is hours for a taxon of any size. **Pass
`--tbl-dir` on the first evaluation run**, every time, whether or not you
think you want the page.

Do not hand-synthesise a feature table to get around this. It has been tried:
the gene-symbol column came out empty and the gap analysis reported
"0 genomes" without erroring.
