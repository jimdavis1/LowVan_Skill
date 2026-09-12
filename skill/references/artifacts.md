# The two reports

Both are published as Artifacts. They answer the two questions anyone reviewing
a module actually asks: *what did the annotation cleanup achieve*, and *does the
module work*.

Load the `artifact-design` skill before writing either page.

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

## Both pages

- put the counts in the page, not in the chat
- name what is missing as prominently as what works
- the colours in `assets/registry.css` are a validated categorical pair
  (`#2a78d6` / `#eb6834` light, `#3987e5` / `#d95926` dark) — reuse rather than
  re-pick, and re-run the `dataviz` validator if you change them
- save the generator alongside the HTML in `Reports/generators/` so the page can
  be regenerated after the next rebuild
