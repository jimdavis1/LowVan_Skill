# Reports

The source HTML for every published page.

**The index of live URLs is [../ARTIFACTS.md](../ARTIFACTS.md)** — generated
from `published.tsv` by `gen_artifacts_index.py`, and the only place the links
are kept. Do not start a second list here; that is how they got lost the first
time.

Artifacts are published to claude.ai, not to GitHub. Opening a file in this
directory renders it locally; the shared copy is behind the URL in the index.
Keeping the source here means any page can be republished without rebuilding
its measurements.

Each page is generated from a facts file by a script in `generators/` or
`../skill/scripts/gen_*.py`. Edit the facts and regenerate — do not hand-edit
the HTML.

## Refreshing the index

```bash
#  in a Claude Code session: Artifact list --scope mine
#  paste the rows into published.tsv, then
python3 reports/gen_artifacts_index.py
```

The script matches each source file to its published page on the file's own
`<title>`, not on its filename, so renaming a file cannot silently unlink it.
It also separates *built but never published* from *never made*, and picks the
current URL for a repeated title by date rather than by listing order.

Titles are matched on the **canonical taxon**. `gen_registry.py` and its
siblings put a `--display` gloss in the `<title>`, so pages go out titled
`<taxon> — <plain-English gloss> <Kind>`. Read literally, the gloss parses as
part of the taxon: fourteen modules were split across two half-empty rows, and
forty-seven pages showed no source file although the file was sitting here, be-
cause the live title and the file's title glossed differently. The gloss is
display text, never identity, on either side of the lookup.

## Recovering a facts file

A page with no facts file can only be hand-edited, which the rule above
forbids. Before concluding one is lost, invert the generator:

```bash
python3 ../skill/scripts/facts_from_html.py <stem>-coverage-audit.html
python3 ../skill/scripts/recover_facts.py          # keeps only what round-trips
```

`gen_coverage_audit.py` emits a deterministic page, so the mapping inverts.
The test is not inspection — `recover_facts.py` regenerates from the recovered
facts and keeps the result only when the page comes back **byte-identical**.
32 were recovered this way on 17 Nov 2026; 73 of 77 audits now regenerate
exactly. Four — allexivirus, alphaflexiviridae, rhabdoviridae, togaviridae —
are a different page design (no `<section>`, bespoke `kpi`/`gloss` classes) and
were left alone rather than silently redesigned.

## Publishing a figure on its own

`gen_genome_map.py` writes a bare `<figure>` fragment, meant to be lifted into
an audit by a `{"type":"figure","from":...}` block. It carries no styles,
because inside the audit it inherits them.

Published as its own page it inherits nothing, and the failure is silent rather
than visible. Every gene arrow is filled `style="fill:var(--built)"`; an
undefined custom property invalidates the declaration, so the fill is dropped
and the browser paints the shape default black. The map still draws at the
right scale with the right labels — only the colour is gone, so *essential
CDS*, *other called CDS* and *mature peptide* become three identical black
shapes under a legend claiming they differ.

```bash
python3 ../skill/scripts/wrap_figure_page.py \
    --figure <stem>-genome-map.html \
    --audit  <stem>-coverage-audit.html \
    --title  "<Taxon> Genome Map" --write
```

It prepends the audit's **own** `<style>` block, so the page cannot drift from
the figure's source, and passes the figure through byte for byte — asserted
before writing, because the audit still lifts it. Re-running is safe. After a
wrap, confirm the audit still regenerates byte-identically.

`Bornaviridae Genome Map` (`PGcDVpmXmmqcJGcZDFN2SC`) is live in the unstyled
state and has not been republished.

## Keeping the runtime readable

Box serves `_infrastructure/Viral_Annotation` as placeholder files. When they
stop hydrating, listings and sizes return instantly while every attempt to READ
a `.pssm` times out, and `tblastn` reports `Unsupported format for PSSM` —
which reads like a corrupt profile and is an I/O failure. Run against a local
copy instead: `../skill/scripts/sync_runtime.py --git --from-box`, which merges
the (behind) upstream repo with whatever Box will serve, resumes, and judges a
module complete on file count rather than on the directory existing.
