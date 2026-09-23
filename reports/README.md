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
