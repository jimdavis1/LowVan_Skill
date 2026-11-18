#!/usr/bin/env python3
"""Regenerate ARTIFACTS.md from reports/published.tsv.

The URLs were getting lost. They lived in a RESUME file, in prose, and in one
half-complete table, and three of the forty-five appeared nowhere at all. This
makes one index, generated rather than maintained by hand, so it cannot drift
from the thing it describes.

Two distinctions the hand-written version got wrong and this one enforces:

  * A page with a SOURCE FILE but no live URL was built and never published.
    That is not the same as a page that was never made, and collapsing them
    hides finished work sitting on disk.
  * When one title has been published more than once, the CURRENT one is the
    one with the latest date -- not the one that happened to come first in the
    listing. The listing is only roughly ordered, and taking its order made
    Rhabdoviridae's 3 September registry current over its 9 September one.

Source files are matched on each file's own <title>, never on its filename, so
renaming a file cannot silently unlink it.

    python3 reports/gen_artifacts_index.py
"""
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

#  Taxa are DERIVED from the published titles, not listed here. The list used
#  to be hardcoded, and every module added after it was written -- ten of them
#  -- had its four pages silently filed under "Not tied to one taxon", which
#  is where the index stopped being usable for finding anything.
TAXA_EXTRA = ["Rhabdoviridae"]        # published before the naming convention settled
#  Titles that parse as "<taxon> <kind>" but whose leading token is not a
#  taxon. "BV-BRC String Collapse" is a cross-taxon page; "Rubella Vocabulary
#  Saturation" is Matonaviridae's own page published under the disease name,
#  and ALIAS already files it there -- listing "Rubella" as well would show
#  one page twice under two headings.
#  "Mammalian Orthoreovirus" is the same taxon as "MammalianOrthoreovirus":
#  three of its four pages were published without the space and the coverage
#  audit with it. ALIAS below files the audit under the closed-up spelling, so
#  listing the spaced one as well would split one module across two headings.
NOT_TAXA = {"BV-BRC", "Rubella", "Mammalian Orthoreovirus",
            #  one page went out with a space where the module, and its
            #  other three pages, use an underscore; ALIAS files it under
            #  the underscore, so listing the spaced form too would split
            #  one module across two rows.
            "Mastadenovirus MultiFiber"}
#  "Genome Map" was added 18 Nov 2026. Every module is supposed to ship a
#  genome map (SKILL.md step 12), but until Bornaviridae they were only ever
#  EMBEDDED in the coverage audit, never published as pages of their own, so
#  the column is empty for every earlier taxon. That is accurate rather than
#  missing data: those maps exist, inside their audits.
KINDS = [("Coverage Audit", "coverage audit"),
         ("PSSM Registry", "PSSM registry"),
         ("String Collapse", "string collapse"),
         ("Vocabulary Saturation", "vocabulary saturation"),
         ("Genome Map", "genome map")]
#  Pages whose title does not follow "<taxon> <kind>".
ALIAS = {("Matonaviridae", "Vocabulary Saturation"): "Rubella Vocabulary Saturation",
         ("MammalianOrthoreovirus", "Coverage Audit"):
             "Mammalian Orthoreovirus Coverage Audit",
         #  same module, two spellings: three pages went out with the
         #  underscore the module directory uses and the audit with a space.
         ("Mastadenovirus_MultiFiber", "Coverage Audit"):
             "Mastadenovirus MultiFiber Coverage Audit"}

#  Many pages are titled "<taxon> — <plain-English gloss> <Kind>", because
#  gen_registry.py and friends take a --display string and put it in the
#  <title>. Left alone that splits one module across two matrix rows, each
#  holding a different half of its pages: Picobirnaviridae showed a coverage
#  audit and three dashes while its collapse, registry and saturation sat
#  under "Picobirnaviridae -- picobirnaviruses". Fourteen modules were split
#  this way on 18 Nov 2026. The gloss is display text, not identity, so the
#  taxon is canonicalised to what precedes it and the full title is kept for
#  the lookup.
GLOSS = re.compile(r"\s+(?:—|--)\s+")


def canonical(taxon):
    return GLOSS.split(taxon, 1)[0].strip()


def load():
    rows = []
    for line in open(os.path.join(HERE, "published.tsv")):
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) == 3:
            rows.append((p[0].strip(), p[1].strip(), p[2].strip()))
    return rows


def sources():
    by = {}
    for f in sorted(os.listdir(HERE)):
        if not f.endswith(".html"):
            continue
        m = re.search(r"<title>([^<]*)</title>",
                      open(os.path.join(HERE, f), errors="replace").read())
        if m and m.group(1).strip():
            by.setdefault(m.group(1).strip(), f)
    return by


def slug(taxon):
    return re.sub(r"[^a-z0-9]+", "", taxon.lower())


def map_source(taxon):
    """Genome maps are written by gen_genome_map.py as <figure> FRAGMENTS for
    embedding in the coverage audit, so they carry no <title> and sources()
    cannot see them -- 25 of the 26 in this repo were invisible to the index.
    Match them on the filename the generator uses instead, so the map column
    says "src only" (the file is here, it was never published as its own page)
    rather than "—" (never made)."""
    s = slug(taxon)
    for f in os.listdir(HERE):
        if f.endswith("-genome-map.html") and slug(f[: -len("-genome-map.html")]) == s:
            return f
    return None


def derive_taxa(titles):
    """Every leading token set that precedes one of the kind names, with any
    plain-English gloss stripped so one module is one row."""
    found = set()
    for t in titles:
        for kind, _label in KINDS:
            if t.endswith(" " + kind):
                found.add(canonical(t[: -len(kind) - 1]))
    for tax, kind in ALIAS:
        found.add(tax)
    return sorted((found | set(TAXA_EXTRA)) - NOT_TAXA, key=str.lower)


def title_index(titles):
    """(canonical taxon, kind) -> the title actually published, gloss and all."""
    idx = {}
    for t in titles:
        for kind, _label in KINDS:
            if t.endswith(" " + kind):
                key = (canonical(t[: -len(kind) - 1]), kind)
                #  a glossed title is the fuller one; prefer it only if the
                #  plain form was never published
                if key not in idx or len(t) > len(idx[key]):
                    idx[key] = t
    return idx


def main():
    rows, src = load(), sources()

    #  Latest date wins. Ties keep the first seen, which only happens when the
    #  same page was republished the same day and either URL is the same page.
    best = {}
    for t, u, d in rows:
        if t not in best or d > best[t][1]:
            best[t] = (u, d)
    superseded = [(t, u, d) for t, u, d in rows if best[t][0] != u]

    TIDX = title_index(set(best) | set(src))

    def title_for(tax, kind):
        if (tax, kind) in ALIAS:
            return ALIAS[(tax, kind)]
        plain = "%s %s" % (tax, kind)
        if plain in best or plain in src:
            return plain
        return TIDX.get((tax, kind), plain)

    TAXA = derive_taxa(set(best) | set(src))

    def row(t):
        u, d = best[t]
        s = "[source](reports/%s)" % src[t] if t in src else "—"
        return "| %s | [open](https://claude.ai/artifact/%s) | %s | %s |" % (t, u, s, d)

    o = ["# Published artifacts", "",
         "Every page this project has published, with its live URL and its source",
         "in this repo. **Artifacts live on claude.ai, not on GitHub** — a `.html`",
         "file here renders locally in a browser; the `open` link is the shared",
         "copy. Pages are private until shared from the page's own share menu, so",
         "a link will not open for someone who has not been given access.", "",
         "Generated by `reports/gen_artifacts_index.py` from",
         "`reports/published.tsv`. Do not edit this file — refresh the TSV from",
         "`Artifact list` in a Claude Code session and re-run the script.", "",
         "Taxa are listed alphabetically. The matrix is the fast way in; the",
         "per-taxon tables below carry the dates and the source paths.", "",
         "## At a glance", "",
         "| taxon | collapse | registry | saturation | audit | map |",
         "|---|---|---|---|---|---|"]

    def cell(tax, kind):
        t = title_for(tax, kind)
        if t in best:
            return "[live](https://claude.ai/artifact/%s)" % best[t][0]
        if t in src:
            return "[src only](reports/%s)" % src[t]
        if kind == "Genome Map":
            f = map_source(tax)
            if f:
                return "[src only](reports/%s)" % f
        return "—"

    TAXA_M = derive_taxa(set(best) | set(src))
    for tax in TAXA_M:
        cells = [cell(tax, k) for k in
                 ("String Collapse", "PSSM Registry",
                  "Vocabulary Saturation", "Coverage Audit", "Genome Map")]
        if all(c == "—" for c in cells):
            continue
        o.append("| **%s** | %s |" % (tax, " | ".join(cells)))
    o += ["",
          "`live` is published and current. `src only` was built and never",
          "published — the file is in this repo. `—` was not made.", "",
          "One row is one module. A page titled `<taxon> — <gloss> <Kind>` is",
          "filed under the taxon, not under the gloss; before 18 Nov 2026 the",
          "gloss made a second row, and fourteen modules were split across two",
          "rows each holding half their pages.", "",
          "Most `map` cells read `src only` because genome maps are written as",
          "`<figure>` fragments to embed in the coverage audit, so they sit in",
          "this repo without ever being published as pages of their own.", "",
          "## Per taxon", ""]

    used = set()
    for tax in TAXA:
        live, built, never = [], [], []
        for kind, label in KINDS:
            t = title_for(tax, kind)
            used.add(t)
            if t in best:
                live.append(t)
            elif t in src:
                built.append((label, src[t]))
            else:
                never.append(label)
        if not (live or built):
            continue
        o += ["### %s" % tax, "", "| page | live | source | updated |", "|---|---|---|---|"]
        o += [row(t) for t in live]
        o.append("")
        if built:
            o.append("Built but never published, source is here: " +
                     ", ".join("[%s](reports/%s)" % (l, f) for l, f in built) + ".")
        if never:
            o.append("Not made: " + ", ".join(never) + ".")
        if built or never:
            o.append("")

    rest = sorted((t for t in best if t not in used), key=lambda t: best[t][1], reverse=True)
    o += ["## Not tied to one taxon", "", "| page | live | source | updated |",
          "|---|---|---|---|"] + [row(t) for t in rest] + [""]

    o += ["## Superseded", "",
          "The same title published more than once. The tables above carry the",
          "current URL; these earlier ones still resolve and should not be",
          "circulated.", "", "| page | live | updated |", "|---|---|---|"]
    o += ["| %s | [open](https://claude.ai/artifact/%s) | %s |" % (t, u, d)
          for t, u, d in sorted(superseded, key=lambda r: (r[0], r[2]))]
    o.append("")

    open(os.path.join(ROOT, "ARTIFACTS.md"), "w").write("\n".join(o) + "\n")
    print("ARTIFACTS.md: %d published, %d current, %d superseded"
          % (len(rows), len(best), len(superseded)))
    nosrc = sorted(t for t in best if t not in src)
    if nosrc:
        print("no source file in reports/ for %d:" % len(nosrc))
        for t in nosrc:
            print("   -", t)
    return 0


sys.exit(main())
