#!/usr/bin/env python3
"""Make a figure file publishable on its own, without changing the figure.

gen_genome_map.py writes a <figure> FRAGMENT, meant to be lifted into the
coverage audit by a {"type":"figure","from":...} block. The fragment carries
no styles, because inside the audit it inherits them.

Published as a page of its own it inherits nothing, and that is not a visible
failure -- it is a silent one. Every gene arrow in the map is filled with
`style="fill:var(--built)"`. An undefined custom property makes the whole
declaration invalid, so the fill is dropped and the browser paints the shape
the default BLACK. The map still renders, at the right scale, with the right
labels; only the colour is gone, so essential CDS, other called CDS and mature
peptide become three identical black shapes under a legend that claims they
differ. Bornaviridae Genome Map was published that way on 18 Nov 2026.

So this prepends the audit's own <style> block -- not a copy written by hand,
which would drift from it -- plus a <title> and the page gutter. The figure
itself is passed through byte for byte, so the audit still lifts exactly what
it lifted before and the page keeps regenerating identically.

    wrap_figure_page.py --figure reports/x-genome-map.html \
                        --audit  reports/x-coverage-audit.html \
                        --title  "X Genome Map" --write

Re-running it is safe: anything before the first <figure is rebuilt, never
stacked up.
"""
import argparse, os, re, sys

#  The standalone page needs three things the audit gives the figure for free:
#  a gutter, a figure that cannot push the page sideways, and a little air.
EXTRA = """
  body { padding-block: 28px;
         padding-left: max(16px, env(safe-area-inset-left, 0px));
         padding-right: max(16px, env(safe-area-inset-right, 0px)); }
  figure { margin: 0 auto; max-width: 1000px; }
  figure svg { width: 100%; height: auto; max-width: 100%; }
"""


def figure_of(text):
    i = text.find("<figure")
    j = text.rfind("</figure>")
    if i < 0 or j < 0:
        sys.exit("no <figure> in that file")
    return text[i:j + len("</figure>")]


def styles_of(text):
    blocks = re.findall(r"<style>(.*?)</style>", text, re.S)
    if not blocks:
        sys.exit("no <style> in the audit; cannot inherit the tokens")
    return "\n".join(blocks)


def build(fig, css, title):
    return ("<title>%s</title>\n"
            '<meta name="claude:viewport-fit-sides" content="cover">\n'
            "<style>%s\n%s</style>\n%s\n" % (title, css, EXTRA, fig))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figure", required=True)
    ap.add_argument("--audit", required=True)
    ap.add_argument("--title", required=True)
    ap.add_argument("--out")
    ap.add_argument("--write", action="store_true",
                    help="rewrite --figure in place")
    a = ap.parse_args()

    raw = open(a.figure).read()
    fig = figure_of(raw)
    page = build(fig, styles_of(open(a.audit).read()), a.title)

    out = a.out or (a.figure if a.write else None)
    if not out:
        sys.stdout.write(page)
        return 0
    open(out, "w").write(page)
    #  The figure is the contract with the audit. Say that it survived.
    assert figure_of(open(out).read()) == fig, "figure changed -- refusing"
    print("%s  %d bytes  (figure unchanged, %d bytes)"
          % (out, os.path.getsize(out), len(fig)))
    return 0


sys.exit(main())
