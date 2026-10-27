#!/usr/bin/env python3
"""Render the string-collapse artifact from collect_synmap.py's output.

Shows how many distinct BV-BRC product strings each controlled annotation string
absorbs -- the headline result the LowVan manuscript reports -- plus the strings
that bind to nothing, which is what the module does not cover.

    python3 collect_synmap.py <workdir>          # -> agg.json unbinned.json typos.json
    python3 gen_collapse.py --taxon Rhabdoviridae --out collapse.html

Reads agg.json / unbinned.json / typos.json from the current directory.
"""

import argparse
import collections
import datetime
import html
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ASSETS = os.path.join(os.path.dirname(HERE), "assets")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--taxon", required=True)
    ap.add_argument("--display", default=None,
                    help="heading text; spells out an abbreviated module name, "
                         "e.g. 'SRLV -- small ruminant lentivirus'. Defaults to --taxon.")
    ap.add_argument("--agg", default="agg.json")
    ap.add_argument("--unbinned", default="unbinned.json")
    ap.add_argument("--typos", default="typos.json")
    ap.add_argument("--standfirst", help="override the opening paragraph (HTML allowed)")
    ap.add_argument("--date")
    ap.add_argument("--out", default="collapse.html")
    args = ap.parse_args()

    agg = json.load(open(args.agg))
    unb = json.load(open(args.unbinned)) if os.path.exists(args.unbinned) else []
    typ = json.load(open(args.typos)) if os.path.exists(args.typos) else []
    date = args.date or datetime.date.today().strftime("%-d %b %Y")

    nbound = sum(r["nstrings"] for r in agg)
    nocc = sum(r["total"] for r in agg)
    nanno = len({r["anno"] for r in agg})
    #  DISTINCT means the union, not the sum. One string can be bound to two
    #  features -- Coltivirus deposits "VP6" on both VP6 and VP7, because the
    #  two labels are swapped between submitters -- and summing per-feature
    #  counts then reports 35 "distinct strings" where there are 33. The sum
    #  is still the right number for the per-feature cards below, so both are
    #  kept and only the headline uses the union.
    bound_strings = {s[0] for r in agg for s in r["srcs"]}
    ndistinct = len(bound_strings | {u[0] for u in unb})
    ratio = (ndistinct / nanno) if nanno else 0
    #  The headline tile has to use the same union as the standfirst. It used
    #  `nbound` -- the SUM of per-feature counts -- under the label "strings that
    #  bind to a feature", so a taxon with any string bound to two features
    #  showed a tile and a sentence that disagreed: HIV-1 rendered "159 strings
    #  that bind" beside "149 of them bind", three generic strings such as
    #  "undefined CDS" being counted once per feature. The sum is still correct
    #  for the per-feature cards, and the pair count is reported on stderr.
    nbound_unique = len(bound_strings)

    cards = []
    for r in agg:
        rows = "".join(
            '<tr><td class="tp">%s</td><td class="m num">%s</td><td class="m">%s</td></tr>'
            % (html.escape(s), "{:,}".format(c), html.escape(", ".join(g))[:90])
            for s, c, g in r["srcs"])
        cards.append(
            '<section class="card">\n'
            ' <div class="chd"><span class="key">%s</span><span class="arrow">&rarr;</span>'
            '<span class="target">%s</span>'
            '<span class="cmeta">%s string%s &middot; %s occurrences</span></div>\n'
            ' <div class="tw"><table><thead><tr><th>BV-BRC product string</th>'
            '<th class="num">Count</th><th>Genera</th></tr></thead><tbody>%s</tbody></table></div>\n'
            '</section>'
            % (html.escape(r["key"]), html.escape(r["anno"]),
               r["nstrings"], "" if r["nstrings"] == 1 else "s",
               "{:,}".format(r["total"]), rows))

    typrows = "".join(
        '<tr><td class="tp">%s</td><td class="m num">%s</td><td class="tp">%s</td>'
        '<td class="m num">%s</td><td class="m">%s</td></tr>'
        % (html.escape(v), "{:,}".format(c), html.escape(top), d,
           html.escape(", ".join(g))[:60])
        for _k, top, v, c, d, g in typ) or \
        '<tr><td colspan="5" class="dash">none found</td></tr>'

    unbrows = "".join(
        '<tr><td class="tp">%s</td><td class="m num">%s</td><td class="m">%s</td></tr>'
        % (html.escape(s), "{:,}".format(c), html.escape(", ".join(g))[:90])
        for s, c, g in unb[:120]) or \
        '<tr><td colspan="3" class="dash">every string binds</td></tr>'

    stand = args.standfirst or (
        "The %s dump uses <b>%d distinct product strings</b> for its proteins. "
        "<b>%d of them bind</b> to a feature and collapse into <b>%d controlled "
        "annotation strings</b> &mdash; a %.1f-fold reduction across %s protein "
        "occurrences. The remaining <b>%d bind to nothing</b>, and are listed at "
        "the end rather than hidden."
        % (html.escape(args.taxon), ndistinct, len(bound_strings), nanno, ratio,
           "{:,}".format(nocc), len(unb)))

    css = open(os.path.join(ASSETS, "collapse.css")).read()

    #  The taxon has to be in the <title>: it is the artifact's name in the
    #  gallery, and five modules all called "BV-BRC String Collapse" cannot be
    #  told apart there.
    page = """<title>%(TAXON)s String Collapse</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Serif:wght@500;600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>%(CSS)s</style>
<div class="wrap">
 <header>
  <div class="eyebrow">LowVan &middot; %(TAXON)s module &middot; %(DATE)s</div>
  <h1>%(TAXON)s</h1>
  <h2 style="margin:8px 0 0;font-weight:500;font-size:20px;color:var(--ink2)">What the annotation cleanup actually collapsed</h2>
  <p class="stand">%(STAND)s</p>
 </header>
 <div class="figs">
  <div class="fig"><span class="n">%(NDISTINCT)d</span><span class="k">distinct BV-BRC strings</span></div>
  <div class="fig"><span class="n">%(NBOUND)d</span><span class="k">strings that bind to a feature</span></div>
  <div class="fig"><span class="n">%(NANNO)d</span><span class="k">controlled strings emitted</span></div>
  <div class="fig"><span class="n">%(NOCC)s</span><span class="k">protein occurrences bound</span></div>
  <div class="fig"><span class="n">%(NUNB)d</span><span class="k">strings that never bind</span></div>
 </div>
%(CARDS)s
 <section class="gap">
  <h2>Spelling variants folded in</h2>
  <p>Source strings within two edits of the most common string for the same feature.
  These are submission typos, not different proteins &mdash; the clearest evidence of
  what the underlying vocabulary looks like.</p>
  <div class="tw"><table><thead><tr><th>Variant</th><th class="num">Count</th>
  <th>Most common form</th><th class="num">Edits</th><th>Genera</th></tr></thead>
  <tbody>%(TYPOS)s</tbody></table></div>
 </section>
 <section class="gap">
  <h2>Strings that never bind</h2>
  <p>Every product string that matched no binning rule. This is the honest measure of
  what the module does not cover: some are proteins deliberately not modelled, some are
  accessories with too few sequences, and some are rules still to be written.</p>
  <div class="tw"><table><thead><tr><th>Unbound product string</th>
  <th class="num">Count</th><th>Genera</th></tr></thead><tbody>%(UNB)s</tbody></table></div>
 </section>
 <div class="src">Generated from <code>synonyms.tsv</code> and the module JSON by
 <code>collect_synmap.py</code> and <code>gen_collapse.py</code>. Counts are protein
 occurrences (features), not unique sequences, so a protein shared by 40 genomes counts 40
 times &mdash; that is the number a database user experiences.</div>
</div>
""" % {"CSS": css, "TAXON": html.escape(args.display or args.taxon), "DATE": html.escape(date),
       "STAND": stand, "NDISTINCT": ndistinct, "NBOUND": nbound_unique, "NANNO": nanno,
       "NOCC": "{:,}".format(nocc), "NUNB": len(unb),
       "CARDS": "\n".join(cards), "TYPOS": typrows, "UNB": unbrows}

    with open(args.out, "w") as fh:
        fh.write(page)
    print("wrote %s (%d bytes)" % (args.out, len(page)))
    print("  %d distinct strings -> %d controlled (%.1fx), %d unbound, %s occurrences"
          % (ndistinct, nanno, ratio, len(unb), "{:,}".format(nocc)))
    if nbound != len(bound_strings):
        print("  (%d string-feature pairs: %d string(s) bind to more than one "
              "feature)" % (nbound, nbound - len(bound_strings)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
