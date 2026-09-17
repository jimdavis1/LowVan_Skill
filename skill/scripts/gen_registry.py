#!/usr/bin/env python3
"""Render the PSSM registry artifact from registry.json.

Every number on the page is computed from registry.json (produced by
collect_registry.py). The taxon-specific commentary is not: pass it in with
--notes so the page explains *this* module rather than reciting counts. A
registry without commentary is still correct, just less useful.

    python3 collect_registry.py --workdir . --out registry.json
    python3 gen_registry.py --registry registry.json --taxon Rhabdoviridae \
        --notes notes.json --out registry.html

notes.json is optional and entirely free-form except for its keys:

    {
      "standfirst": "one paragraph under the headline, HTML allowed",
      "modules": {"Alpharhabdovirinae": "why this module looks the way it does"},
      "gaps":    "prose for the 'what was missed' section",
      "derived": "prose about mat_peptide/derived features",
      "why":     {"SRIPU_MX": "single sequence; positional label, not a homolog"}
    }
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

SLABEL = {"built": "built", "derived": "derived", "special": "external",
          "too-few": "too few", "no-collection": "no collection"}
# which statuses count as a gap the reader should be told about
GAPS = ("too-few", "no-collection")


def row_html(r):
    st = r["status"]
    pss = " ".join("<code>%s</code>" % html.escape(p.rsplit(".", 2)[-2])
                   for p in r["pssms"]) or '<span class="dash">&mdash;</span>'
    n, called = r["n"], r["called"]
    if st in ("built", "derived") and called is not None and n:
        pc = 100.0 * called / n
        cov = ('<div class="cov" title="%d of %d callable"><i style="width:%.2f%%"></i>'
               '<b style="width:%.2f%%"></b></div><span class="covn">%d%%</span>'
               % (called, n, pc, 100 - pc, round(pc)))
    else:
        cov = '<span class="dash">&mdash;</span>'
    bounds = ("%s&ndash;%s" % (r["minl"], r["maxl"])) if r["minl"] is not None else "&mdash;"
    return ('<tr data-status="%s"><td class="k">%s%s</td><td class="anno">%s</td>'
            '<td class="m">%s</td><td class="m ty">%s</td><td class="m num">%s</td>'
            '<td class="m num">%s</td><td class="m num">%s</td><td class="pss">%s</td>'
            '<td class="covcell">%s</td><td><span class="chip s-%s">%s</span></td></tr>'
            % (st, html.escape(r["key"]),
               ' <span class="nd" title="non-default build parameters">&#9702;</span>'
               if r["nondefault"] else "",
               html.escape(r["anno"]),
               html.escape(r["gene"]) or "&mdash;",
               "mat_pep" if r["ftype"] == "mat_peptide" else html.escape(r["ftype"] or "CDS"),
               r["bit"] if r["bit"] is not None else "&mdash;",
               bounds,
               n if n else '<span class="dash">&mdash;</span>',
               pss, cov, st, SLABEL.get(st, st)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--registry", default="registry.json")
    ap.add_argument("--taxon", required=True, help="family or group name for the title")
    ap.add_argument("--json-name", help="module JSON filename to cite (default <taxon>_Viral_PSSM.json)")
    ap.add_argument("--notes", help="JSON of taxon-specific prose")
    ap.add_argument("--date", help="date shown in the eyebrow (default: today)")
    ap.add_argument("--out", default="registry.html")
    args = ap.parse_args()

    with open(args.registry) as fh:
        d = json.load(fh, object_pairs_hook=collections.OrderedDict)
    notes = {}
    if args.notes:
        with open(args.notes) as fh:
            notes = json.load(fh)
    jname = args.json_name or ("%s_Viral_PSSM.json" % args.taxon)
    date = args.date or datetime.date.today().strftime("%-d %b %Y")

    allf = [(m, r) for m, v in d.items() for r in v["features"]]
    counts = collections.Counter(r["status"] for _m, r in allf)
    nall = len(allf)
    withp = [r for _m, r in allf if r["pssms"]]
    npssm = sum(len(r["pssms"]) for _m, r in allf)
    scored = [r for _m, r in allf if r["called"] is not None and r["n"]]
    nseq = sum(r["n"] for r in scored)
    ncall = sum(r["called"] for r in scored)
    pct = round(100.0 * ncall / nseq) if nseq else 0
    ngap = sum(counts[g] for g in GAPS)
    nnd = sum(1 for _m, r in allf if r["nondefault"])
    annos = len({r["anno"] for _m, r in allf})
    reach = len({r["anno"] for _m, r in allf if r["pssms"]})

    # modules ordered by how much of the work they carry
    order = sorted(d, key=lambda m: -sum(len(r["pssms"]) for r in d[m]["features"]))

    sections = []
    for m in order:
        rows = d[m]["features"]
        b = [r for r in rows if r["status"] == "built"]
        s = sum(r["n"] for r in rows if r["called"] is not None and r["n"])
        k = sum(r["called"] for r in rows if r["called"] is not None and r["n"])
        note = notes.get("modules", {}).get(m, "")
        sections.append(
            '<section class="tax">\n <div class="taxhd">'
            '<h2>%s</h2><span class="taxmeta">%d of %d features built &middot; %d PSSMs%s</span></div>\n'
            '%s'
            ' <div class="tw"><table><thead><tr><th>Key</th><th>Annotation string</th><th>Gene</th>'
            '<th>Type</th><th class="num">Bit</th><th class="num">Length</th><th class="num">Seqs</th>'
            '<th>PSSMs</th><th>Callable</th><th>Status</th></tr></thead>\n <tbody>%s</tbody></table></div>\n'
            '</section>'
            % (html.escape(m), len(b), len(rows),
               sum(len(r["pssms"]) for r in rows),
               (" &middot; %d%% callable" % round(100.0 * k / s)) if s else "",
               ('<p class="taxnote">%s</p>\n' % note) if note else "",
               "\n  ".join(row_html(r) for r in rows)))

    # gap tables
    nc = collections.defaultdict(list)
    tf = []
    for m, r in allf:
        if r["status"] == "no-collection":
            nc[m].append(r["key"])
        elif r["status"] == "too-few":
            tf.append((m, r))
    gapnc = "".join(
        '<tr><td>%s</td><td class="m num">%d</td><td class="keys">%s</td></tr>'
        % (html.escape(m), len(v), ", ".join("<code>%s</code>" % html.escape(x) for x in sorted(v)))
        for m, v in sorted(nc.items(), key=lambda x: -len(x[1])))
    gaptf = "".join(
        '<tr><td>%s</td><td class="k">%s</td><td class="m num">%d</td><td>%s</td></tr>'
        % (html.escape(m), html.escape(r["key"]), r["n"],
           html.escape(notes.get("why", {}).get(r["key"], "")) or
           "<span class='dash'>needs a look</span>")
        for m, r in sorted(tf, key=lambda x: (x[0], x[1]["key"])))

    stand = notes.get("standfirst") or (
        "Every feature declared in <code>%s</code>, the annotation string it emits, and the "
        "profiles built for it. <b>%d of %d declared features have a profile</b>, giving "
        "<b>%d PSSMs</b> over %s sequences &mdash; of which <b>%d%% are callable</b> by their own "
        "module. The remaining %d cannot be satisfied from this BV-BRC dump."
        % (html.escape(jname), len(withp), nall, npssm, "{:,}".format(nseq), pct, ngap))

    css = open(os.path.join(ASSETS, "registry.css")).read()
    js = open(os.path.join(ASSETS, "registry.js")).read()

    page = """<title>%(TAXON)s PSSM Registry</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Serif:wght@500;600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>%(CSS)s</style>
<div class="wrap">
 <header>
  <div class="eyebrow">LowVan &middot; %(TAXON)s module &middot; %(DATE)s</div>
  <h1>Annotation strings and the PSSMs that call them</h1>
  <p class="stand">%(STAND)s</p>
 </header>
 <div class="figs">
  <div class="fig"><span class="n">%(NPSSM)d</span><span class="k">PSSMs built</span></div>
  <div class="fig"><span class="n">%(NWITHP)d/%(NALL)d</span><span class="k">features with a profile</span></div>
  <div class="fig"><span class="n">%(REACH)d</span><span class="k">of %(ANNOS)d annotation strings reachable</span></div>
  <div class="fig"><span class="n">%(PCT)d%%</span><span class="k">of sequences callable</span></div>
  <div class="fig"><span class="n">%(NGAP)d</span><span class="k">features unsatisfiable</span></div>
 </div>
 <div class="legend">
  <span><i class="sw" style="background:var(--built)"></i> callable by its own module</span>
  <span><i class="sw" style="background:var(--missed)"></i> in the collection but not called</span>
  <span><i class="nd">&#9702;</i> non-default build parameters (%(NND)d)</span>
 </div>
 <div class="ctrl"><span class="eyebrow" style="margin-right:4px">show</span>
  <button aria-pressed="true" data-f="all">all %(NALL)d</button>
  <button aria-pressed="false" data-f="built">built (%(NBUILT)d)</button>
  <button aria-pressed="false" data-f="gap">gaps (%(NGAP)d)</button>
  <button aria-pressed="false" data-f="derived">derived (%(NDERIV)d)</button>
 </div>
%(SECTIONS)s
 <section class="gap">
  <h2>What was missed, and why</h2>
  <p>%(GAPPROSE)s</p>
  <div class="tw"><table><thead><tr><th>Module</th><th class="num">Features</th><th>Keys with no collection</th></tr></thead>
  <tbody>%(GAPNC)s</tbody></table></div>
  <p style="margin-top:14px"><b>%(NTOOFEW)d features have a collection too small to profile.</b></p>
  <div class="tw"><table><thead><tr><th>Module</th><th>Feature</th><th class="num">Seqs</th><th>Why not buildable</th></tr></thead>
  <tbody>%(GAPTF)s</tbody></table></div>
  %(DERIVED)s
 </section>
 <div class="src">Generated from <code>%(JNAME)s</code> and the built
 <code>Alignments/&lt;Module&gt;/&lt;Feature&gt;/pssms/</code> trees by <code>collect_registry.py</code> and
 <code>gen_registry.py</code>. Callable measured by running each feature&rsquo;s own PSSMs back against its own
 collection with <code>psiblast -in_pssm</code> and testing every sequence against that feature&rsquo;s
 <code>bit_cutoff</code> &mdash; self-recall, a floor rather than a validation. PSSM labels are cluster ids;
 files are <code>&lt;Module&gt;.&lt;Feature&gt;.&lt;cluster&gt;.pssm</code>.</div>
</div>
%(JS)s
""" % {
        "TAXON": html.escape(args.taxon), "DATE": html.escape(date),
        "CSS": css, "JS": js, "STAND": stand, "JNAME": html.escape(jname),
        "NPSSM": npssm, "NWITHP": len(withp), "NALL": nall,
        "REACH": reach, "ANNOS": annos, "PCT": pct, "NGAP": ngap, "NND": nnd,
        "NBUILT": counts["built"], "NDERIV": counts["derived"] + counts["special"],
        "SECTIONS": "\n".join(sections),
        "GAPNC": gapnc or '<tr><td colspan="3" class="dash">none</td></tr>',
        "GAPTF": gaptf or '<tr><td colspan="4" class="dash">none</td></tr>',
        "NTOOFEW": counts["too-few"],
        "GAPPROSE": notes.get("gaps") or (
            "<b>%d features have no collection at all</b> &mdash; the protein does not appear in "
            "BV-BRC under any product string the binning rules match, so there is nothing to align."
            % counts["no-collection"]),
        "DERIVED": ('<p style="margin-top:14px">%s</p>' % notes["derived"]) if notes.get("derived") else "",
    }
    with open(args.out, "w") as fh:
        fh.write(page)
    print("wrote %s (%d bytes)" % (args.out, len(page)))
    print("  %d PSSMs, %d/%d features with a profile, %d%% of %d sequences callable"
          % (npssm, len(withp), nall, pct, nseq))
    return 0


if __name__ == "__main__":
    sys.exit(main())
