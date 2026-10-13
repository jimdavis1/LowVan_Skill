#!/usr/bin/env python3
"""Lint a rendered artifact page before it is published.

Run it on every page, every time. These are not hypothetical defects: each
one below shipped to a published URL at least once, and two of them shipped
twice, because the only check being applied was reading the page.

    python3 check_artifact_html.py artifacts/*.html
    python3 check_artifact_html.py --quiet artifacts/*.html && publish

Exit status is 1 if anything is reported, so it can gate a publish step.

WHAT IT CHECKS, AND WHY EACH ONE EXISTS

1. DOUBLE-ESCAPED ENTITY -- "&amp;times;" renders as the literal text
   "&times;". A facts file writes entities throughout because most of what it
   emits is prose, which renderers insert raw; but titles, eyebrows, tile
   values, table headers and table cells go through html.escape(), which
   turns the entity into its own source text. Shipped on the Rotavirus and
   Orbivirus coverage audits as "10.8&times;" in a headline tile,
   "LowVan &middot; Orbivirus module" in the eyebrow, and "contig &ge; 5 kb"
   in the first table.

2. PLACEHOLDER LEFT IN TEXT -- "%s", "%d", "nan", "None", "(not in JSON)",
   "undefined" standing alone. A format string that was never filled, or a
   lookup that missed. "undefined CDS" is a real BV-BRC annotation and is not
   flagged; a bare "undefined" is.

3. ANGLE-BRACKET PLACEHOLDER IN PROSE -- "the explicit VP<n> forms" reads as
   broken markup to anyone who does not already know it is a placeholder.
   Write the range instead: VP1-VP7. The same construct inside a <code> span
   and surrounded by path separators -- Alignments/<Module>/<Feature>/ -- is
   the universal convention and is reported as a note, not a failure.

4. UNKNOWN TAG -- anything outside the set a generated page is expected to
   emit, which catches (3) in its parsed form plus any stray construct.

5. EMPTY SECTION -- a <section> with no text, meaning a prose key was missing
   from the facts file and nothing said so.
"""
import argparse, html, os, re, sys

ENTITY = re.compile(r"&amp;(#\d{1,6}|#x[0-9A-Fa-f]{1,6}|[A-Za-z][A-Za-z0-9]{1,31});")
PLACEHOLDER = re.compile(
    r"%[sdfr]\b|\bnan\b|\bNone\b|\(not in JSON\)|(?<![\w-])undefined(?!\s+CDS)\b")
ANGLE = re.compile(r"&lt;[A-Za-z_][A-Za-z0-9_ -]{0,20}&gt;")
TAGS_OK = {
    "a","abbr","aside","b","blockquote","br","button","canvas","caption","circle",
    "clippath","code","col","colgroup","dd","defs","details","div","dl","dt","em",
    "figcaption","figure","footer","g","h1","h2","h3","h4","header","hr","i","img",
    "input","label","li","line","lineargradient","link","main","mark","mask","nav",
    "ol","p","path","polygon","polyline","pre","rect","script","section","small",
    "span","stop","strong","style","sub","summary","sup","svg","table","tbody","td",
    "text","textpath","tspan","th","thead","time","title","tr","ul","use",
}


def check(path):
    raw = open(path, errors="replace").read()
    notes = []
    body = re.sub(r"<style.*?</style>|<script.*?</script>|<!--.*?-->", " ", raw, flags=re.S)
    text = re.sub(r"<[^>]+>", " ", body)
    out = []

    ent = sorted(set(ENTITY.findall(body)))
    if ent:
        out.append(("double-escaped entity", ["&%s;" % e for e in ent][:8]))

    ph = sorted({m.group(0) for m in PLACEHOLDER.finditer(text)})
    if ph:
        out.append(("placeholder left in text", ph[:8]))

    #  A path template inside a <code> span -- Alignments/<Module>/<Feature>/ --
    #  is the universal convention and the slashes make it obviously a
    #  template. The same construct in a SENTENCE is what misleads: "the
    #  explicit VP<n> forms" reads as broken markup. Flag prose, note code.
    incode = set()
    for m in re.finditer(r"<code>(.*?)</code>", body, re.S):
        incode.update(ANGLE.findall(m.group(1)))
    ang = sorted(set(ANGLE.findall(body)))
    prose_ang = [a for a in ang if a not in incode]
    if prose_ang:
        out.append(("angle-bracket placeholder in prose",
                    [html.unescape(a) for a in prose_ang][:8]))
    elif ang:
        notes.append(("angle placeholder in a code span (allowed)",
                      [html.unescape(a) for a in sorted(incode)][:6]))

    tags = {t.lower() for t in re.findall(r"</?([A-Za-z][A-Za-z0-9]*)[\s/>]", body)}
    unknown = sorted(tags - TAGS_OK)
    if unknown:
        out.append(("unknown tag", unknown[:8]))

    #  A section always has its number and title in the <h2>, so "no text at
    #  all" never fires. What matters is a section with a heading and NO BODY:
    #  a prose key missing from the facts file, which renders as a lone
    #  heading and says nothing about itself.
    empty = []
    for m in re.finditer(r"<section[^>]*>(.*?)</section>", body, re.S):
        inner = re.sub(r"<h[1-6][^>]*>.*?</h[1-6]>", " ", m.group(1), flags=re.S)
        if not re.sub(r"<[^>]+>|\s|&[A-Za-z#0-9]+;", "", inner):
            empty.append(m.start())
    if empty:
        out.append(("empty section", ["%d of them" % len(empty)]))
    return out, notes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pages", nargs="+")
    ap.add_argument("--quiet", action="store_true", help="print only failures")
    a = ap.parse_args()
    bad = 0
    for p in a.pages:
        if not os.path.exists(p):
            print("  %-44s MISSING" % os.path.basename(p)); bad += 1; continue
        f, notes = check(p)
        if f:
            bad += 1
            print("  %-44s %s" % (os.path.basename(p), f[0][0]))
            for name, items in f:
                print("      %-26s %s" % (name, ", ".join(items)))
        elif not a.quiet:
            print("  %-44s clean" % os.path.basename(p))
        if notes and not a.quiet:
            for name, items in notes:
                print("      note: %-20s %s" % (name, ", ".join(items)))
    if bad:
        print("\n  %d page(s) with a problem -- do NOT publish until these are fixed." % bad)
    elif not a.quiet:
        print("\n  all %d page(s) clean" % len(a.pages))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
