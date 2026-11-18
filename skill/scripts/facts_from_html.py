#!/usr/bin/env python3
"""Recover a coverage-audit facts file from its rendered page.

37 of the 77 coverage audits have no facts file in the repo, so they cannot
be regenerated at all: gen_coverage_audit.py is the only supported way to
change one, and without its input the page can only be hand-edited, which
reports/README forbids. The numbers behind those pages came from
gen_facts.py over measurement files, but gen_facts also needs a hand-written
prose file, and for these modules that prose exists nowhere except the
rendered page. Re-deriving is therefore not possible; recovering is.

gen_coverage_audit.py emits a fully deterministic page from the facts, so
the mapping is invertible. This walks it backwards. The test that it worked
is not inspection: regenerate from the recovered facts and diff against the
page it came from. Byte-identical or it is wrong.

Figures are lifted at render time out of a named file, so a page carrying a
figure whose source is gone cannot be regenerated. Each embedded <figure> is
written beside the page as its own fragment and the block points there,
which makes the facts file self-contained.
"""
import argparse, html, json, os, re, sys

#  The inverse of gen_coverage_audit.esc, applied in reverse order: esc runs
#  html.escape, then unescapes the named-entity form, then unescapes a
#  whitelist of inline tags.
_TAGS = "br|i|b|em|strong|code|sub|sup|small"
_R_TAG = re.compile(r"<(/?(?:%s))\s*/?>" % _TAGS)
_R_ENT = re.compile(r"&(#\d{1,6}|#x[0-9A-Fa-f]{1,6}|[A-Za-z][A-Za-z0-9]{1,31});")


def unesc(r):
    """Recover x from esc(x). Verified by the caller re-running esc."""
    s = _R_TAG.sub(r"&lt;\1&gt;", r)
    s = _R_ENT.sub(r"&amp;\1;", s)
    return (s.replace("&amp;lt;", "&lt;").replace("&amp;gt;", "&gt;")
             .replace("&amp;quot;", '&quot;').replace("&amp;#x27;", "&#x27;")
             .replace("&amp;amp;", "&amp;"))


def _esc(x):
    """gen_coverage_audit.esc, copied so the round trip can be checked here."""
    t = re.compile(r"&lt;(/?(?:%s))\s*/?&gt;" % _TAGS)
    e = re.compile(r"&amp;(#\d{1,6}|#x[0-9A-Fa-f]{1,6}|[A-Za-z][A-Za-z0-9]{1,31});")
    return t.sub(r"<\1>", e.sub(r"&\1;", html.escape(str(x))))


def un(r):
    v = unesc(r)
    if _esc(v) != r:          # fall back to the literal where inversion drifts
        return r if _esc(r) == r else v
    return v


def one(pat, s, flags=re.S):
    m = re.search(pat, s, flags)
    return m.group(1) if m else None


def parse_tiles(frag):
    return [[un(v), un(l)] for v, l in re.findall(
        r'<div class="tile"><div class="v">(.*?)</div>'
        r'<div class="l">(.*?)</div></div>', frag, re.S)]


def parse_bars(frag):
    rows, unit = [], ""
    for lab, w, val, note in re.findall(
            r'<div class="bar"><div class="bl">(.*?)</div>'
            r'<div class="bt"><i style="width:([\d.]+)%"></i></div>'
            r'<div class="bv">(.*?)</div><div class="bn">(.*?)</div></div>',
            frag, re.S):
        rows.append([un(lab), un(val), float(w), un(note)])
    #  bars() appends a shared unit to every value; recover it only if every
    #  row ends with the same non-numeric tail.
    tails = {r[1][-1] for r in rows if r[1]}
    if len(tails) == 1 and tails and not list(tails)[0].isdigit():
        t = list(tails)[0]
        if all(len(r[1]) > 1 for r in rows):
            unit = t
            for r in rows:
                r[1] = r[1][:-1]
    return rows, unit


def parse_table(frag):
    headers = [un(h) for h in re.findall(r"<th>(.*?)</th>", frag, re.S)]
    body = one(r"<tbody>(.*)</tbody>", frag) or ""
    rows = []
    for tr in re.findall(r"<tr>(.*?)</tr>", body, re.S):
        rows.append([c if c.startswith("<") else un(c)
                     for c in re.findall(r"<td>(.*?)</td>", tr, re.S)])
    return headers, rows


def blocks_of(body, stem, outdir, figs):
    """Walk a section body in document order and type each block."""
    out, i = [], 0
    pat = re.compile(
        r'<div class="tiles">|<div class="bars">|<div class="scroll">'
        r'|<pre class="mono">|<p class="note">|<div class="verdict'
        r'|<figure')
    while True:
        m = pat.search(body, i)
        if not m:
            break
        k, j = m.group(0), m.start()
        if k == '<div class="tiles">':
            e = body.index("</div>", body.index('</div></div>', j)) + 6
            e = body.index("</div>", j + len('<div class="tiles">'))
            e = _close_div(body, j)
            out.append({"type": "tiles", "items": parse_tiles(body[j:e])})
        elif k == '<div class="bars">':
            e = _close_div(body, j)
            r, u = parse_bars(body[j:e])
            b = {"type": "bars", "rows": r}
            if u:
                b["unit"] = u
            out.append(b)
        elif k == '<div class="scroll">':
            e = _close_div(body, j)
            frag = body[j:e]
            #  Nine parvovirus pages have a genome map spliced in between
            #  </table> and the closing </div>, so the figure sits INSIDE the
            #  table's overflow-x:auto wrapper. That was my splice, and it is
            #  wrong: a figure block belongs after the scroll container, not
            #  within it. Hoist it, which both recovers it and re-seats it.
            hoisted = None
            fi = frag.find("<figure")
            if fi >= 0:
                fj = frag.rfind("</figure>") + len("</figure>")
                hoisted = frag[fi:fj]
                frag = frag[:fi] + frag[fj:]
            h, r = parse_table(frag)
            out.append({"type": "table", "headers": h, "rows": r})
            if hoisted is not None:
                name = _existing(hoisted, stem, outdir)
                if name is None:
                    name = ("%s-saturation-figure.html" % stem
                            if "chartbox" in hoisted
                            else "%s-figure%d.html" % (stem, len(figs) + 1))
                    figs.append((name, hoisted))
                out.append({"type": "figure", "from": name})
        elif k == '<pre class="mono">':
            e = body.index("</pre>", j) + 6
            out.append({"type": "pre",
                        "text": un(body[j + len(k):e - 6])})
        elif k == '<p class="note">':
            e = body.index("</p>", j) + 4
            out.append({"type": "prose", "html": body[j + len(k):e - 4]})
        elif k.startswith('<div class="verdict'):
            e = _close_div(body, j)
            cls = one(r'<div class="verdict([^"]*)">', body[j:e])
            kind = {"": "ok", " warn": "warn", " stop": "stop"}[cls]
            inner = one(r'<div class="verdict[^"]*">(.*)</div>$', body[j:e])
            out.append({"type": "verdict", "kind": kind, "html": inner})
        else:                                   # <figure>
            e = body.index("</figure>", j) + len("</figure>")
            frag = body[j:e]
            #  Prefer an existing file that already holds this exact figure --
            #  the genome maps are published in their own right, and pointing
            #  at a second copy would leave two files to keep in step.
            name = _existing(frag, stem, outdir)
            if name is None:
                #  The saturation curve is the only other figure in the series;
                #  name it for what it is rather than by position.
                name = ("%s-saturation-figure.html" % stem
                        if "chartbox" in frag
                        else "%s-figure%d.html" % (stem, len(figs) + 1))
                figs.append((name, frag))
            out.append({"type": "figure", "from": name})
        i = e
    return out


def _existing(frag, stem, outdir):
    """A file already beside the page holding exactly this <figure>, or None."""
    for cand in ("%s-genome-map.html" % stem,
                 "%s-saturation-figure.html" % stem):
        p = os.path.join(outdir, cand)
        if not os.path.exists(p):
            continue
        h = open(p).read()
        i, j = h.find("<figure"), h.rfind("</figure>")
        if i >= 0 and j >= 0 and h[i:j + len("</figure>")] == frag:
            return cand
    return None


def _close_div(s, start):
    """Index just past the <div> opened at `start`, matching nesting."""
    i, depth = start, 0
    while True:
        m = re.compile(r"<div\b|</div>").search(s, i)
        if not m:
            raise ValueError("unbalanced div")
        if m.group(0) == "</div>":
            depth -= 1
            if depth == 0:
                return m.end()
        else:
            depth += 1
        i = m.end()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("page")
    ap.add_argument("--out", required=True)
    ap.add_argument("--figdir", default=None,
                    help="where figure fragments are written (default: page's dir)")
    a = ap.parse_args()

    h = open(a.page).read()
    stem = os.path.basename(a.page)[:-len("-coverage-audit.html")]
    figdir = a.figdir or os.path.dirname(os.path.abspath(a.page))

    F = {}
    F["title"] = un(one(r"<title>(.*?)</title>", h))
    F["eyebrow"] = un(one(r'<div class="eyebrow">(.*?)</div>', h))
    F["h1"] = un(one(r"<h1>(.*?)</h1>", h))
    sub = one(r'<h2 style="margin:8px 0 0;[^"]*">(.*?)</h2>', h)
    if sub is not None:
        F["subhead"] = un(sub)
    F["standfirst"] = one(r'<p class="stand">(.*?)</p>', h)
    hdr = one(r"<header>(.*?)</header>", h)
    F["headline_tiles"] = parse_tiles(hdr)
    F["source_note"] = one(r'<p class="src">(.*?)</p></div>\s*$', h)

    figs, secs = [], []
    for sec in re.findall(r"<section>(.*?)</section>", h, re.S):
        num = one(r'<h2><span class="num">(.*?)</span>', sec)
        title = one(r'<h2><span class="num">.*?</span>(.*?)</h2>', sec)
        lede = one(r'<p class="lede">(.*?)</p>', sec)
        body = sec[sec.index('</p>', sec.index('<p class="lede">')) + 4:]
        s = {"num": un(num), "title": un(title),
             "blocks": blocks_of(body, stem, figdir, figs)}
        if lede:
            s["lede"] = lede
        secs.append(s)
    F["sections"] = secs

    for name, frag in figs:
        open(os.path.join(figdir, name), "w").write(frag + "\n")

    json.dump(F, open(a.out, "w"), indent=1, sort_keys=True)
    print("%-26s %d sections, %d figure(s)" % (stem, len(secs), len(figs)))


main()
