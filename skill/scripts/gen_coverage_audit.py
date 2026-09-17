#!/usr/bin/env python3
"""Render the coverage audit page for a module.

The fourth report, and the one the kit had no generator for: Rhabdoviridae and
Togaviridae each shipped a hand-authored `<taxon>-coverage-audit.html` and
nothing produced them, so Matonaviridae never got one and the Hepeviridae build
did not know to make one either.

Driven by an explicit facts file rather than by scraping a working directory,
because the numbers on this page come from six different tools across steps 10,
11 and 12, and a page that silently picks up whichever files happen to be lying
around is worse than one whose inputs are written down. Every figure in the
output appears in the facts JSON, so the page is auditable against the build
notes.

    python3 gen_coverage_audit.py --facts hepeviridae_audit.json \
            --out hepeviridae-coverage-audit.html

Schema: see the two shipped facts files in example-hepeviridae/ and
example-matonaviridae/. Palette is the report series' validated categorical
pair (#2a78d6/#eb6834 light, #3987e5/#d95926 dark).
"""
import json, argparse, html

CSS = """
  :root {
    color-scheme: light;
    --plane:#f9f9f7; --surface:#fcfcfb; --ink:#0b0b0b; --ink2:#52514e; --muted:#898781;
    --grid:#e1e0d9; --base:#c3c2b7; --hair:rgba(11,11,11,.10);
    --built:#2a78d6; --missed:#eb6834; --warn:#b07d18;
    --serif:"IBM Plex Serif",Georgia,serif;
    --sans:"IBM Plex Sans",system-ui,-apple-system,"Segoe UI",sans-serif;
    --mono:"IBM Plex Mono",ui-monospace,Menlo,monospace;
  }
  @media (prefers-color-scheme:dark) { :root:not([data-theme="light"]) {
    color-scheme: dark;
    --plane:#0d0d0d; --surface:#1a1a19; --ink:#fff; --ink2:#c3c2b7; --muted:#898781;
    --grid:#2c2c2a; --base:#4a4a46; --hair:rgba(255,255,255,.10);
    --built:#3987e5; --missed:#d95926; --warn:#d3a445;
  } }
  :root[data-theme="dark"] {
    color-scheme: dark;
    --plane:#0d0d0d; --surface:#1a1a19; --ink:#fff; --ink2:#c3c2b7; --muted:#898781;
    --grid:#2c2c2a; --base:#4a4a46; --hair:rgba(255,255,255,.10);
    --built:#3987e5; --missed:#d95926; --warn:#d3a445;
  }
  * { box-sizing:border-box; }
  body { margin:0; background:var(--plane); color:var(--ink);
         font-family:var(--sans); font-size:15px; line-height:1.62;
         -webkit-font-smoothing:antialiased; }
  .wrap { max-width:860px; margin:0 auto; padding:0 26px 90px; }
  header { padding:54px 0 30px; border-bottom:2px solid var(--ink); margin-bottom:34px; }
  .eyebrow { font-family:var(--mono); font-size:11px; letter-spacing:.13em;
              text-transform:uppercase; color:var(--muted); margin-bottom:16px; }
  h1 { font-family:var(--serif); font-weight:600; font-size:clamp(30px,5.2vw,46px);
        line-height:1.1; margin:0 0 16px; letter-spacing:-.015em; text-wrap:balance; }
  .stand { font-size:17px; line-height:1.58; color:var(--ink2); max-width:64ch; margin:0; }
  .stand b { color:var(--ink); font-weight:600; }
  .tiles { display:grid; grid-template-columns:repeat(auto-fit,minmax(146px,1fr));
            gap:1px; background:var(--hair); border:1px solid var(--hair); margin:30px 0 0; }
  .tile { background:var(--surface); padding:15px 16px 14px; }
  .tile .v { font-family:var(--mono); font-size:25px; font-weight:500; line-height:1.1;
              font-variant-numeric:tabular-nums; }
  .tile .l { font-size:11.5px; color:var(--muted); margin-top:5px; line-height:1.35; }
  section { margin:46px 0 0; }
  h2 { font-family:var(--serif); font-weight:600; font-size:23px; line-height:1.25;
        margin:0 0 4px; letter-spacing:-.01em; text-wrap:balance; }
  h2 .num { font-family:var(--mono); font-size:12px; color:var(--muted);
             letter-spacing:.08em; display:block; margin-bottom:7px; font-weight:400; }
  h3 { font-family:var(--sans); font-weight:600; font-size:14px; letter-spacing:.02em;
        margin:30px 0 8px; color:var(--ink); }
  p { max-width:68ch; margin:14px 0; }
  .lede { color:var(--ink2); }
  code, .mono { font-family:var(--mono); font-size:.89em; }
  code { background:var(--surface); border:1px solid var(--hair); padding:.08em .32em; border-radius:2px; }
  .scroll { overflow-x:auto; margin:20px 0; border:1px solid var(--hair); background:var(--surface); }
  table { border-collapse:collapse; width:100%; font-size:13px; }
  th { text-align:left; font-family:var(--mono); font-weight:500; font-size:10.5px;
        letter-spacing:.09em; text-transform:uppercase; color:var(--muted);
        padding:10px 11px; border-bottom:1px solid var(--base); white-space:nowrap; }
  td { padding:8px 11px; border-bottom:1px solid var(--hair); vertical-align:top; }
  tr:last-child td { border-bottom:none; }
  td.k { font-family:var(--mono); font-weight:500; white-space:nowrap; }
  td.n { font-family:var(--mono); font-variant-numeric:tabular-nums; white-space:nowrap; }
  td.a { color:var(--ink2); }
  td.src { font-size:11.5px; color:var(--muted); line-height:1.5; }
  .dim { color:var(--muted); }
  figure { margin:22px 0; }
  figure svg { width:100%; height:auto; display:block; }
  figcaption { font-size:12.5px; color:var(--muted); margin-top:10px; max-width:68ch; line-height:1.5; }
  .legend { display:flex; flex-wrap:wrap; gap:16px; font-size:12px; color:var(--ink2);
             font-family:var(--mono); margin-top:12px; }
  .legend i { width:9px; height:9px; display:inline-block; border-radius:2px; margin-right:6px;
               vertical-align:baseline; }
  .verdict { border-left:3px solid var(--built); background:var(--surface);
              padding:15px 18px; margin:22px 0; }
  .verdict.warn { border-left-color:var(--warn); }
  .verdict.stop { border-left-color:var(--missed); }
  .verdict h4 { margin:0 0 6px; font-size:13px; font-weight:600; letter-spacing:.02em; }
  .verdict p { margin:0; font-size:14px; color:var(--ink2); max-width:none; }
  .verdict p + p { margin-top:9px; }
  pre { font-family:var(--mono); font-size:12px; line-height:1.55; background:var(--surface);
         border:1px solid var(--hair); padding:14px 16px; overflow-x:auto; margin:18px 0;
         color:var(--ink2); }
  ul { max-width:68ch; padding-left:19px; }
  li { margin:7px 0; }
  footer { margin-top:60px; padding-top:20px; border-top:1px solid var(--base);
            font-size:12.5px; color:var(--muted); }
  a { color:var(--built); }
  a:focus-visible, [tabindex]:focus-visible { outline:2px solid var(--built); outline-offset:2px; }
  @media (prefers-reduced-motion:reduce) { * { animation:none!important; transition:none!important; } }

  .bars { display:flex; flex-direction:column; gap:9px; margin:18px 0 0; }
  .bar { display:grid; grid-template-columns:132px 1fr 62px; gap:11px; align-items:center; }
  .bl { font-size:13px; color:var(--ink2); }
  .bt { height:9px; background:var(--grid); position:relative; overflow:hidden; }
  .bt i { display:block; height:100%; background:var(--built); }
  .bv { font-family:var(--mono); font-size:13px; font-variant-numeric:tabular-nums;
        text-align:right; color:var(--ink); }
  .bn { grid-column:2/4; font-size:11.5px; color:var(--muted); margin-top:-4px; }
  pre.mono { font-family:var(--mono); font-size:12.5px; line-height:1.55; color:var(--ink2);
             background:var(--surface); border:1px solid var(--hair); padding:14px 15px;
             overflow-x:auto; margin:18px 0 0; }
  table { border-collapse:collapse; width:100%; font-size:13.5px; margin:18px 0 0; }
  thead th { text-align:left; font-size:11px; letter-spacing:.07em; text-transform:uppercase;
             color:var(--muted); font-weight:500; border-bottom:1px solid var(--grid);
             padding:0 12px 7px 0; white-space:nowrap; }
  tbody td { padding:7px 12px 7px 0; border-bottom:1px solid var(--hair); color:var(--ink2);
             vertical-align:top; }
  tbody td:first-child { color:var(--ink); font-family:var(--mono); font-size:12.5px;
                         white-space:nowrap; }
  .scroll { overflow-x:auto; }
  @media (max-width:560px) { .bar { grid-template-columns:104px 1fr 56px; } }
"""

def esc(x): return html.escape(str(x))

def tiles(items):
    out = ['<div class="tiles">']
    for v, l in items:
        out.append('<div class="tile"><div class="v">%s</div><div class="l">%s</div></div>'
                   % (esc(v), esc(l)))
    out.append('</div>')
    return "".join(out)

def bars(rows, unit=""):
    """rows: (label, value, pct, note). pct drives the bar width."""
    out = ['<div class="bars">']
    for lab, val, pct, note in rows:
        w = max(0.0, min(100.0, float(pct)))
        out.append(
          '<div class="bar"><div class="bl">%s</div>'
          '<div class="bt"><i style="width:%.2f%%"></i></div>'
          '<div class="bv">%s%s</div><div class="bn">%s</div></div>'
          % (esc(lab), w, esc(val), esc(unit), esc(note or "")))
    out.append('</div>')
    return "".join(out)

def table(headers, rows, scroll=True):
    out = ['<div class="scroll">' if scroll else '', '<table><thead><tr>']
    out += ['<th>%s</th>' % esc(h) for h in headers]
    out.append('</tr></thead><tbody>')
    for r in rows:
        out.append('<tr>' + "".join('<td>%s</td>' % (c if str(c).startswith('<') else esc(c))
                                    for c in r) + '</tr>')
    out.append('</tbody></table>')
    out.append('</div>' if scroll else '')
    return "".join(out)

def section(num, title, lede, body):
    return ('<section><h2><span class="num">%s</span>%s</h2>'
            '<p class="lede">%s</p>%s</section>'
            % (esc(num), esc(title), lede, body))

def verdict(kind, text):
    return '<div class="verdict%s">%s</div>' % (
        "" if kind == "ok" else (" warn" if kind == "warn" else " stop"), text)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--facts", required=True)
    ap.add_argument("--out", default="coverage-audit.html")
    a = ap.parse_args()
    F = json.load(open(a.facts))

    S = []
    for s in F["sections"]:
        body = []
        for b in s.get("blocks", []):
            t = b["type"]
            if   t == "tiles":   body.append(tiles(b["items"]))
            elif t == "bars":    body.append(bars(b["rows"], b.get("unit", "")))
            elif t == "table":   body.append(table(b["headers"], b["rows"]))
            elif t == "pre":     body.append('<pre class="mono">%s</pre>' % esc(b["text"]))
            elif t == "prose":   body.append('<p class="note">%s</p>' % b["html"])
            elif t == "verdict": body.append(verdict(b.get("kind", "ok"), b["html"]))
        S.append(section(s["num"], s["title"], s["lede"], "".join(body)))

    page = ("<title>%s</title>\n<style>%s</style>\n"
            '<div class="wrap"><header>'
            '<div class="eyebrow">%s</div><h1>%s</h1><p class="stand">%s</p>%s</header>'
            "%s"
            '<p class="src">%s</p></div>'
            % (esc(F["title"]), CSS, esc(F["eyebrow"]), esc(F["h1"]),
               F["standfirst"], tiles(F["headline_tiles"]), "".join(S), F["source_note"]))
    open(a.out, "w").write(page)
    print("wrote %s (%d bytes, %d sections)" % (a.out, len(page), len(S)))

main()
