#!/usr/bin/env python3
"""Render the lowvan-module skill's documentation as one browsable page.

Embeds each markdown file verbatim and renders it client-side with marked, so
the page stays in step with the sources -- rerun this after editing any of them.

    python3 gen_handbook.py --skill ~/.claude/skills/lowvan-module --out handbook.html
"""

import argparse
import json
import os
import re
import sys

# SKILL.md is the spine; the references are depth behind particular steps, so
# the nav says which step each one belongs to rather than listing them flat.
REFS = [
    ("inputs.md",                "1", "Taking delivery of the BV-BRC dump"),
    ("module-partitioning.md",   "2", "Deciding what the modules are"),
    ("annotation-triage.md",     "2", "Binning strings into collections"),
    ("annotation-vocabulary.md", "6", "Choosing annotation strings"),
    ("pipeline.md",              "3", "Clustering, alignment, PSSMs"),
    ("curation.md",              "4", "Curating the alignments"),
    ("json-schema.md",           "6", "The Viral_PSSM.json block"),
    ("install-and-test.md",     "8-9", "Installing and scoring"),
    ("artifacts.md",            "10", "The two reports"),
]


def h2s(md):
    """Top-level sections, for the sidebar's second level."""
    out = []
    fence = False
    for line in md.split("\n"):
        if line.startswith("```"):
            fence = not fence
        elif not fence and line.startswith("## "):
            out.append(line[3:].strip())
    return out


def slug(s):
    return re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skill", required=True)
    ap.add_argument("--out", default="handbook.html")
    args = ap.parse_args()

    skill = os.path.abspath(os.path.expanduser(args.skill))
    docs = []

    p = os.path.join(skill, "SKILL.md")
    md = open(p, encoding="utf8").read()
    # strip the YAML frontmatter; it is metadata for the loader, not prose
    body = re.sub(r"\A---\n.*?\n---\n", "", md, flags=re.S)
    docs.append({"id": "skill", "file": "SKILL.md", "step": "",
                 "title": "The workflow", "sub": "start here",
                 "md": body, "h2": h2s(body)})

    for fname, step, sub in REFS:
        fp = os.path.join(skill, "references", fname)
        if not os.path.exists(fp):
            print("  missing: %s" % fp, file=sys.stderr)
            continue
        m = open(fp, encoding="utf8").read()
        first = next((l[2:].strip() for l in m.split("\n") if l.startswith("# ")), fname)
        docs.append({"id": slug(fname[:-3]), "file": "references/" + fname,
                     "step": step, "title": first, "sub": sub,
                     "md": m, "h2": h2s(m)})

    nscripts = len([f for f in os.listdir(os.path.join(skill, "scripts"))
                    if f.endswith(".py")]) if os.path.isdir(os.path.join(skill, "scripts")) else 0
    nlines = sum(len(d["md"].split("\n")) for d in docs)

    nav = []
    for d in docs:
        subs = "".join(
            '<a class="h2" href="#%s--%s">%s</a>' % (d["id"], slug(h), h) for h in d["h2"])
        nav.append(
            '<div class="navgroup"><a class="doc" href="#%s" data-doc="%s">'
            '%s<span class="dt">%s</span><span class="ds">%s</span></a>'
            '<div class="subs">%s</div></div>'
            % (d["id"], d["id"],
               ('<span class="step">%s</span>' % d["step"]) if d["step"] else
               '<span class="step spine">&#9679;</span>',
               d["title"], d["sub"], subs))

    payload = json.dumps([{k: d[k] for k in ("id", "file", "title", "md", "h2")}
                          for d in docs])

    # NB token replacement, not %-formatting: the CSS is full of literal "%"
    html = TEMPLATE
    for k, v in (("@NAV@", "\n".join(nav)), ("@DATA@", payload),
                 ("@NDOCS@", str(len(docs))), ("@NLINES@", "{:,}".format(nlines)),
                 ("@NSCRIPTS@", str(nscripts))):
        html = html.replace(k, v)
    with open(args.out, "w", encoding="utf8") as fh:
        fh.write(html)
    print("wrote %s (%d bytes, %d documents, %s lines)"
          % (args.out, len(html), len(docs), "{:,}".format(nlines)))
    return 0


TEMPLATE = r"""<title>LowVan Module Handbook</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Serif:wght@500;600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<script src="https://cdnjs.cloudflare.com/ajax/libs/marked/12.0.2/marked.min.js"></script>
<style>
:root{
  color-scheme:light;
  --plane:#f9f9f7; --surface:#fcfcfb; --sink:#f2f2ee;
  --ink:#111110; --ink2:#4f4e4a; --muted:#8a8880;
  --grid:#e3e2db; --base:#c6c5ba; --hair:rgba(17,17,16,.09);
  --accent:#2a78d6; --accent-soft:rgba(42,120,214,.10);
  --serif:"IBM Plex Serif",Georgia,serif;
  --sans:"IBM Plex Sans",system-ui,-apple-system,"Segoe UI",sans-serif;
  --mono:"IBM Plex Mono",ui-monospace,Menlo,monospace;
}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
  color-scheme:dark;
  --plane:#0d0d0c; --surface:#171716; --sink:#121211;
  --ink:#f7f7f5; --ink2:#c3c2b9; --muted:#8a8880;
  --grid:#2b2b28; --base:#3a3a36; --hair:rgba(255,255,255,.09);
  --accent:#5c9dea; --accent-soft:rgba(92,157,234,.14);
}}
:root[data-theme="dark"]{
  color-scheme:dark;
  --plane:#0d0d0c; --surface:#171716; --sink:#121211;
  --ink:#f7f7f5; --ink2:#c3c2b9; --muted:#8a8880;
  --grid:#2b2b28; --base:#3a3a36; --hair:rgba(255,255,255,.09);
  --accent:#5c9dea; --accent-soft:rgba(92,157,234,.14);
}
*{box-sizing:border-box}
body{background:var(--plane);color:var(--ink);font-family:var(--sans);line-height:1.62;margin:0}

.masthead{border-bottom:1px solid var(--base);background:var(--surface)}
.mi{max-width:1240px;margin:0 auto;padding:22px 26px 20px;display:flex;
    align-items:flex-end;justify-content:space-between;gap:24px;flex-wrap:wrap}
.eyebrow{font-family:var(--mono);font-size:10.5px;letter-spacing:.13em;
         text-transform:uppercase;color:var(--muted)}
h1.top{font-family:var(--serif);font-weight:600;font-size:26px;margin:6px 0 0;letter-spacing:-.01em}
.blurb{color:var(--ink2);font-size:14px;margin:7px 0 0;max-width:62ch}
.facts{display:flex;gap:26px;font-family:var(--mono);font-size:11.5px;color:var(--muted);white-space:nowrap}
.facts b{display:block;font-family:var(--sans);font-size:19px;font-weight:600;color:var(--ink);
         line-height:1.2;font-variant-numeric:tabular-nums}

.shell{max-width:1240px;margin:0 auto;padding:0 26px;display:grid;
       grid-template-columns:292px minmax(0,1fr);gap:38px;align-items:start}

nav{position:sticky;top:0;max-height:100vh;overflow-y:auto;padding:24px 0 60px;
    border-right:1px solid var(--grid)}
.navgroup{margin-bottom:3px}
a.doc{display:grid;grid-template-columns:26px 1fr;gap:0 9px;text-decoration:none;
      padding:7px 12px 7px 4px;border-radius:5px;align-items:baseline}
a.doc:hover{background:var(--sink)}
a.doc.on{background:var(--accent-soft)}
.step{font-family:var(--mono);font-size:10px;color:var(--muted);text-align:right;
      padding-top:2px;grid-row:span 2}
.step.spine{color:var(--accent)}
.dt{font-size:13.5px;font-weight:500;color:var(--ink);line-height:1.35}
a.doc.on .dt{color:var(--accent)}
.ds{font-size:11.5px;color:var(--muted);grid-column:2}
.subs{display:none;padding:2px 0 8px 39px;border-left:1px solid var(--grid);margin-left:16px}
.navgroup.open .subs{display:block}
a.h2{display:block;font-size:12px;color:var(--ink2);text-decoration:none;padding:3.5px 8px;
     border-radius:4px;line-height:1.4}
a.h2:hover{color:var(--accent);background:var(--sink)}

main{padding:30px 0 110px;min-width:0}
article{background:var(--surface);border:1px solid var(--grid);border-radius:9px;
        padding:34px 40px 40px;margin-bottom:26px}
.fname{font-family:var(--mono);font-size:11px;color:var(--muted);margin-bottom:14px;
       padding-bottom:12px;border-bottom:1px solid var(--grid)}
article h1{font-family:var(--serif);font-weight:600;font-size:25px;margin:0 0 4px;
           text-wrap:balance;letter-spacing:-.01em}
article h2{font-family:var(--serif);font-weight:600;font-size:18.5px;
           margin:34px 0 10px;padding-top:16px;border-top:1px solid var(--grid);text-wrap:balance}
article h1+h2{border-top:0;padding-top:0;margin-top:22px}
article h3{font-size:14.5px;font-weight:600;margin:22px 0 6px;color:var(--ink)}
article p,article li{font-size:14.6px;color:var(--ink2);max-width:74ch}
article strong{color:var(--ink);font-weight:600}
article a{color:var(--accent)}
article ul,article ol{padding-left:22px}
article li{margin:4px 0}
article code{font-family:var(--mono);font-size:12.5px;background:var(--sink);
             border:1px solid var(--hair);border-radius:4px;padding:1px 5px;color:var(--ink)}
article pre{background:var(--sink);border:1px solid var(--grid);border-radius:7px;
            padding:14px 16px;overflow-x:auto;margin:14px 0}
article pre code{background:none;border:0;padding:0;font-size:12.4px;line-height:1.6;color:var(--ink2)}
article blockquote{margin:14px 0;padding:2px 0 2px 16px;border-left:3px solid var(--accent);
                   color:var(--ink2);font-style:italic}
.tw{overflow-x:auto;margin:14px 0;border:1px solid var(--grid);border-radius:7px}
article table{border-collapse:collapse;width:100%;font-size:13.4px}
article th{background:var(--sink);text-align:left;font-weight:600;color:var(--ink);
           padding:9px 13px;border-bottom:1px solid var(--base);white-space:nowrap}
article td{padding:8px 13px;border-bottom:1px solid var(--grid);color:var(--ink2);vertical-align:top}
article tr:last-child td{border-bottom:0}
article hr{border:0;border-top:1px solid var(--grid);margin:26px 0}

.foot{max-width:1240px;margin:0 auto;padding:0 26px 60px;color:var(--muted);font-size:12.5px}
@media (max-width:900px){
  .shell{grid-template-columns:1fr;gap:0}
  nav{position:static;max-height:none;border-right:0;border-bottom:1px solid var(--grid);
      padding-bottom:16px}
  .subs{display:none!important}
  article{padding:24px 20px 28px}
}
@media (prefers-reduced-motion:reduce){*{scroll-behavior:auto!important}}
html{scroll-behavior:smooth}
</style>

<div class="masthead"><div class="mi">
  <div>
    <div class="eyebrow">Claude Code skill &middot; lowvan-module</div>
    <h1 class="top">LowVan Module Handbook</h1>
    <p class="blurb">How to build a LowVan annotation module for a viral taxon: triage the BV-BRC
    strings, cluster and curate the alignments, build the PSSMs, write the JSON, then install it
    into a Viral_Annotation checkout and run the annotator against held-out genomes.</p>
  </div>
  <div class="facts">
    <div><b>@NDOCS@</b>documents</div>
    <div><b>@NLINES@</b>lines</div>
    <div><b>@NSCRIPTS@</b>scripts</div>
  </div>
</div></div>

<div class="shell">
  <nav>@NAV@</nav>
  <main id="main"></main>
</div>
<div class="foot">Rendered from the skill's own markdown by <code>gen_handbook.py</code>.
The files on disk are the source of truth &mdash; rerun it after editing any of them.</div>

<script>
const DOCS = @DATA@;
const slug = s => s.toLowerCase().replace(/[^a-z0-9]+/g,"-").replace(/^-|-$/g,"");

// If the CDN is unreachable the page must still be readable, so fall back to
// the raw markdown rather than rendering nothing.
const HAVE_MARKED = typeof marked !== "undefined";
if (HAVE_MARKED) marked.use({ mangle:false, headerIds:false });
const esc = s => s.replace(/[&<>]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;"}[c]));

const main = document.getElementById("main");
main.innerHTML = DOCS.map(d => {
  if (!HAVE_MARKED) {
    return `<article id="${d.id}"><div class="fname">${d.file}</div>` +
           `<pre><code>${esc(d.md)}</code></pre></article>`;
  }
  let html = marked.parse(d.md);
  // anchor every h2 so the sidebar can jump into a document
  html = html.replace(/<h2>(.*?)<\/h2>/g,
    (m,t) => `<h2 id="${d.id}--${slug(t.replace(/<[^>]+>/g,""))}">${t}</h2>`);
  // wide tables scroll inside their own box rather than the page
  html = html.replace(/<table>[\s\S]*?<\/table>/g, m => `<div class="tw">${m}</div>`);
  return `<article id="${d.id}"><div class="fname">${d.file}</div>${html}</article>`;
}).join("");

// open the nav group for whichever document is in view, and mark it active
const groups = [...document.querySelectorAll(".navgroup")];
const byId = {};
groups.forEach(g => { byId[g.querySelector("a.doc").dataset.doc] = g; });

function setActive(id){
  groups.forEach(g => {
    const on = g === byId[id];
    g.classList.toggle("open", on);
    g.querySelector("a.doc").classList.toggle("on", on);
  });
}
setActive(DOCS[0].id);

const io = new IntersectionObserver(entries => {
  const vis = entries.filter(e => e.isIntersecting)
                     .sort((a,b) => a.boundingClientRect.top - b.boundingClientRect.top);
  if (vis.length) setActive(vis[0].target.id);
}, { rootMargin:"-15% 0px -70% 0px" });
DOCS.forEach(d => { const el = document.getElementById(d.id); if (el) io.observe(el); });
</script>
"""

if __name__ == "__main__":
    sys.exit(main())
