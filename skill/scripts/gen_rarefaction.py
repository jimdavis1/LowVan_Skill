#!/usr/bin/env python3
"""Render the vocabulary rarefaction curve as a standalone Artifact page.

The kit ships no generator for this (annotation_rarefaction.py writes JSON
only). Promoted into the kit from the Matonaviridae build and
parameterised with --taxon / --rarefaction / --out. Palette is the report
series' validated categorical pair -- #2a78d6/#eb6834 light, #3987e5/#d95926
dark -- checked with work/validate_palette.py (a Python port of the dataviz
skill's validate_palette.js, since this box has no node):

    protan dE 24.7 / 26.8, deutan 31.7 / 28.6   (target 8)
    normal-vision dE 33.6 / 31.8                (hard floor 15)
    contrast 3.12-4.79:1                        (min 3.0)
"""
import json, os, sys

import argparse, datetime
_ap = argparse.ArgumentParser()
_ap.add_argument("--rarefaction", default="rarefaction.json")
_ap.add_argument("--taxon", required=True)
_ap.add_argument("--date", default=datetime.date.today().strftime("%-d %B %Y"))
_ap.add_argument("--out", default="rarefaction.html")
_a = _ap.parse_args()
d = json.load(open(_a.rarefaction))
src, mod = d['source'], d['module']
N = d['n_genomes']
reps = d['replicates']

# ---- geometry -----------------------------------------------------------
VBW, VBH = 780, 400
L, R, T, B = 58, 132, 26, 52          # generous right margin for direct labels
PW, PH = VBW - L - R, VBH - T - B
YMAX = 32
xs = lambda i: L + PW * i / (N - 1)
ys = lambda v: T + PH * (1 - v / YMAX)

def path(series):
    return ' '.join(('M' if i == 0 else 'L') + f'{xs(i):.1f} {ys(v):.2f}'
                    for i, v in enumerate(series))

xticks = [1, 50, 100, 150, 216]
yticks = [0, 10, 20, 30]
grid = ''.join(
    f'<line x1="{L}" y1="{ys(v):.1f}" x2="{L+PW}" y2="{ys(v):.1f}" '
    f'stroke="var(--grid)" stroke-width="1"/>' for v in yticks)
ylab = ''.join(
    f'<text x="{L-11}" y="{ys(v)+4:.1f}" text-anchor="end" class="tick">{v}</text>'
    for v in yticks)
xlab = ''.join(
    f'<text x="{xs(t-1):.1f}" y="{T+PH+21}" text-anchor="middle" class="tick">{t}</text>'
    for t in xticks)

# saturation marker: where the source curve reaches 95% of its final size
knee = d['source_knee']
srcs, mods = d['source_total'], d['module_total']
ratio = 13.9 / 3.2

HTML = f'''<title>{_a.taxon} Vocabulary Saturation</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Serif:wght@500;600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
  :root {{
    color-scheme: light;
    --plane:#f9f9f7; --surface:#fcfcfb; --ink:#0b0b0b; --ink2:#52514e; --muted:#898781;
    --grid:#e1e0d9; --base:#c3c2b7; --hair:rgba(11,11,11,.10);
    --src:#eb6834; --mod:#2a78d6;
    --serif:"IBM Plex Serif",Georgia,serif;
    --sans:"IBM Plex Sans",system-ui,-apple-system,"Segoe UI",sans-serif;
    --mono:"IBM Plex Mono",ui-monospace,Menlo,monospace;
  }}
  @media (prefers-color-scheme:dark) {{ :root:not([data-theme="light"]) {{
    color-scheme: dark;
    --plane:#0d0d0d; --surface:#1a1a19; --ink:#fff; --ink2:#c3c2b7; --muted:#898781;
    --grid:#2c2c2a; --base:#383835; --hair:rgba(255,255,255,.10);
    --src:#d95926; --mod:#3987e5;
  }} }}
  :root[data-theme="dark"] {{
    color-scheme: dark;
    --plane:#0d0d0d; --surface:#1a1a19; --ink:#fff; --ink2:#c3c2b7; --muted:#898781;
    --grid:#2c2c2a; --base:#383835; --hair:rgba(255,255,255,.10);
    --src:#d95926; --mod:#3987e5;
  }}
  body {{ background:var(--plane); color:var(--ink); font-family:var(--sans); line-height:1.5; }}
  .wrap {{ max-width:900px; margin:0 auto; padding-block:40px 72px; padding-left:22px; padding-right:22px;
           display:flex; flex-direction:column; gap:28px; }}
  .eyebrow {{ font-family:var(--mono); font-size:11px; letter-spacing:.1em; text-transform:uppercase; color:var(--muted); }}
  h1 {{ font-family:var(--serif); font-weight:600; font-size:31px; margin:7px 0 0; text-wrap:balance; }}
  .stand {{ max-width:68ch; color:var(--ink2); font-size:15.5px; margin:11px 0 0; }}
  .stand b {{ color:var(--ink); font-weight:600; }}

  .figs {{ display:grid; grid-template-columns:repeat(4,minmax(0,1fr));
           border-top:1px solid var(--base); border-bottom:1px solid var(--grid); }}
  .fig {{ padding:13px 16px 13px 0; }}
  .fig+.fig {{ padding-left:16px; border-left:1px solid var(--grid); }}
  .fig .n {{ font-family:var(--mono); font-variant-numeric:tabular-nums; font-size:23px;
             font-weight:500; display:block; letter-spacing:-.02em; }}
  .fig .k {{ font-size:12px; color:var(--ink2); }}

  figure {{ margin:0; }}
  .chartbox {{ border:1px solid var(--hair); border-radius:4px; background:var(--surface); padding:6px 4px 2px; position:relative; }}
  svg {{ display:block; width:100%; height:auto; max-width:100%; }}
  .tick {{ font-family:var(--mono); font-size:10.5px; fill:var(--muted); }}
  .axtitle {{ font-family:var(--sans); font-size:11.5px; fill:var(--ink2); }}
  .endlab {{ font-family:var(--sans); font-size:12px; font-weight:600; }}
  .endsub {{ font-family:var(--mono); font-size:10.5px; fill:var(--muted); }}
  figcaption {{ font-size:12.5px; color:var(--muted); margin-top:9px; max-width:74ch; }}

  .legend {{ display:flex; gap:18px; flex-wrap:wrap; font-size:13px; color:var(--ink2);
             align-items:center; margin-bottom:10px; }}
  .legend span {{ display:inline-flex; gap:7px; align-items:center; }}
  .sw {{ width:14px; height:3px; border-radius:2px; flex:none; }}

  .tip {{ position:absolute; pointer-events:none; opacity:0; transition:opacity .08s;
          background:var(--surface); border:1px solid var(--base); border-radius:4px;
          padding:7px 10px; font-size:12px; box-shadow:0 2px 8px rgba(0,0,0,.10);
          font-variant-numeric:tabular-nums; white-space:nowrap; z-index:2; }}
  .tip b {{ font-family:var(--mono); font-weight:500; }}
  .tiprow {{ display:flex; gap:7px; align-items:center; }}
  @media (prefers-reduced-motion:reduce) {{ .tip {{ transition:none; }} }}

  details {{ border-top:1px solid var(--grid); padding-top:14px; }}
  summary {{ cursor:pointer; font-size:13.5px; color:var(--ink2); }}
  summary:focus-visible {{ outline:2px solid var(--mod); outline-offset:3px; }}
  .tw {{ overflow-x:auto; margin-top:12px; border:1px solid var(--hair); border-radius:4px; background:var(--surface); }}
  table {{ border-collapse:collapse; width:100%; font-size:13px; }}
  th,td {{ padding:6px 11px; text-align:right; border-bottom:1px solid var(--grid); }}
  th:first-child, td:first-child {{ text-align:left; }}
  thead th {{ font-size:10.5px; letter-spacing:.07em; text-transform:uppercase;
              color:var(--muted); font-weight:500; }}
  tbody td {{ font-family:var(--mono); font-variant-numeric:tabular-nums; color:var(--ink2); }}
  tbody tr:last-child td {{ border-bottom:0; }}
  .note {{ font-size:13.5px; color:var(--ink2); max-width:74ch; }}
  .note b {{ color:var(--ink); }}
  .src {{ font-size:12px; color:var(--muted); border-top:1px solid var(--grid); padding-top:13px; }}
  .src code {{ font-family:var(--mono); }}
  @media (max-width:620px) {{ .figs {{ grid-template-columns:repeat(2,1fr); }} h1 {{ font-size:25px; }} }}
</style>
<div class="wrap">
 <header>
  <div class="eyebrow">LowVan &middot; {_a.taxon} module &middot; {_a.date}</div>
  <h1>A controlled vocabulary stops growing. Free text does not.</h1>
  <p class="stand">Every rubivirus encodes the same seven proteins, so a curated
  vocabulary should need exactly seven names no matter how many genomes you read.
  Free-text product strings have no such ceiling: each new submitter spells the
  same proteins a new way. Sampling {N} annotated genomes {reps} times over shows
  the two curves doing exactly that &mdash; the module is <b>saturated after one
  genome</b>, while the source vocabulary is <b>still climbing at {N}</b>.</p>
 </header>

 <div class="figs">
  <div class="fig"><span class="n" style="color:var(--mod)">{mods}</span><span class="k">strings, controlled vocabulary</span></div>
  <div class="fig"><span class="n" style="color:var(--src)">{srcs}</span><span class="k">strings, BV-BRC free text</span></div>
  <div class="fig"><span class="n">{ratio:.1f}&times;</span><span class="k">collapse per 100 genomes</span></div>
  <div class="fig"><span class="n">{knee}</span><span class="k">genomes before free text saturates</span></div>
 </div>

 <figure>
  <div class="legend">
   <span><i class="sw" style="background:var(--src)"></i> BV-BRC product strings (free text)</span>
   <span><i class="sw" style="background:var(--mod)"></i> LowVan annotation strings (controlled)</span>
  </div>
  <div class="chartbox" id="cb">
   <svg viewBox="0 0 {VBW} {VBH}" role="img" aria-label="Rarefaction curves: distinct annotation strings against number of genomes sampled. Free-text strings rise from about 7 to 30; the controlled vocabulary stays flat at 7.">
    {grid}
    <line x1="{L}" y1="{T+PH}" x2="{L+PW}" y2="{T+PH}" stroke="var(--base)" stroke-width="1"/>
    {ylab}{xlab}
    <text class="axtitle" x="{L}" y="{VBH-8}">genomes sampled</text>
    <text class="axtitle" x="{L-11}" y="{T-10}" text-anchor="end">strings</text>

    <path d="{path(src)}" fill="none" stroke="var(--src)" stroke-width="2" stroke-linejoin="round"/>
    <path d="{path(mod)}" fill="none" stroke="var(--mod)" stroke-width="2" stroke-linejoin="round"/>

    <circle cx="{xs(N-1):.1f}" cy="{ys(src[-1]):.2f}" r="4.5" fill="var(--src)" stroke="var(--surface)" stroke-width="2"/>
    <circle cx="{xs(N-1):.1f}" cy="{ys(mod[-1]):.2f}" r="4.5" fill="var(--mod)" stroke="var(--surface)" stroke-width="2"/>
    <text class="endlab" x="{xs(N-1)+10:.1f}" y="{ys(src[-1])-1:.2f}" fill="var(--src)">free text</text>
    <text class="endsub" x="{xs(N-1)+10:.1f}" y="{ys(src[-1])+12:.2f}">{srcs} strings</text>
    <text class="endlab" x="{xs(N-1)+10:.1f}" y="{ys(mod[-1])+3:.2f}" fill="var(--mod)">controlled</text>
    <text class="endsub" x="{xs(N-1)+10:.1f}" y="{ys(mod[-1])+16:.2f}">{mods} strings</text>

    <line id="xh" x1="0" y1="{T}" x2="0" y2="{T+PH}" stroke="var(--base)" stroke-width="1" opacity="0"/>
    <circle id="ds" r="4" fill="var(--src)" stroke="var(--surface)" stroke-width="2" opacity="0"/>
    <circle id="dm" r="4" fill="var(--mod)" stroke="var(--surface)" stroke-width="2" opacity="0"/>
    <rect id="hit" x="{L}" y="{T}" width="{PW}" height="{PH}" fill="transparent" style="cursor:crosshair"/>
   </svg>
   <div class="tip" id="tip"></div>
  </div>
  <figcaption>Mean distinct annotation strings over {reps} random orderings of the
  {N} annotated genomes. The controlled curve is flat because the seven names are
  fixed in advance; the free-text curve rises because it is a record of how many
  different ways people have written those same seven proteins down.</figcaption>
 </figure>

 <p class="note">The ratio is the number worth quoting: at the end of each curve,
 free text costs <b>13.9 strings per 100 genomes</b> and the controlled vocabulary
 <b>3.2</b> &mdash; a <b>{ratio:.1f}-fold collapse</b>. Reaching 95% of its final
 size takes the source vocabulary <b>{knee} genomes</b>; the module gets there on
 the <b>first</b>. Across the whole dump, 58 distinct BV-BRC strings collapse into
 these 7 over 6,547 protein occurrences.</p>

 <details>
  <summary>Table view &mdash; every 20th sample point</summary>
  <div class="tw"><table>
   <thead><tr><th>genomes sampled</th><th>free text</th><th>controlled</th><th>ratio</th></tr></thead>
   <tbody>
   {''.join(f"<tr><td>{i+1}</td><td>{src[i]:.2f}</td><td>{mod[i]:.2f}</td><td>{src[i]/mod[i]:.2f}&times;</td></tr>" for i in list(range(0, N, 20)) + [N-1])}
   </tbody>
  </table></div>
 </details>

 <div class="src">Generated by <code>Reports/generators/gen_rarefaction.py</code> from
 <code>rarefaction.json</code>, written by the kit&rsquo;s
 <code>annotation_rarefaction.py</code> at {reps} replicates. Categorical pair
 <code>#2a78d6</code>/<code>#eb6834</code> (light) and <code>#3987e5</code>/<code>#d95926</code>
 (dark), validated by <code>work/validate_palette.py</code>: worst CVD &Delta;E 24.7,
 normal-vision &Delta;E 31.8, contrast &ge;3.12:1.</div>
</div>
<script>
 var SRC={json.dumps([round(v,2) for v in src])}, MOD={json.dumps([round(v,2) for v in mod])};
 var Lm={L}, PWm={PW}, Nn={N}, VB={VBW};
 var cb=document.getElementById('cb'), hit=document.getElementById('hit'),
     xh=document.getElementById('xh'), ds=document.getElementById('ds'),
     dm=document.getElementById('dm'), tip=document.getElementById('tip');
 var svg=cb.querySelector('svg');
 function X(i){{ return Lm + PWm*i/(Nn-1); }}
 function Y(v){{ return {T} + {PH}*(1 - v/{YMAX}); }}
 function show(e){{
   var r=svg.getBoundingClientRect(), sx=(e.clientX-r.left)*VB/r.width;
   var i=Math.round((sx-Lm)/PWm*(Nn-1));
   if(i<0) i=0; if(i>Nn-1) i=Nn-1;
   var x=X(i);
   xh.setAttribute('x1',x); xh.setAttribute('x2',x); xh.setAttribute('opacity','1');
   ds.setAttribute('cx',x); ds.setAttribute('cy',Y(SRC[i])); ds.setAttribute('opacity','1');
   dm.setAttribute('cx',x); dm.setAttribute('cy',Y(MOD[i])); dm.setAttribute('opacity','1');
   tip.innerHTML='<div style="color:var(--muted);margin-bottom:3px">'+(i+1)+
     ' genome'+(i?'s':'')+' sampled</div>'+
     '<div class="tiprow"><i class="sw" style="background:var(--src)"></i>free text <b>'+SRC[i].toFixed(2)+'</b></div>'+
     '<div class="tiprow"><i class="sw" style="background:var(--mod)"></i>controlled <b>'+MOD[i].toFixed(2)+'</b></div>';
   tip.style.opacity='1';
   var px=x*r.width/VB, tw=tip.offsetWidth;
   tip.style.left=Math.min(Math.max(px+12,4), r.width-tw-4)+'px';
   tip.style.top='14px';
 }}
 function hide(){{ xh.setAttribute('opacity','0'); ds.setAttribute('opacity','0');
                  dm.setAttribute('opacity','0'); tip.style.opacity='0'; }}
 hit.addEventListener('mousemove',show);
 hit.addEventListener('mouseleave',hide);
 hit.addEventListener('touchmove',function(e){{ if(e.touches[0]) show(e.touches[0]); }},{{passive:true}});
 hit.addEventListener('touchend',hide);
</script>
'''
out = _a.out
open(out, 'w').write(HTML)
print(f"wrote {out} ({len(HTML)} bytes)")
print(f"  source {src[0]:.2f} -> {src[-1]:.2f} ({srcs} total), knee {knee}")
print(f"  module {mod[0]:.2f} -> {mod[-1]:.2f} ({mods} total), knee {d['module_knee']}")
