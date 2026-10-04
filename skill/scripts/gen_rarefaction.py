#!/usr/bin/env python3
"""Render the vocabulary rarefaction curve as a standalone Artifact page.

Curves are over PROTEINS sampled and start at 0 (see annotation_rarefaction.py).
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
if d.get("unit") != "proteins":
    raise SystemExit("%s is a genome-level rarefaction from the old annotation_rarefaction.py; "
                     "re-run it -- the curve must be over proteins and start at 0" % _a.rarefaction)
S, M = d['source'], d['module']
G = d['n_genomes']
reps = d['replicates']
XMAX = max(S['n'], M['n'])

# ---- geometry -----------------------------------------------------------
VBW, VBH = 780, 400
L, R, T, B = 58, 132, 26, 52          # generous right margin for direct labels
PW, PH = VBW - L - R, VBH - T - B
#  Axes, tick values and every rate on this page are derived from the input.
def _nice(hi):
    """A round axis maximum at or above hi, with a sensible tick step."""
    import math
    if hi <= 0:
        return 1, [0, 1]
    mag = 10 ** math.floor(math.log10(hi))
    for m in (1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10):
        top = m * mag
        if top >= hi * 1.06:
            break
    step = top / 4.0
    ticks = [int(round(step * k)) if step >= 1 else round(step * k, 2)
             for k in range(5)]
    return top, ticks

YMAX, yticks = _nice(max(max(S['y']), max(M['y'])))
XTOP, xticks = _nice(XMAX)
xs = lambda x: L + PW * x / XTOP
ys = lambda v: T + PH * (1 - v / YMAX)

def path(c):
    return ' '.join(('M' if i == 0 else 'L') + f'{xs(x):.1f} {ys(v):.2f}'
                    for i, (x, v) in enumerate(zip(c['x'], c['y'])))

grid = ''.join(
    f'<line x1="{L}" y1="{ys(v):.1f}" x2="{L+PW}" y2="{ys(v):.1f}" '
    f'stroke="var(--grid)" stroke-width="1"/>' for v in yticks)
ylab = ''.join(
    f'<text x="{L-11}" y="{ys(v)+4:.1f}" text-anchor="end" class="tick">{v}</text>'
    for v in yticks)
xlab = ''.join(
    f'<text x="{xs(t):.1f}" y="{T+PH+21}" text-anchor="middle" class="tick">{t:,}</text>'
    for t in xticks)

knee = d['source_knee']
srcs, mods = d['source_total'], d['module_total']
CN = d['common_n']
ratio = d['source_at_common'] / d['module_at_common'] if d['module_at_common'] else float('nan')
src_end, mod_end = S['y'][-1], M['y'][-1]
#  one marker where both curves have been sampled to the same depth
cx = xs(CN)
#  label positions: keep the two end labels from overlapping
sy, my = ys(src_end), ys(mod_end)
if abs(sy - my) < 30:
    if sy <= my: my = sy + 30
    else: sy = my + 30

def at(c, x):
    """value of a thinned curve at protein count x (linear between kept points)"""
    import bisect
    if x > c['n']: return None
    i = bisect.bisect_left(c['x'], x)
    if c['x'][i] == x: return c['y'][i]
    x0, x1, y0, y1 = c['x'][i-1], c['x'][i], c['y'][i-1], c['y'][i]
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
rows = sorted({int(round(CN * f)) for f in (0.002, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0)} - {0})
table = ''.join(f"<tr><td>{k:,}</td><td>{at(S,k):.2f}</td><td>{at(M,k):.2f}</td>"
                f"<td>{at(S,k)/at(M,k):.2f}&times;</td></tr>" for k in rows)

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
  <h1>{_a.taxon}</h1>
  <h2 style="margin:8px 0 0;font-weight:500;font-size:20px;color:var(--ink2)">A controlled vocabulary stops growing. Free text does not.</h2>
  <p class="stand">Take every annotated protein from {G:,} {_a.taxon} genomes, shuffle
  them, and draw them one at a time, counting the distinct annotation strings seen so far.
  A curated vocabulary needs exactly {mods} names however many proteins you read, so its
  curve flattens; free-text product strings have no ceiling, because each submitter spells
  the same proteins a new way. Over {reps} random orderings the module reaches 95% of its
  vocabulary after <b>{d['module_knee']:,} proteins</b>; the source after
  <b>{knee:,}</b> of its {S['n']:,}.</p>
 </header>

 <div class="figs">
  <div class="fig"><span class="n" style="color:var(--mod)">{mods}</span><span class="k">strings, controlled vocabulary</span></div>
  <div class="fig"><span class="n" style="color:var(--src)">{srcs}</span><span class="k">strings, BV-BRC free text</span></div>
  <div class="fig"><span class="n">{ratio:.1f}&times;</span><span class="k">collapse at {CN:,} proteins sampled</span></div>
  <div class="fig"><span class="n">{knee:,}</span><span class="k">proteins before free text reaches 95%</span></div>
 </div>

 <figure>
  <div class="legend">
   <span><i class="sw" style="background:var(--src)"></i> BV-BRC product strings (free text)</span>
   <span><i class="sw" style="background:var(--mod)"></i> LowVan annotation strings (controlled)</span>
  </div>
  <div class="chartbox" id="cb">
   <svg viewBox="0 0 {VBW} {VBH}" role="img" aria-label="Rarefaction curves: distinct annotation strings against number of proteins sampled. Both start at zero. Free-text strings rise to {srcs} over {S['n']} proteins; the controlled vocabulary flattens at {mods}.">
    {grid}
    <line x1="{L}" y1="{T+PH}" x2="{L+PW}" y2="{T+PH}" stroke="var(--base)" stroke-width="1"/>
    {ylab}{xlab}
    <text class="axtitle" x="{L}" y="{VBH-8}">proteins sampled</text>
    <text class="axtitle" x="{L-11}" y="{T-10}" text-anchor="end">strings</text>

    <line x1="{cx:.1f}" y1="{T}" x2="{cx:.1f}" y2="{T+PH}" stroke="var(--base)" stroke-width="1" stroke-dasharray="3 3"/>
    <path d="{path(S)}" fill="none" stroke="var(--src)" stroke-width="2" stroke-linejoin="round"/>
    <path d="{path(M)}" fill="none" stroke="var(--mod)" stroke-width="2" stroke-linejoin="round"/>

    <circle cx="{xs(S['n']):.1f}" cy="{ys(src_end):.2f}" r="4.5" fill="var(--src)" stroke="var(--surface)" stroke-width="2"/>
    <circle cx="{xs(M['n']):.1f}" cy="{ys(mod_end):.2f}" r="4.5" fill="var(--mod)" stroke="var(--surface)" stroke-width="2"/>
    <text class="endlab" x="{min(xs(S['n'])+10, L+PW+8):.1f}" y="{sy-1:.2f}" fill="var(--src)">free text</text>
    <text class="endsub" x="{min(xs(S['n'])+10, L+PW+8):.1f}" y="{sy+12:.2f}">{srcs} strings</text>
    <text class="endlab" x="{min(xs(M['n'])+10, L+PW+8):.1f}" y="{my+3:.2f}" fill="var(--mod)">controlled</text>
    <text class="endsub" x="{min(xs(M['n'])+10, L+PW+8):.1f}" y="{my+16:.2f}">{mods} strings</text>

    <line id="xh" x1="0" y1="{T}" x2="0" y2="{T+PH}" stroke="var(--base)" stroke-width="1" opacity="0"/>
    <circle id="ds" r="4" fill="var(--src)" stroke="var(--surface)" stroke-width="2" opacity="0"/>
    <circle id="dm" r="4" fill="var(--mod)" stroke="var(--surface)" stroke-width="2" opacity="0"/>
    <rect id="hit" x="{L}" y="{T}" width="{PW}" height="{PH}" fill="transparent" style="cursor:crosshair"/>
   </svg>
   <div class="tip" id="tip"></div>
  </div>
  <figcaption>Mean distinct annotation strings over {reps} random orderings of every
  annotated protein in the {G:,} genomes both sources annotate: {S['n']:,} source
  product records and {M['n']:,} module calls. Both curves start at zero. The dashed line
  marks {CN:,} proteins, the depth both curves reach, where the two are compared.</figcaption>
 </figure>

 <p class="note">Compared at the same sampling depth, {CN:,} proteins, free text has
 used <b>{d['source_at_common']:.1f} distinct strings</b> and the controlled vocabulary
 <b>{d['module_at_common']:.1f}</b> &mdash; a <b>{ratio:.1f}-fold collapse</b>. The
 source vocabulary reaches 95% of its final size after {knee:,} of its {S['n']:,} proteins;
 the module after {d['module_knee']:,} of {M['n']:,}.</p>

 <details>
  <summary>Table view &mdash; proteins sampled</summary>
  <div class="tw"><table>
   <thead><tr><th>proteins sampled</th><th>free text</th><th>controlled</th><th>ratio</th></tr></thead>
   <tbody>
   {table}
   </tbody>
  </table></div>
 </details>

 <div class="src">Generated by <code>skill/scripts/gen_rarefaction.py</code> from
 <code>rarefaction.json</code>, written by the kit&rsquo;s
 <code>annotation_rarefaction.py</code> at {reps} replicates. Categorical pair
 <code>#2a78d6</code>/<code>#eb6834</code> (light) and <code>#3987e5</code>/<code>#d95926</code>
 (dark), validated by <code>validate_palette.py</code>: worst CVD &Delta;E 24.7,
 normal-vision &Delta;E 31.8, contrast &ge;3.12:1.</div>
</div>
<script>
 var SX={json.dumps(S['x'])}, SY={json.dumps(S['y'])}, MX={json.dumps(M['x'])}, MY={json.dumps(M['y'])};
 var Lm={L}, PWm={PW}, XT={XTOP}, VB={VBW};
 var cb=document.getElementById('cb'), hit=document.getElementById('hit'),
     xh=document.getElementById('xh'), ds=document.getElementById('ds'),
     dm=document.getElementById('dm'), tip=document.getElementById('tip');
 var svg=cb.querySelector('svg');
 function X(x){{ return Lm + PWm*x/XT; }}
 function at(xa,ya,x){{ if(x>xa[xa.length-1]) return null; var i=0; while(xa[i]<x) i++;
   if(xa[i]===x||i===0) return ya[i]; return ya[i-1]+(ya[i]-ya[i-1])*(x-xa[i-1])/(xa[i]-xa[i-1]); }}
 function Y(v){{ return {T} + {PH}*(1 - v/{YMAX}); }}
 function show(e){{
   var r=svg.getBoundingClientRect(), sx=(e.clientX-r.left)*VB/r.width;
   var k=Math.round((sx-Lm)/PWm*XT); var top=Math.max(SX[SX.length-1],MX[MX.length-1]);
   if(k<0) k=0; if(k>top) k=top;
   var x=X(k), s=at(SX,SY,k), m=at(MX,MY,k);
   xh.setAttribute('x1',x); xh.setAttribute('x2',x); xh.setAttribute('opacity','1');
   if(s!==null){{ ds.setAttribute('cx',x); ds.setAttribute('cy',Y(s)); ds.setAttribute('opacity','1'); }} else ds.setAttribute('opacity','0');
   if(m!==null){{ dm.setAttribute('cx',x); dm.setAttribute('cy',Y(m)); dm.setAttribute('opacity','1'); }} else dm.setAttribute('opacity','0');
   tip.innerHTML='<div style="color:var(--muted);margin-bottom:3px">'+k.toLocaleString()+
     ' protein'+(k===1?'':'s')+' sampled</div>'+
     '<div class="tiprow"><i class="sw" style="background:var(--src)"></i>free text <b>'+(s===null?'&mdash;':s.toFixed(2))+'</b></div>'+
     '<div class="tiprow"><i class="sw" style="background:var(--mod)"></i>controlled <b>'+(m===null?'&mdash;':m.toFixed(2))+'</b></div>';
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
print(f"  source 0 -> {src_end:.2f} over {S['n']} proteins ({srcs} total), knee {knee}")
print(f"  module 0 -> {mod_end:.2f} over {M['n']} proteins ({mods} total), knee {d['module_knee']}")
print(f"  at {CN} proteins: {d['source_at_common']:.1f} vs {d['module_at_common']:.1f} = {ratio:.1f}x")
