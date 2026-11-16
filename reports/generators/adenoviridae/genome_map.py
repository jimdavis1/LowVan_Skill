#!/usr/bin/env python3
"""To-scale genome map of one reference genome as the module annotates it, as an
inline SVG <figure> for the coverage audit (section 01).
  + strand above the axis, - strand below; genes as arrows labelled by gene symbol
  essential (copy_num) = --mod, other called CDS = --muted, mature peptides = --warn
  spliced features: exons joined by a thin line
usage: genome_map.py <Module> <genome_id> <out.html>"""
import json, sys, os, html

m, gid, out = sys.argv[1:4]
W = os.path.expanduser('~/lowvan_scratch/adeno_v3')
G = json.load(open(f'{W}/panel/{m}.rules.anno/{gid}.qual.gto'))
F = json.load(open(f'{W}/build/{m}/{m}_Viral_PSSM.json'))[m]['features']
L = sum(len(c['dna']) for c in G['contigs'])

X0, XW = 64, 900                      # drawing span (px)
sx = lambda nt: X0 + XW * (nt - 1) / L
CH = 6.2                              # approx px per label character at 10.5px
BOX, LANE = 10, 25                    # arrow height, lane pitch

feats = []
for f in G['features']:
    fa = f.get('family_assignments')
    if not fa:
        continue
    key = fa[0][1].split('.')[1]
    spec = F.get(key, {})
    exons = []
    for c, s, st, n in f['location']:
        s, n = int(s), int(n)
        a, b = (s, s + n - 1) if st == '+' else (s - n + 1, s)
        exons.append((a, b))
    strand = f['location'][0][2]
    kind = 'mat' if f['type'] == 'mat_peptide' else ('ess' if spec.get('copy_num') else 'oth')
    sym = spec.get('gene_symbol') or key
    feats.append(dict(key=key, sym=sym, exons=sorted(exons), strand=strand, kind=kind,
                      lo=min(a for a, _ in exons), hi=max(b for _, b in exons)))

# lanes per strand: occupancy = box extent or label extent, whichever is wider
lanes = {'+': [], '-': []}
for f in sorted(feats, key=lambda f: (f['kind'] == 'mat', f['lo'])):
    x1, x2 = sx(f['lo']), sx(f['hi'])
    cx = (x1 + x2) / 2
    lw = len(f['sym']) * CH
    occ = (min(x1, cx - lw / 2) - 3, max(x2, cx + lw / 2) + 3)
    ls = lanes[f['strand']]
    for i, used in enumerate(ls):
        if all(occ[1] < u0 or occ[0] > u1 for u0, u1 in used):
            used.append(occ); f['lane'] = i; break
    else:
        ls.append([occ]); f['lane'] = len(ls) - 1

np_, nm_ = len(lanes['+']), len(lanes['-'])
TOP = 58                               # legend band
AX = TOP + np_ * LANE + 18             # axis y
H = AX + 34 + nm_ * LANE + 14
COL = {'ess': 'var(--mod)', 'oth': 'var(--muted)', 'mat': 'var(--warn)'}

def arrow(a, b, y, strand, col, tip):
    x1, x2 = sx(a), sx(b) + 1
    h = BOX if col != COL['mat'] else BOX - 4
    y0 = y - h / 2
    t = min(5, (x2 - x1) * 0.6) if tip else 0
    if strand == '+':
        pts = [(x1, y0), (x2 - t, y0), (x2, y), (x2 - t, y0 + h), (x1, y0 + h)]
    else:
        pts = [(x1 + t, y0), (x2, y0), (x2, y0 + h), (x1 + t, y0 + h), (x1, y)]
    return '<polygon points="%s" style="fill:%s"/>' % (' '.join('%.1f,%.1f' % p for p in pts), col)

svg = []
for f in feats:
    y = (AX - 18 - f['lane'] * LANE) if f['strand'] == '+' else (AX + 30 + f['lane'] * LANE)
    col = COL[f['kind']]
    ex = f['exons']
    if len(ex) > 1:
        svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" style="stroke:%s" stroke-width="1"/>'
                   % (sx(ex[0][1]), y, sx(ex[-1][0]), y, col))
    last = ex[-1] if f['strand'] == '+' else ex[0]
    for a, b in ex:
        svg.append(arrow(a, b, y, f['strand'], col, (a, b) == last))
    cx = (sx(f['lo']) + sx(f['hi'])) / 2
    ty = (y - BOX / 2 - 3) if f['strand'] == '+' else (y + BOX / 2 + 11)
    svg.append('<text x="%.1f" y="%.1f" text-anchor="middle" font-size="10.5" fill="currentColor">%s</text>'
               % (cx, ty, html.escape(f['sym'])))

ax = ['<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="currentColor" stroke-width="1.5"/>' % (X0, AX, X0 + XW, AX)]
for k in range(0, L + 1, 5000):
    x = sx(max(k, 1))
    ax.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="currentColor"/>' % (x, AX - 4, x, AX + 4))
    ax.append('<text x="%.1f" y="%d" text-anchor="middle" font-size="10.5" fill="currentColor">%d kb</text>' % (x, AX + 17, k // 1000))
ax.append('<text x="4" y="%d" font-size="10.5" fill="currentColor">+ strand</text>' % (AX - 14))
ax.append('<text x="4" y="%d" font-size="10.5" fill="currentColor">&#8722; strand</text>' % (AX + 34))

leg, lx = [], 4
for kind, label in (('ess', 'essential CDS'), ('oth', 'other called CDS'), ('mat', 'mature peptide')):
    leg.append('<rect x="%d" y="14" width="18" height="10" style="fill:%s"/>' % (lx, COL[kind]))
    leg.append('<text x="%d" y="23" font-size="11" fill="currentColor">%s</text>' % (lx + 24, label))
    lx += 24 + len(label) * 6.6 + 26
leg.append('<line x1="%d" y1="19" x2="%d" y2="19" stroke="currentColor"/>' % (lx, lx + 18))
leg.append('<text x="%d" y="23" font-size="11" fill="currentColor">intron (spliced)</text>' % (lx + 24))

name = G.get('scientific_name', gid)
ne = sum(1 for f in feats if F.get(f['key'], {}).get('copy_num')); nm = sum(f['kind'] == 'mat' for f in feats)
claim = ('%s (%s, %.1f kb) as the %s module annotates it: %d features (%d essential, %d mature peptides).'
         % (name, gid, L / 1000, m, len(feats), ne, nm))
page = ('<figure><svg viewBox="0 0 980 %d" role="img" aria-label="%s">%s%s%s</svg>'
        '<figcaption>%s Drawn to scale from this genome&rsquo;s own calls; forward-strand genes above the axis, '
        'reverse-strand below.</figcaption></figure>'
        % (H, html.escape(claim), ''.join(leg), ''.join(ax), ''.join(svg), html.escape(claim)))
open(out, 'w').write('<!doctype html><meta charset="utf-8">' + page)
print(out, len(feats), 'features', np_, '+lanes', nm_, '-lanes', 'height', H)
