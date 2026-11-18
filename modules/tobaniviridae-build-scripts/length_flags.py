#!/usr/bin/env python3
"""For every 'Feature is too short/long' flag in the whole-taxon run: what did
the module declare, and what was actually called? Step 11's question 3 -- when
that flag dominates, the bounds are usually wrong rather than the genomes."""
import collections, glob, json, os
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
MODS = ['Torovirinae', 'Piscanivirinae', 'Remotovirinae', 'Serpentovirinae']
decl = {}
for M in MODS:
    for k, v in json.load(open('%s/%s_Viral_PSSM.json' % (M, M)))[M]['features'].items():
        decl[(M, v.get('anno'))] = (v.get('min_len'), v.get('max_len'), k)
rows = collections.defaultdict(list)
for q in sorted(glob.glob('cov/work/Quality_GTOs/*.qual.gto')):
    g = os.path.basename(q).split('.qual')[0]
    d = json.load(open(q))
    mod = None
    for f in d.get('features') or []:
        fa = (f.get('family_assignments') or [[None, '']])[0]
        k = fa[1] if len(fa) > 1 else ''
        if k and '.' in k: mod = k.split('.')[0]
    for f in d.get('features') or []:
        for fl in f.get('feature_quality_flags') or []:
            if 'too short' in fl or 'too long' in fl:
                aa = len(f.get('protein_translation') or f.get('aa_sequence') or '')
                lo, hi, key = decl.get((mod, f.get('function')), (None, None, '?'))
                rows[(mod, key)].append((g, aa, lo, hi, fl))
for (mod, key), v in sorted(rows.items(), key=lambda x: -len(x[1])):
    lens = sorted(x[1] for x in v)
    lo, hi = v[0][2], v[0][3]
    print('%-16s %-8s declared %s-%s aa; %d flagged call(s), lengths %d-%d'
          % (mod, key, lo, hi, len(v), lens[0], lens[-1]))
    for g, aa, _lo, _hi, fl in v[:4]:
        print('      %-12s called %4d aa  %s' % (g, aa, fl))
