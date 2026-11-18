#!/usr/bin/env python3
"""Summarise the near-complete re-run, which was made against the module as it
ships (raised bit cutoffs). Reads cov/nc_work."""
import collections, glob, json, os
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
meta = {l.split('\t')[0]: (l.split('\t')[1], int(l.rstrip('\n').split('\t')[6] or 0))
        for l in open('tob.full.tsv')}
rows = sorted(glob.glob('cov/nc_work/Quality_GTOs/*.qual.gto'))
qual, routed, feats, flags = {}, {}, {}, {}
for q in rows:
    g = os.path.basename(q).split('.qual')[0]
    d = json.load(open(q))
    mods, ks = set(), set()
    for f in d.get('features') or []:
        fa = (f.get('family_assignments') or [[None, '']])[0]
        k = fa[1] if len(fa) > 1 else ''
        if k and '.' in k:
            mod, _, rest = k.partition('.')
            head, dot, tail = rest.rpartition('.')
            mods.add(mod); ks.add(head if dot and tail.replace('_', '').isdigit() else rest)
        elif (f.get('function') or '').startswith('ORF1ab'): ks.add('ORF1AB')
    routed[g] = sorted(mods); qual[g] = (d.get('quality') or {}).get('genome_quality')
    feats[g] = ks
    fl = list(d.get('genome_quality_flags') or [])
    for f in d.get('features') or []:
        if f.get('feature_quality_flags'):
            fl.append('%s: %s' % ((f.get('function') or '?')[:26], '; '.join(f['feature_quality_flags'])))
    flags[g] = fl
n = len(qual)
good = sum(1 for g in qual if qual[g] == 'Good')
r = sum(1 for g in qual if routed[g])
print('near-complete exemplars re-run against the shipped module: %d' % n)
print('  routed %d (%.1f%%), good %d (%.1f%%), poor %d'
      % (r, 100.0 * r / max(1, n), good, 100.0 * good / max(1, n), n - good))
bym = collections.defaultdict(list)
for g in qual:
    for m in routed[g] or ['UNROUTED']: bym[m].append(g)
for m, gs in sorted(bym.items(), key=lambda x: -len(x[1])):
    c = collections.Counter()
    for g in gs:
        for k in feats[g]: c[k] += 1
    print('  %-16s %d genome(s): %s' % (m, len(gs),
          ', '.join('%s %d%%' % (k, round(100.0 * v / len(gs))) for k, v in c.most_common())))
poor = [g for g in qual if qual[g] != 'Good']
if poor:
    print('\n  poor, with the reason:')
    for g in sorted(poor):
        print('    %-12s %-40s %s' % (g, meta.get(g, ('', 0))[0][:40],
                                      ' | '.join(flags[g])[:100] or '(no flag)'))
