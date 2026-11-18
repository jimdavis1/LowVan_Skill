#!/usr/bin/env python3
"""Routing and calls for the exemplars the GTO pipeline could not annotate.

Those 26 failed at rast-create-genome with HTTP 403 from the BV-BRC ID server,
so they were never offered to a module. The coverage report counts them as
UNROUTED, which reads as a routing failure; this measures what actually
happens to them with the flat annotator, which needs no ID server."""
import collections, glob, os
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
meta = {l.split('\t')[0]: (l.split('\t')[1], int(l.rstrip('\n').split('\t')[6] or 0))
        for l in open('tob.full.tsv')}
staged = [os.path.basename(p).rsplit('.fna', 1)[0] for p in glob.glob('cov/work/ex_retry3/*')]
mods, feats = {}, {}
for g in staged:
    t = 'cov/flat/%s.feature.tbl' % g
    if not os.path.exists(t): continue
    rows = [l.rstrip('\n').split('\t') for l in open(t) if l.count('\t') > 6]
    mods[g] = sorted({r[11] for r in rows if len(r) > 11})
    feats[g] = [r[6] for r in rows]
print('%d of %d ID-server-failed exemplars annotated with the flat annotator'
      % (len(mods), len(staged)))
routed = [g for g in mods if mods[g]]
print('routed: %d/%d' % (len(routed), len(mods)))
bym = collections.Counter(m for g in routed for m in mods[g])
print('modules they route to:', dict(bym))
print()
for g in sorted(mods):
    print('%-12s %-40s %6d nt  %-16s %s'
          % (g, meta.get(g, ('', 0))[0][:40], meta.get(g, ('', 0))[1],
             ','.join(mods[g]) or 'UNROUTED', ' '.join(feats[g])))
