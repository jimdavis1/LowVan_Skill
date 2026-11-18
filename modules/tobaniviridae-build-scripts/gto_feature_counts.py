#!/usr/bin/env python3
"""Count, straight out of the Quality GTOs, how many genomes carry each feature
by its family_assignments key -- the same source apply_copy_num.py reads, so a
disagreement between this and the feature tables is worth knowing about.
usage: gto_feature_counts.py <Module>"""
import collections, glob, json, os, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
M = sys.argv[1]
n = 0
cnt = collections.Counter()
for q in sorted(glob.glob('gto_out/%s/*.qual.gto' % M)):
    d = json.load(open(q))
    if not d.get('features'): continue
    n += 1
    keys = set()
    for f in d['features']:
        fa = (f.get('family_assignments') or [[None, '']])[0]
        k = fa[1] if len(fa) > 1 else ''
        if k: keys.add(k)
        else: keys.add('(no family_assignment) ' + (f.get('function') or '')[:40])
    for k in keys: cnt[k] += 1
print('%s: %d genome(s)' % (M, n))
for k, v in cnt.most_common(): print('   %-46s %3d  %5.1f%%' % (k, v, 100.0 * v / max(n, 1)))
