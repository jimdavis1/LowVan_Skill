#!/usr/bin/env python3
"""Why is a near-complete exemplar Poor? A genome can be Poor on FEATURE flags
with no genome flag at all, so both have to be read."""
import json, os, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
meta = {l.split('\t')[0]: (l.split('\t')[1], int(l.rstrip('\n').split('\t')[6] or 0))
        for l in open('tob.full.tsv')}
for g in sys.argv[1:]:
    d = json.load(open('cov/work/Quality_GTOs/%s.qual.gto' % g))
    parts = list(d.get('genome_quality_flags') or [])
    for f in d.get('features') or []:
        if f.get('feature_quality_flags'):
            parts.append('%s: %s' % ((f.get('function') or '?')[:30],
                                     '; '.join(f['feature_quality_flags'])))
    print('%-12s %-40s %6d nt  %s'
          % (g, meta.get(g, ('', 0))[0][:40], meta.get(g, ('', 0))[1],
             ' | '.join(parts)[:120] or '(no flag recorded)'))
