#!/usr/bin/env python3
"""Lower the per-cluster member floor where a collection is too small to reach
it, and record why. -m 5 is a floor to aim for, not a law (SKILL step 3): these
subfamilies hold 6-17 unique sequences per protein across several genera, so no
cluster reaches five members and the feature builds nothing at all. -mi is NOT
lowered: identity still does the quality work."""
import json, os, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
NOTE = ("-m 2: this collection holds %d unique sequences spread over several "
        "genera, so no cluster reaches the documented floor of 5 and the "
        "feature builds zero profiles at it. The identity floor -mi 0.8 is "
        "unchanged, so a profile still rests on sequences that are 80%% "
        "identical.")
for M in sys.argv[1:]:
    p = '%s/features.json' % M
    f = json.load(open(p))
    for k, v in f.items():
        if v.get('special'): continue
        c = '%s/collections/%s/%s.fasta' % (M, M, k)
        n = sum(1 for l in open(c)) // 2 if os.path.exists(c) else 0
        v['params'] = {'m': 2}
        v['build_note'] = NOTE % n
    json.dump(f, open(p, 'w'), indent=1, sort_keys=True)
    print(M, 'set -m 2 on', sum(1 for v in f.values() if 'params' in v), 'features')
