#!/usr/bin/env python3
"""Put back the generator's cutoff on features with fewer than five observed
calls.

Half the lowest accepted score is a defensible cutoff when there are enough
calls to know what "lowest" means. UNC4 was raised on two observations and
UNC5-UNC9 on one each: a single lineage-specific accessory seen once says
nothing about how a divergent member of the same orthogroup would score, and a
cutoff set from it would reject that member silently. Recall matters more than
precision here -- these features are already constrained by being
lineage-specific.
"""
import json, os
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
PRIOR = {('Serpentovirinae', 'UNC4'): 60, ('Serpentovirinae', 'UNC5'): 60,
         ('Serpentovirinae', 'UNC6'): 60, ('Serpentovirinae', 'UNC7'): 60,
         ('Serpentovirinae', 'UNC8'): 30, ('Serpentovirinae', 'UNC9'): 60}
for (M, k), v in PRIOR.items():
    p = '%s/%s_Viral_PSSM.json' % (M, M)
    d = json.load(open(p))
    cur = d[M]['features'][k].get('bit_cutoff')
    d[M]['features'][k]['bit_cutoff'] = v
    json.dump(d, open(p, 'w'), indent=3, sort_keys=True)
    print('%-16s %-6s %s -> %d (fewer than five observed calls)' % (M, k, cur, v))
