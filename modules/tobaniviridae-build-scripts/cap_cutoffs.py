#!/usr/bin/env python3
"""Cap every raised bit_cutoff at 500.

Half the lowest accepted score is the right SHAPE for a cutoff and the wrong
SIZE for a 4,500-residue polyprotein: Piscanivirinae ORF1A's panel calls score
6,960-9,747, so half is 3,000 -- and a divergent fish tobanivirus whose pp1a
matches over half its length scores far less than that. Measured on the
near-complete exemplars, raising ORF1A to 3,000 took Piscanivirinae ORF1A from
70% to 33%, rejecting partial matches that are informative calls (and that the
quality script already flags as short).

The false positive this whole exercise is aimed at was a 230-residue profile
scoring 30.3 bits inside a replicase ORF. 500 bits is more than an order of
magnitude above that and still well below any genuine call measured here, so
it keeps the protection and gives the divergent genomes back.
"""
import json, os
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
CAP = 500
for M in ('Torovirinae', 'Piscanivirinae', 'Remotovirinae', 'Serpentovirinae'):
    p = '%s/%s_Viral_PSSM.json' % (M, M)
    d = json.load(open(p))
    for k, v in d[M]['features'].items():
        if v.get('bit_cutoff', 0) > CAP:
            print('%-16s %-8s %5d -> %d' % (M, k, v['bit_cutoff'], CAP))
            v['bit_cutoff'] = CAP
    json.dump(d, open(p, 'w'), indent=3, sort_keys=True)
