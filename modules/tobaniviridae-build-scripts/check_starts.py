#!/usr/bin/env python3
"""Compare the ORF1a start my reference builder picked with the start the
submitter deposited for pp1a/pp1ab on the same genome."""
import collections, os, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
M = sys.argv[1]
dep = collections.defaultdict(list)
for l in open('cls/feature_locs.tsv'):
    p = l.rstrip('\n').split('\t')
    try: aa = int(p[6])
    except: continue
    if aa >= 2500: dep[p[1]].append((int(p[3]), int(p[4]), aa, p[7]))
mine = {}
for l in open('%s/Transcript-Editing/%s/ORF1AB.build.tsv' % (M, M)):
    p = l.rstrip('\n').split('\t')
    if p[0] == 'genome' or len(p) < 5 or not p[3]: continue
    mine[p[0]] = (int(p[3]), p[2])
v = collections.Counter()
for g, (st, verdict) in sorted(mine.items()):
    if g not in dep: v['no deposit'] += 1; continue
    d = min(x[0] for x in dep[g])
    v['start %+d nt' % (st - d)] += 1
for k, n in v.most_common(): print('  %-16s %d' % (k, n))
