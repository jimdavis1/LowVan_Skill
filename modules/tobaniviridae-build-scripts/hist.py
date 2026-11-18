#!/usr/bin/env python3
"""Histogram: annotation string x genus -> count, length range, median (unique proteins)."""
import collections, os, sys, re
M = sys.argv[1]
T = {'Torovirinae': 'Toro', 'Serpentovirinae': 'Serp', 'Piscanivirinae': 'Pisc', 'Remotovirinae': 'Remo'}[M]
d = os.path.expanduser('~/lowvan_scratch/tobani/dumps/' + M)
gen = {l.split('\t')[0]: l.rstrip('\n').split('\t')[3] for l in open('%s/%s.id_name_gs' % (d, T))}
ann = [l.rstrip('\n').split('\t') for l in open('%s/%s.uniq.id_ann' % (d, T))]
seq = [l.rstrip('\n').split('\t')[1] for l in open('%s/%s.uniq.seq' % (d, T))]
h = collections.defaultdict(list)
for (f, a), s in zip(ann, seq):
    g = gen[re.match(r'fig\|(\d+\.\d+)\.', f).group(1)]
    h[(a.lower(), g)].append(len(s))
for (a, g), v in sorted(h.items(), key=lambda x: (x[0][0], x[0][1])):
    v.sort()
    print('%4d  %-14s %5d-%5d med %5d  %s' % (len(v), g, v[0], v[-1], v[len(v) // 2], a))
