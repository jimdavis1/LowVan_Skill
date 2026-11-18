#!/usr/bin/env python3
"""Print each in-scope genome's deposited protein features in genome order.
usage: genome_layouts.py <Module> [min_genome_len]"""
import collections, os, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
M = sys.argv[1]; minlen = int(sys.argv[2]) if len(sys.argv) > 2 else 15000
T = {'Torovirinae': 'Toro', 'Serpentovirinae': 'Serp', 'Piscanivirinae': 'Pisc', 'Remotovirinae': 'Remo'}[M]
gs = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open('dumps/%s/%s.id_name_gs' % (M, T))}
meta = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open('tob.full.tsv')}
feats = collections.defaultdict(list)
for l in open('cls/feature_locs.tsv'):
    p = l.rstrip('\n').split('\t')
    if p[1] in gs: feats[p[1]].append(p)
for g in sorted(gs, key=lambda x: (gs[x][3], gs[x][1])):
    if int(meta[g][6] or 0) < minlen or not feats[g]: continue
    print('%s  %s  [%s] %s nt' % (g, gs[g][1], gs[g][3], meta[g][6]))
    for p in sorted(feats[g], key=lambda p: (p[2], int(p[3]))):
        print('     %-12s %6s-%6s %s %5s aa  %s' % (p[2][:12], p[3], p[4], p[5], p[6], p[7]))
