#!/usr/bin/env python3
"""Per orthogroup: members vs DISTINCT isolates. BV-BRC holds many isolates
twice (GenBank + RefSeq), so a 2-member cluster is often one sequence."""
import collections, os, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
M = sys.argv[1]
T = {'Torovirinae': 'Toro', 'Serpentovirinae': 'Serp', 'Piscanivirinae': 'Pisc', 'Remotovirinae': 'Remo'}[M]
gs = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open('dumps/%s/%s.id_name_gs' % (M, T))}
locs = {}
for l in open('cls/feature_locs.tsv'):
    p = l.rstrip('\n').split('\t')
    if p[1] in gs: locs[p[0]] = p
seq = {l.split('\t')[0]: l.rstrip('\n').split('\t')[1] for l in open('Tob.uniq.seq')}
cl = collections.defaultdict(list)
for l in open('orth/%s/c_cluster.tsv' % M):
    a, b = l.split(); cl[a].append(b)
for i, (rep, mem) in enumerate(sorted(cl.items(), key=lambda x: -len(x[1])), 1):
    iso = {gs[locs[m][1]][1].replace(' isolate', '').strip() for m in mem}
    uniq = {seq[locs[m][8]] for m in mem}
    names = collections.Counter(locs[m][7] for m in mem)
    print('C%-3d members=%-3d isolates=%-2d uniq_seqs=%-2d  %s' %
          (i, len(mem), len(iso), len(uniq), '; '.join('%s x%d' % x for x in names.most_common(3))))
