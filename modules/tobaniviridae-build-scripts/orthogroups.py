#!/usr/bin/env python3
"""Cluster a module's deposited proteins by homology and show each group with
its members' genome ORDER (which ORF it is, counting CDS 5'->3'), length and
the names submitters gave it. Position is used to SPLIT groups, never to call
them. usage: orthogroups.py <Module> [min_seq_id] [cov]"""
import collections, os, subprocess, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
M = sys.argv[1]
mi = sys.argv[2] if len(sys.argv) > 2 else '0.25'
cov = sys.argv[3] if len(sys.argv) > 3 else '0.5'
T = {'Torovirinae': 'Toro', 'Serpentovirinae': 'Serp', 'Piscanivirinae': 'Pisc', 'Remotovirinae': 'Remo'}[M]
gs = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open('dumps/%s/%s.id_name_gs' % (M, T))}
locs = {}
feats = collections.defaultdict(list)
for l in open('cls/feature_locs.tsv'):
    p = l.rstrip('\n').split('\t')
    if p[1] in gs and p[8]:
        locs[p[0]] = p; feats[p[1]].append(p)
rank = {}
for g, fs in feats.items():
    fs2 = [f for f in fs if int(f[6] or 0) < 2500]      # rank the structural block only
    for i, f in enumerate(sorted(fs2, key=lambda x: int(x[3]))): rank[f[0]] = i + 1
seq = {l.split('\t')[0]: l.rstrip('\n').split('\t')[1] for l in open('Tob.uniq.seq')}
d = 'orth/%s' % M; os.makedirs(d, exist_ok=True)
with open('%s/in.faa' % d, 'w') as o:
    for fid, p in sorted(locs.items()):
        if int(p[6] or 0) < 2500: o.write('>%s\n%s\n' % (fid, seq[p[8]]))
subprocess.run('mmseqs easy-cluster %s/in.faa %s/c %s/tmp --min-seq-id %s -c %s --cov-mode 0 '
               '-e 1e-3 >/dev/null 2>&1' % (d, d, d, mi, cov), shell=True, check=True)
cl = collections.defaultdict(list)
for l in open('%s/c_cluster.tsv' % d):
    a, b = l.split(); cl[a].append(b)
for i, (rep, mem) in enumerate(sorted(cl.items(), key=lambda x: -len(x[1])), 1):
    lens = sorted(len(seq[locs[m][8]]) for m in mem)
    ranks = collections.Counter(rank.get(m, '?') for m in mem)
    names = collections.Counter(locs[m][7] for m in mem)
    gen = collections.Counter(gs[locs[m][1]][3] for m in mem)
    print('C%-3d n=%-3d len %4d-%-4d med %4d  ORF#%s  %s' %
          (i, len(mem), lens[0], lens[-1], lens[len(lens) // 2],
           dict(sorted(ranks.items(), key=lambda x: str(x[0]))), dict(gen)))
    print('      names: %s' % '; '.join('%s x%d' % (a, c) for a, c in names.most_common(5)))
