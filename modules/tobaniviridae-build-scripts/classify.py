#!/usr/bin/env python3
"""Genome-to-genome replicase protein identity (blastx HSP-length-weighted), then
single-linkage groups at several identity thresholds. Labels are BV-BRC's; the
grouping is from sequence."""
import collections, os
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
gen, name, length = {}, {}, {}
for l in open('tob.full.tsv'):
    p = l.rstrip('\n').split('\t'); gen[p[0]] = p[3]; name[p[0]] = p[1]; length[p[0]] = p[6]
# accumulate identity over aligned length per (query genome, subject genome)
agg = collections.defaultdict(lambda: [0.0, 0])
for l in open('cls/nc_vs_long.tsv'):
    q, s, pid, ln, bit, qs, qe = l.split('\t')
    sg = s.split('|')[0]
    if sg == q: continue
    a = agg[(q, sg)]; a[0] += float(pid) * int(ln); a[1] += int(ln)
best = {}
for (q, sg), (w, ln) in agg.items():
    if ln < 1000: continue
    pid = w / ln
    if q not in best or pid > best[q][1]: best[q] = (sg, pid, ln)
gids = sorted(set(l.split('\t')[0] for l in open('cls/nc_vs_long.tsv')) | set(x.split('\t')[0] for x in open('cls/nc.fna') if False))
allq = [l[1:].strip() for l in open('cls/nc.fna') if l.startswith('>')]
sim = {}
for (q, sg), (w, ln) in agg.items():
    if ln >= 1000: sim[(q, sg)] = w / ln
def groups(th):
    par = {g: g for g in allq}
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for (q, s), v in sim.items():
        if v >= th and q in par and s in par: par[f(q)] = f(s)
    c = collections.defaultdict(list)
    for g in allq: c[f(g)].append(g)
    return list(c.values())
for th in (40, 60, 80, 90, 95):
    gs = groups(th)
    print('identity >= %d%%: %d groups' % (th, len(gs)))
print()
gs = sorted(groups(40), key=len, reverse=True)
with open('cls/groups40.tsv', 'w') as o:
    for i, g in enumerate(gs):
        lab = collections.Counter(gen[x] or '-' for x in g)
        print('G%d n=%d %s' % (i, len(g), dict(lab)))
        for x in g:
            b = best.get(x, ('-', 0, 0))
            o.write('G%d\t%s\t%s\t%s\t%s\t%s\t%.1f\n' % (i, x, gen[x], name[x], length[x], b[0], b[1]))
            print('    %s  %-14s %-55s %6s  best %s %.1f' % (x, gen[x], name[x][:55], length[x], b[0], b[1]))
