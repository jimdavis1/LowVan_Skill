#!/usr/bin/env python3
"""Where a genome deposits a long pp1ab protein, find EVERY single-base
duplication near the ORF1a stop that reproduces that protein exactly, and print
the implied slip position and its heptamer context. This asks the deposits
where the ribosome slips instead of assuming it."""
import os, sys, collections
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
M = sys.argv[1]
CODON = {}; bb = 'TCAG'; aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
i = 0
for a in bb:
    for b in bb:
        for c in bb:
            CODON[a + b + c] = aas[i]; i += 1
tr = lambda s: ''.join(CODON.get(s[j:j + 3], 'X') for j in range(0, len(s) - 2, 3))
meta = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open('budget/mod.full.tsv')}
seq = {l.split('\t')[0]: l.rstrip('\n').split('\t')[1] for l in open('Tob.uniq.seq')}
dep = collections.defaultdict(list)
for l in open('cls/feature_locs.tsv'):
    p = l.rstrip('\n').split('\t')
    if p[8] in seq: dep[p[1]].append((p[7], seq[p[8]]))
hept = collections.Counter()
for g in sorted(meta):
    p = meta[g]
    if p[3] != M or not os.path.exists('Contigs/%s.fna' % g): continue
    longd = [(a, d) for a, d in dep.get(g, []) if len(d) > 4000]
    if not longd: continue
    a, d = max(longd, key=lambda x: len(x[1]))
    s = ''.join(x.strip() for x in open('Contigs/%s.fna' % g) if not x.startswith('>')).upper()
    hits = []
    for f in range(3):
        aa = tr(s[f:]); pos = 0
        for orf in aa.split('*'):
            if len(orf) >= 2500:
                start_nt = f + 3 * pos
                m = orf.find('M')
                if m >= 0:
                    st = start_nt + 3 * m
                    stop = start_nt + 3 * len(orf)
                    for k in range(max(st, stop - 200), stop):
                        cons = s[st:k + 1] + s[k:]
                        pr = tr(cons); e = pr.find('*')
                        if e > 0 and pr[:e] == d:
                            hits.append((k + 1, s[k - 6:k + 8]))
            pos += len(orf) + 1
    if hits:
        for k, ctx in hits[:2]:
            print('%-12s %-45s %-28s slip after nt %d   ...%s...' % (g, p[1][:45], a[:28], k, ctx))
            hept[ctx] += 1
    else:
        print('%-12s %-45s %-28s NO single-base duplication reproduces the deposit' % (g, p[1][:45], a[:28]))
print()
for k, v in hept.most_common(): print('%4d  %s' % (v, k))
