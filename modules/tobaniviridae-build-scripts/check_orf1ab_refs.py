#!/usr/bin/env python3
"""Compare each ORF1ab reference's translation to that genome's own deposited
pp1ab protein, where one exists. The junction is only proven when the product
matches the deposit."""
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
seq = {l.split('\t')[0]: l.rstrip('\n').split('\t')[1] for l in open('Tob.uniq.seq')}
dep = collections.defaultdict(list)
for l in open('cls/feature_locs.tsv'):
    p = l.rstrip('\n').split('\t')
    if p[8] in seq: dep[p[1]].append((p[7], seq[p[8]]))
refs = {}
h = None
for l in open('%s/Transcript-Editing/%s/ORF1AB.fasta' % (M, M)):
    if l.startswith('>'): h = l.split('|')[1].split('/')[0]; refs[h] = ''
    else: refs[h] += l.strip()
v = collections.Counter()
for g, s in sorted(refs.items()):
    prot = tr(s).rstrip('*')
    cands = [(a, d) for a, d in dep.get(g, []) if len(d) > 4000]
    if not cands: v['no_long_deposit'] += 1; continue
    a, d = max(cands, key=lambda x: len(x[1]))
    if d == prot: v['identical'] += 1
    elif d in prot or prot in d: v['nested len %+d' % (len(prot) - len(d))] += 1
    else: v['differs len %+d (%s)' % (len(prot) - len(d), a)] += 1
for k, n in v.most_common(): print('%-40s %d' % (k, n))
