#!/usr/bin/env python3
"""Where a reference's translation and the genome's deposited pp1ab are the same
length but not identical, show the mismatching positions."""
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
refs = {}; h = None
for l in open('%s/Transcript-Editing/%s/ORF1AB.fasta' % (M, M)):
    if l.startswith('>'): h = l.split('|')[1].split('/')[0]; refs[h] = ''
    else: refs[h] += l.strip()
n = 0
for g, s in sorted(refs.items()):
    prot = tr(s).rstrip('*')
    cands = [(a, d) for a, d in dep.get(g, []) if len(d) > 4000 and len(d) == len(prot)]
    for a, d in cands:
        mis = [(j, d[j], prot[j]) for j in range(len(d)) if d[j] != prot[j]]
        if not mis: continue
        n += 1
        if n > 6: sys.exit()
        print('%s  %s  %d aa  %d mismatches  %s' % (g, a, len(d), len(mis), mis[:12]))
