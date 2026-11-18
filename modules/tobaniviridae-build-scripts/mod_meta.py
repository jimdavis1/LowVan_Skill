#!/usr/bin/env python3
"""tob.full.tsv with the genus column replaced by the sequence-assigned module.
In-scope records only (intake rule + the three species-assigned exceptions)."""
import os, re
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
SP = {'3071291.3', '3071292.3', '3071293.3'}
G2M = {'G0': 'Torovirinae', 'G1': 'Serpentovirinae', 'G4': 'Serpentovirinae', 'G5': 'Serpentovirinae',
       'G2': 'Piscanivirinae', 'G3': 'Remotovirinae', 'G6': 'TABLED_Glycera', 'G7': 'TABLED_Salmo'}
seed = {l.split('\t')[1]: G2M[l.split('\t')[0]] for l in open('cls/groups40.tsv')}
mod = {l.split('\t')[0]: l.split('\t')[1] for l in open('cls/genome_module.tsv')}
os.makedirs('budget', exist_ok=True)
n = 0
with open('budget/mod.full.tsv', 'w') as o:
    for l in open('tob.full.tsv'):
        p = l.rstrip('\n').split('\t')
        if re.search(r'\bsp\.($|\s)', p[1]) and p[0] not in SP: continue
        m = seed.get(p[0]) or mod.get(p[0], '')
        if m in ('NO_PROTEINS', 'UNPLACED', ''): m = 'UNPLACED'
        p[3] = m; o.write('\t'.join(p) + '\n'); n += 1
print(n, 'records')
