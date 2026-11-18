#!/usr/bin/env python3
"""Contigs.metadata split per module, so a module's held-out panel holds only
its own genomes."""
import os
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
mod = {l.split('\t')[0]: l.split('\t')[3] for l in open('budget/mod.full.tsv')}
lines = open('Contigs.metadata').read().splitlines()
hdr, rows = lines[0], lines[1:]
for M in ('Torovirinae', 'Piscanivirinae', 'Remotovirinae', 'Serpentovirinae'):
    sel = [r for r in rows if mod.get(r.split('\t')[0]) == M]
    open('panel/%s.metadata' % M, 'w').write('\n'.join([hdr] + sel) + '\n')
    print('%-16s %d genomes' % (M, len(sel)))
