#!/usr/bin/env python3
"""One metadata file per module: --meta SELECTS, so each rep-contig run must be
given only its own module's genomes."""
import os, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
rows = [l for l in open('budget/mod.full.tsv')]
for M in ('Torovirinae', 'Piscanivirinae', 'Remotovirinae', 'Serpentovirinae'):
    sel = [l for l in rows if l.split('\t')[3] == M]
    open('budget/%s.meta.tsv' % M, 'w').writelines(sel)
    print('%-16s %d records' % (M, len(sel)))
