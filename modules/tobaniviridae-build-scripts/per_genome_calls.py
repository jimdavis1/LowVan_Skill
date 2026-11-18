#!/usr/bin/env python3
"""Per-genome called gene symbols, and the true per-FEATURE call rate.

run_gto_eval counts a feature as called when its ANNOTATION STRING appears in
the GTO, so every feature sharing "Uncharacterized lineage-specific protein"
is counted as called whenever any one of them is -- which is why a module with
nine such features reports the same rate for all nine. The feature tables
carry the gene symbol, so the honest rate is computed from those.
usage: per_genome_calls.py <Module>"""
import collections, glob, os, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
M = sys.argv[1]
meta = {l.split('\t')[0]: l.split('\t')[1] for l in open('Contigs.metadata')}
tot = collections.Counter(); n = 0
for f in sorted(glob.glob('tbl/%s/*.feature.tbl' % M)):
    g = os.path.basename(f)[:-len('.feature.tbl')]
    syms = [l.split('\t')[6] for l in open(f) if l.count('\t') > 6]
    n += 1
    for s in set(syms): tot[s] += 1
    print('%-12s %-42s %s' % (g, meta.get(g, '')[:42], ' '.join(syms)))
print('\n%d genome(s); per-feature call rate from the feature tables:' % n)
for k, v in tot.most_common():
    print('   %-10s %5.1f%%  (%d/%d)' % (k, 100.0 * v / n, v, n))
