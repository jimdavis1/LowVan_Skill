#!/usr/bin/env python3
"""Classify every panel genome where ORF1AB was not called.

Three causes, and they are not the same thing: ambiguous bases in the region
(the annotator demotes to partial_cds by design, because the gap-fill draws
from the subject); no reference within the 95% gate; or no pp1a region on the
contig at all. usage: te_misses.py <Module>"""
import glob, os, re, subprocess, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
M = sys.argv[1]
meta = {l.split('\t')[0]: l.split('\t')[1] for l in open('Contigs.metadata')}
ref = '%s/Transcript-Editing/%s/ORF1AB.fasta' % (M, M)
miss = []
for f in sorted(glob.glob('tbl/%s/*.feature.tbl' % M)):
    g = os.path.basename(f)[:-len('.feature.tbl')]
    syms = [l.split('\t')[6] for l in open(f) if l.count('\t') > 6]
    if 'ORF1ab' not in syms: miss.append((g, 'ORF1a' in syms))
print('%s: %d panel genome(s) without an ORF1AB call' % (M, len(miss)))
for g, has1a in miss:
    seq = ''.join(l.strip() for l in open('Contigs/%s.fna' % g) if not l.startswith('>')).upper()
    namb = sum(1 for c in seq if c not in 'ACGT')
    out = subprocess.run('blastn -query %s -subject Contigs/%s.fna -word_size 28 -outfmt '
                         '"6 qlen qstart qend pident"' % (ref, g), shell=True,
                         capture_output=True, text=True).stdout.splitlines()
    best = 0.0; bestcov = 0.0
    for l in out:
        qlen, qs, qe, pid = l.split('\t')
        cov = 100.0 * (int(qe) - int(qs) + 1) / int(qlen)
        if float(pid) > best: best, bestcov = float(pid), cov
    if not has1a: why = 'no pp1a called on this contig (partial genome)'
    elif namb: why = '%d ambiguous base(s) in the genome -> annotator demotes to partial_cds' % namb
    elif best < 95: why = 'best reference only %.1f%% identical (gate is 95%%)' % best
    else: why = 'reference at %.1f%%/%.0f%% coverage, cause not one of the three' % (best, bestcov)
    print('   %-12s %-42s %s' % (g, meta.get(g, '')[:42], why))
