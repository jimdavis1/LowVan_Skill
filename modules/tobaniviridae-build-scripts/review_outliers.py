#!/usr/bin/env python3
"""Read every <FEAT>.outliers.fasta the window produced, by genus and by length
CLASS. A length repeated across unrelated genomes is a FORM and must not be
dropped; a ragged spread of singletons is unmarked fragments. Then move the
files to <module>/collections/outliers/ so the length audit stops reading them
as collections. usage: review_outliers.py <Module> [--move]"""
import collections, glob, os, re, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
M = sys.argv[1]
T = {'Torovirinae': 'Toro', 'Serpentovirinae': 'Serp', 'Piscanivirinae': 'Pisc', 'Remotovirinae': 'Remo'}[M]
gs = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open('dumps/%s/%s.id_name_gs' % (M, T))}
for path in sorted(glob.glob('%s/collections/%s/*.outliers.fasta' % (M, M))):
    feat = os.path.basename(path).split('.')[0]
    rows, h = [], None
    for l in open(path):
        if l.startswith('>'): h = l[1:].strip(); rows.append([h, ''])
        else: rows[-1][1] += l.strip()
    cls = collections.Counter()
    for h, s in rows:
        g = re.match(r'fig\|(\d+\.\d+)\.', h)
        gen = gs.get(g.group(1), ['', '', '', '?'])[3] if g else '?'
        cls[(len(s), gen)] += 1
    print('%s %s: %d sequences' % (M, feat, len(rows)))
    multi = [(k, v) for k, v in cls.items() if v >= 3]
    for (ln, gen), v in sorted(multi, key=lambda x: -x[1])[:8]:
        print('   length %4d %-14s x%d  <- repeated length: check this is not a FORM' % (ln, gen, v))
    if not multi: print('   no length repeated 3+ times: a ragged spread, i.e. unmarked fragments')
if '--move' in sys.argv:
    d = '%s/collections/%s/outliers' % (M, M)
    os.makedirs(d, exist_ok=True)
    for path in glob.glob('%s/collections/%s/*.outliers.fasta' % (M, M)):
        os.rename(path, os.path.join(d, os.path.basename(path)))
    print('moved outliers to', d)
