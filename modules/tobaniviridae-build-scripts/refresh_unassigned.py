#!/usr/bin/env python3
"""Rewrite UNASSIGNED_TRACKING.tsv after the homology rescue.

The engine writes that file during binning; `rescue_unassigned.py --write`
then appends the sequences it adopts to the collections but does not rewrite
the tracking file, so the file still lists strings that are now binned. A list
of unfinished work that overstates itself is as misleading as one that
understates. This recomputes it from what is actually in the collections now.
usage: refresh_unassigned.py <Module> <prefix>"""
import collections, glob, os, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
M, T = sys.argv[1], sys.argv[2]
binned = set()
for p in glob.glob('%s/collections/%s/*.fasta' % (M, M)):
    if '.outliers.' in p: continue
    for l in open(p):
        if l.startswith('>'): binned.add(l[1:].split()[0])
rows = collections.Counter(); uniq = collections.defaultdict(set)
md5s = [l.strip() for l in open('%s/%s.uniq.md5' % (M, T))]
ann = [l.rstrip('\n').split('\t') for l in open('%s/%s.uniq.id_ann' % (M, T))]
seq = [l.rstrip('\n').split('\t')[1] for l in open('%s/%s.uniq.seq' % (M, T))]
per = collections.Counter()
for l in open('%s/%s.id_md5' % (M, T)):
    p = l.rstrip('\n').split('\t')
    if len(p) > 1 and p[1].strip(): per[p[1].strip()] += 1
for m, (fid, a), s in zip(md5s, ann, seq):
    if fid in binned: continue
    rows[a] += per.get(m, 1); uniq[a].add(s)
out = '%s/collections/UNASSIGNED_TRACKING.tsv' % M
with open(out, 'w') as fh:
    fh.write('features\tuniq_seqs\tnot_already_called\tannotation\n')
    for a, c in rows.most_common():
        fh.write('%d\t%d\t%d\t%s\n' % (c, len(uniq[a]), len(uniq[a]), a))
print('%s: %d string(s) still unassigned after the rescue, %d sequence(s)'
      % (M, len(rows), sum(len(v) for v in uniq.values())))
