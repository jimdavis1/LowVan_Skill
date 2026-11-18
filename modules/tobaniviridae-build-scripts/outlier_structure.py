#!/usr/bin/env python3
"""Are the repeated-length outliers a distinct FORM or amplicon fragments?

A form is a real protein: it would start at a Met and either share the parent's
C-terminus (an alternative start) or start at the parent's residue 1. An
amplicon aligns to the MIDDLE of the parent. Each outlier is blasted against
its own feature's collection and the alignment's position in the subject is
reported. usage: outlier_structure.py <Module> <FEAT>"""
import collections, os, subprocess, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
M, F = sys.argv[1], sys.argv[2]
d = '%s/collections/%s' % (M, M)
sub = '%s/%s.fasta' % (d, F)
q = '%s/%s.outliers.fasta' % (d, F)
subprocess.run('makeblastdb -in %s -dbtype prot >/dev/null' % sub, shell=True, check=True)
out = subprocess.run('blastp -query %s -db %s -evalue 1e-5 -max_target_seqs 1 -num_threads 8 '
                     '-outfmt "6 qseqid qlen qstart qend sseqid slen sstart send pident"'
                     % (q, sub), shell=True, capture_output=True, text=True).stdout
best = {}
for l in out.splitlines():
    p = l.split('\t')
    if p[0] not in best: best[p[0]] = p
v = collections.Counter()
met = 0
seqs, h = {}, None
for l in open(q):
    if l.startswith('>'): h = l[1:].strip(); seqs[h] = ''
    else: seqs[h] += l.strip()
for qid, p in best.items():
    qlen, qs, qe, slen, ss, se = int(p[1]), int(p[2]), int(p[3]), int(p[5]), int(p[6]), int(p[7])
    if seqs.get(qid, 'X')[0] == 'M': met += 1
    if ss <= 3 and se >= slen - 3: v['spans the whole parent'] += 1
    elif ss <= 3: v["starts at the parent's N-terminus (truncated parent)"] += 1
    elif se >= slen - 3: v["ends at the parent's C-terminus (possible alternative start)"] += 1
    else: v['INTERNAL to the parent (an amplicon fragment)'] += 1
print('%s %s: %d outliers, %d with a hit to the collection, %d begin with Met'
      % (M, F, len(seqs), len(best), met))
for k, n in v.most_common(): print('   %-55s %d' % (k, n))
