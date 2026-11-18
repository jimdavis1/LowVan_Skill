#!/usr/bin/env python3
"""Prepare genome classification: long replicase proteins (labelled by genome) and
the near-complete genome contigs, for blastx-based subfamily assignment."""
import re, os
W = os.path.expanduser('~/lowvan_scratch/tobani')
os.chdir(W)
os.makedirs('cls', exist_ok=True)
gen, name, length = {}, {}, {}
for l in open('tob.full.tsv'):
    p = l.rstrip('\n').split('\t')
    gen[p[0]] = p[3]; name[p[0]] = p[1]
    try: length[p[0]] = int(p[6])
    except: length[p[0]] = 0
md5s = [l.strip() for l in open('Tob.uniq.md5')]
ann = [l.rstrip('\n').split('\t') for l in open('Tob.uniq.id_ann')]
md5seq = {l.split('\t')[0]: l.rstrip('\n').split('\t')[1] for l in open('Tob.uniq.seq')}
annm = {m: a for m, (f, a) in zip(md5s, ann)}
n = 0
with open('cls/long.faa', 'w') as o:
    for l in open('Tob.id_md5'):
        f, m = l.rstrip('\n').split('\t')
        if not m: continue
        s = md5seq[m]
        if len(s) >= 1500 and re.search(r'1b|1ab|rdrp|replicase|polyprotein|pp1b|orf1b|nonstructural|1a', annm[m].lower()):
            gid = re.match(r'fig\|(\d+\.\d+)\.', f).group(1)
            o.write('>%s|%s\n%s\n' % (gid, f.split('|')[1], s)); n += 1
print('long proteins', n)
k = 0
with open('cls/nc.fna', 'w') as o:
    for g in sorted(length):
        if length[g] >= 15000 and os.path.exists('Contigs/%s.fna' % g):
            seq = ''.join(x.strip() for x in open('Contigs/%s.fna' % g) if not x.startswith('>'))
            o.write('>%s\n%s\n' % (g, seq)); k += 1
print('near-complete genomes', k)
