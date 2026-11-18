#!/usr/bin/env python3
"""Do the divergent Serpentovirinae orthogroups belong to S, M or N?

Anchors are the three proteins the literature NAMES in ball python nidovirus
NC_024709 (Stenglein 2014, PMID 25205093): ORF2 = S, ORF4 = M, ORF5a = N.
Every cluster representative is blasted against them; nothing is assigned on
position or on a submitter's string here.
"""
import collections, os, subprocess
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
gs = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open('dumps/Serpentovirinae/Serp.id_name_gs')}
locs = {}
for l in open('cls/feature_locs.tsv'):
    p = l.rstrip('\n').split('\t')
    if p[1] in gs: locs[p[0]] = p
seq = {l.split('\t')[0]: l.rstrip('\n').split('\t')[1] for l in open('Tob.uniq.seq')}
#  NC_024709 = genome 1986118.12; its S / M / N by the paper's ORF numbering,
#  identified here by their deposited coordinates, which match the record.
BP = {'S': (25340, 28219), 'M': (28977, 29624), 'N': (29637, 30095)}
anch = {}
for fid, p in locs.items():
    if p[1] == '1986118.12':
        for k, (s, e) in BP.items():
            if int(p[3]) == s and int(p[4]) == e: anch[k] = seq[p[8]]
print('anchors:', {k: len(v) for k, v in anch.items()})
os.makedirs('orth/anchor', exist_ok=True)
with open('orth/anchor/anchors.faa', 'w') as o:
    for k, s in anch.items(): o.write('>%s\n%s\n' % (k, s))
subprocess.run('makeblastdb -in orth/anchor/anchors.faa -dbtype prot >/dev/null', shell=True, check=True)
cl = collections.defaultdict(list)
for l in open('orth/Serpentovirinae/c_cluster.tsv'):
    a, b = l.split(); cl[a].append(b)
with open('orth/anchor/q.faa', 'w') as o:
    for i, (rep, mem) in enumerate(sorted(cl.items(), key=lambda x: -len(x[1])), 1):
        for m in mem: o.write('>C%d|%s\n%s\n' % (i, m, seq[locs[m][8]]))
subprocess.run('blastp -query orth/anchor/q.faa -db orth/anchor/anchors.faa -evalue 10 '
               '-outfmt "6 qseqid sseqid pident length evalue bitscore" -max_target_seqs 5 '
               '-num_threads 8 -out orth/anchor/hits.tsv', shell=True, check=True)
best = {}
for l in open('orth/anchor/hits.tsv'):
    q, s, pid, ln, ev, bit = l.rstrip('\n').split('\t')
    c = q.split('|')[0]
    if c not in best or float(bit) > best[c][2]: best[c] = (s, float(ev), float(bit), pid, ln)
for i, (rep, mem) in enumerate(sorted(cl.items(), key=lambda x: -len(x[1])), 1):
    c = 'C%d' % i
    b = best.get(c)
    lens = sorted(len(seq[locs[m][8]]) for m in mem)
    names = collections.Counter(locs[m][7] for m in mem).most_common(2)
    print('%-4s n=%-3d med %4d  %-28s  %s' %
          (c, len(mem), lens[len(lens) // 2], '; '.join('%s x%d' % x for x in names),
           ('best anchor %s  e=%.1g  %s%% over %s aa' % (b[0], b[1], b[3], b[4])) if b else 'NO HIT to S/M/N'))
