#!/usr/bin/env python3
"""Assign every genome in the dump to a module (subfamily) FROM SEQUENCE.

Seeds: near-complete genomes grouped by replicase identity (cls/groups40.tsv).
Every unique protein is blastp'd against the seed genomes' proteins; a genome
takes the module its proteins' best hits agree on. BV-BRC genus labels are
reported beside the call but never used to make it.
"""
import collections, os, re, subprocess
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
GROUP2MOD = {'G0': 'Torovirinae', 'G1': 'Serpentovirinae', 'G4': 'Serpentovirinae',
             'G5': 'Serpentovirinae', 'G2': 'Piscanivirinae', 'G3': 'Remotovirinae',
             'G6': 'TABLED_Glycera', 'G7': 'TABLED_Salmo'}
seedmod = {}
for l in open('cls/groups40.tsv'):
    p = l.rstrip('\n').split('\t'); seedmod[p[1]] = GROUP2MOD[p[0]]
meta = {}
for l in open('tob.full.tsv'):
    p = l.rstrip('\n').split('\t'); meta[p[0]] = p
md5s = [l.strip() for l in open('Tob.uniq.md5')]
seq = {l.split('\t')[0]: l.rstrip('\n').split('\t')[1] for l in open('Tob.uniq.seq')}
g_of_md5 = collections.defaultdict(set); md5_of_g = collections.defaultdict(set)
for l in open('Tob.id_md5'):
    f, m = l.rstrip('\n').split('\t')
    if not m: continue
    g = re.match(r'fig\|(\d+\.\d+)\.', f).group(1)
    g_of_md5[m].add(g); md5_of_g[g].add(m)
# seed db: proteins of seed genomes, labelled by module
with open('cls/seedprot.faa', 'w') as o, open('cls/allprot.faa', 'w') as a:
    for m in md5s:
        a.write('>%s\n%s\n' % (m, seq[m]))
        mods = {seedmod[g] for g in g_of_md5[m] if g in seedmod}
        if len(mods) == 1:
            o.write('>%s|%s\n%s\n' % (m, mods.pop(), seq[m]))
subprocess.run('makeblastdb -in cls/seedprot.faa -dbtype prot -out cls/seedprot > /dev/null', shell=True, check=True)
subprocess.run('blastp -query cls/allprot.faa -db cls/seedprot -evalue 1e-5 -outfmt "6 qseqid sseqid pident length bitscore" '
               '-max_target_seqs 5 -num_threads 8 -out cls/allprot_vs_seed.tsv', shell=True, check=True)
best = {}
for l in open('cls/allprot_vs_seed.tsv'):
    q, s, pid, ln, bit = l.split('\t')
    if q not in best or float(bit) > best[q][1]:
        best[q] = (s.split('|')[1], float(bit), float(pid))
gmod = {}
rows = []
for g, ms in md5_of_g.items():
    if g in seedmod:
        gmod[g] = (seedmod[g], 'seed'); continue
    votes = collections.Counter(best[m][0] for m in ms if m in best)
    if not votes:
        gmod[g] = ('UNPLACED', 'no protein hits a seed'); continue
    top, n = votes.most_common(1)[0]
    gmod[g] = (top if n == sum(votes.values()) else 'CONFLICT:' + str(dict(votes)), '%d/%d proteins' % (n, len(ms)))
with open('cls/genome_module.tsv', 'w') as o:
    for g in sorted(meta):
        mod, why = gmod.get(g, ('NO_PROTEINS', ''))
        p = meta[g]
        o.write('\t'.join([g, mod, why, p[3], p[4], p[1], p[6]]) + '\n')
c = collections.Counter((gmod.get(g, ('NO_PROTEINS',))[0], meta[g][3] or '-') for g in meta)
for k, v in sorted(c.items()): print('%-22s %-14s %d' % (k[0], k[1], v))
