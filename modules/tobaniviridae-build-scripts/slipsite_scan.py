#!/usr/bin/env python3
"""Exploratory: in each near-complete in-scope genome, find the ORF1a reading
frame from sequence (blastx of the genome against that module's pp1a
deposits gives the frame and approximate start), find ORF1a's stop codon in
that frame, then list every slippery-heptamer candidate X_XXY_YYZ in phase 2 of
the ORF1a frame within 120 nt upstream of the stop, and test the -1 construct:
does duplicating the heptamer's last base open an ORF that runs >= 1500 codons?

Output: cls/slipsites.tsv  genome module heptamer pos(1-based, heptamer start)
        dist_to_ORF1a_stop  ORF1b_len_codons  context
"""
import collections, os, re, subprocess, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
CODON = {}
bases = 'TCAG'; aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
i = 0
for a in bases:
    for b in bases:
        for c in bases:
            CODON[a + b + c] = aas[i]; i += 1
def tr(s):
    return ''.join(CODON.get(s[j:j + 3], 'X') for j in range(0, len(s) - 2, 3))
mod = {}
for l in open('budget/mod.full.tsv'):
    p = l.rstrip('\n').split('\t'); mod[p[0]] = (p[3], p[1], int(p[6] or 0))
# pp1a-like deposits per module as blastx db (any >= 2500 aa protein, all deposits)
seq = {l.split('\t')[0]: l.rstrip('\n').split('\t')[1] for l in open('Tob.uniq.seq')}
with open('cls/rep_long.faa', 'w') as o:
    for l in open('cls/feature_locs.tsv'):
        p = l.rstrip('\n').split('\t')
        if p[8] in seq and len(seq[p[8]]) >= 2500:
            o.write('>%s\n%s\n' % (p[0], seq[p[8]]))
subprocess.run('makeblastdb -in cls/rep_long.faa -dbtype prot -out cls/rep_long >/dev/null', shell=True, check=True)
rows = []
for g, (m, name, ln) in sorted(mod.items()):
    if ln < 15000 or m.startswith('UNPL') or not os.path.exists('Contigs/%s.fna' % g): continue
    s = ''.join(x.strip() for x in open('Contigs/%s.fna' % g) if not x.startswith('>')).upper()
    # frame of the longest + strand ORF in the first 60% of the genome
    best = None
    for f in range(3):
        aa = tr(s[f:])
        pos = 0
        for orf in aa.split('*'):
            if best is None or len(orf) > best[2]:
                best = (f, pos, len(orf))
            pos += len(orf) + 1
    f, aapos, olen = best
    stop_nt = f + 3 * (aapos + olen)          # 0-based start of the stop codon
    if olen < 2500:
        rows.append((g, m, name, '-', '-', '-', '-', 'no ORF >= 2500 codons; longest %d' % olen)); continue
    cands = []
    for k in range(max(f, stop_nt - 150), stop_nt):
        if (k - f) % 3 != 2: continue
        h = s[k:k + 7]
        if len(h) == 7 and h[0] == h[1] == h[2] and h[3] == h[4] == h[5]:
            cons = s[f + 3 * aapos:k + 7] + s[k + 6:]
            prot = tr(cons)
            l1b = len(prot.split('*')[0]) - (k + 7 - (f + 3 * aapos)) // 3
            cands.append((h, k + 1, stop_nt - k, l1b, s[k - 9:k + 16]))
    if not cands:
        rows.append((g, m, name, '-', '-', '-', '-', 'no X_XXY_YYZ heptamer in phase 2 within 150 nt of the ORF1a stop'))
    for c in cands:
        rows.append((g, m, name) + c)
with open('cls/slipsites.tsv', 'w') as o:
    for r in rows: o.write('\t'.join(map(str, r)) + '\n')
summ = collections.Counter()
for r in rows:
    summ[(r[1], r[3], 'ORF1b>=1500' if str(r[6]).isdigit() and int(r[6]) >= 1500 else 'short/none')] += 1
for k, v in sorted(summ.items()): print('%-18s %-8s %-12s %d' % (k + (v,)))
