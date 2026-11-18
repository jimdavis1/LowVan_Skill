#!/usr/bin/env python3
"""Derive bafinivirus Mpro by anchoring on the protein whose boundaries were
sequenced, then snapping each end to the nearest Gln|(Ala,Thr,Ser) junction.

A bare motif search does not transfer across this subfamily: the white bream
virus N1 context ITRQ|ALQR is not conserved in Oncotshavirus, so a regex either
misses it or matches something else in a 4,500-residue polyprotein. The anchor
is the sequence the paper measured -- NC_008516 pp1a residues 3443-3709, the
N1-C1 marker protein of PMID 21068254 -- located in each pp1a by blastp, with
each boundary then moved to the nearest Gln whose next residue is Ala, Thr or
Ser, which is that protease's published consensus Gln|(Ala,Thr).

usage: derive_mpro_anchor.py [--write]
"""
import collections, os, subprocess, sys, tempfile
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
M = 'Piscanivirinae'
WRITE = '--write' in sys.argv
ANCHOR_GENOME = '405554.5'          # NC_008516, white bream virus DF24/00
N1, C1 = 3443, 3709                 # pp1a residue numbers, PMID 21068254
SNAP = 12                           # residues either side a boundary may move
OK_PRIME = set('ATS')


def read_fasta(p):
    h, s = None, []
    for l in open(p):
        if l.startswith('>'):
            if h: yield h, ''.join(s)
            h, s = l[1:].strip(), []
        else: s.append(l.strip())
    if h: yield h, ''.join(s)


rows = list(read_fasta('%s/collections/%s/ORF1A.fasta' % (M, M)))
anchor = None
for h, s in rows:
    if ANCHOR_GENOME in h:
        anchor = s[N1 - 1:C1]       # 3443..3709 inclusive
if anchor is None: sys.exit('anchor pp1a not in the collection')
print('anchor: white bream virus pp1a %d-%d, %d aa, starts %s ends %s'
      % (N1, C1, len(anchor), anchor[:8], anchor[-8:]))
tmp = tempfile.mkdtemp(prefix='mpro.')
open('%s/a.faa' % tmp, 'w').write('>anchor\n%s\n' % anchor)
with open('%s/q.faa' % tmp, 'w') as fh:
    for i, (h, s) in enumerate(rows): fh.write('>%d\n%s\n' % (i, s))
hits = subprocess.run('blastp -query %s/q.faa -subject %s/a.faa -evalue 1e-3 '
                      '-outfmt "6 qseqid qstart qend sstart send pident" -max_target_seqs 1'
                      % (tmp, tmp), shell=True, capture_output=True, text=True).stdout
best = {}
for l in hits.splitlines():
    q, qs, qe, ss, se, pid = l.split('\t')
    if q not in best: best[q] = (int(qs), int(qe), int(ss), int(se), float(pid))
out, chem, refused = [], collections.Counter(), collections.Counter()
for i, (h, s) in enumerate(rows):
    b = best.get(str(i))
    if not b: refused['no blastp hit to the measured protein'] += 1; continue
    qs, qe, ss, se, pid = b
    n0 = qs - (ss - 1)              # projected first residue of the mature protein
    c0 = qe + (len(anchor) - se)    # projected last residue
    def snap(pos, want_start):
        cands = []
        for d in range(0, SNAP + 1):
            for p in ((pos - d, pos + d) if d else (pos,)):
                if want_start:
                    if 2 <= p <= len(s) and s[p - 2] == 'Q' and s[p - 1] in OK_PRIME:
                        cands.append((d, p))
                else:
                    if 1 <= p < len(s) and s[p - 1] == 'Q' and s[p] in OK_PRIME:
                        cands.append((d, p))
            if cands: break
        return cands[0][1] if cands else None
    n = snap(n0, True); c = snap(c0, False)
    if n is None or c is None:
        refused['no Gln|(Ala,Thr,Ser) within %d residues of the projected end' % SNAP] += 1; continue
    prod = s[n - 1:c]
    if not (230 <= len(prod) <= 300):
        refused['derived length %d outside 230-300' % len(prod)] += 1; continue
    chem[('P1', s[n - 2])] += 1; chem[('P1prime', s[n - 1])] += 1
    chem[('P1_C', s[c - 1])] += 1; chem[('P1prime_C', s[c] if c < len(s) else '-')] += 1
    out.append((h, prod, n, c, pid))
print('%d pp1a sequences, %d derived, %d refused' % (len(rows), len(out), sum(refused.values())))
for k, v in refused.most_common(): print('   refused: %-55s %d' % (k, v))
n = max(1, len(out))
for pos in ('P1', 'P1prime', 'P1_C', 'P1prime_C'):
    d = {k[1]: v for k, v in chem.items() if k[0] == pos}
    print('   %-10s %s' % (pos, ', '.join('%s %d%%' % (a, round(100 * b / n))
                                          for a, b in sorted(d.items(), key=lambda x: -x[1]))))
for h, p, a, b, pid in out:
    print('   %-28s pp1a %4d-%4d  %3d aa  %.0f%% identical to the measured protein'
          % (h.split()[0][:28], a, b, len(p), pid))
if WRITE and out:
    path = '%s/collections/%s/MPRO.fasta' % (M, M)
    with open(path, 'w') as fh:
        for h, p, _a, _b, _pid in out: fh.write('>%s\n%s\n' % (h.split()[0], p))
    print('   wrote', path)
