#!/usr/bin/env python3
"""Derive the main protease from pp1a at its published cleavage sites.

Torovirus (PMID 16571831, EToV Berne, GenBank DQ310701, numbering from Met1 of
pp1a): N terminus Gln3251|Ser3252 in the context SNFSFQ|SVFSKA, C terminus
Gln3543|Ser3544 in QPFKKQ|SVSNNV, both established by Edman degradation of the
cleavage products; mature Mpro = 292 residues. The same two sites were mapped
again by Edman in porcine torovirus (PMID 32727876, MH603532): Gln3095|Ser3096
in SNFSFQ|SVFSKV and Gln3382|Ser3383 in QPFKKQ|SVSTNV, mature 3CLpro = 287
residues.

Bafinivirus (PMID 21068254, white bream virus DF24/00, NC_008516, same
numbering convention): the protease uses alternative sites at both ends, all
four mapped by Edman. The predominant pair in infected cells is N1
Gln3442|Ala3443 in ITRQ|ALQR and C1 Gln3709|Thr3710 in VTQQ|TSVT, and the
paper's own N1-C1 marker protein runs Ala3443 to Gln3709.

So both ends of this feature are protease cuts with a published, conserved
context -- the case deriving-mature-peptides.md says to cut rather than model.
This script measures the junction chemistry across the whole collection BEFORE
cutting, and refuses any member where either motif is absent or ambiguous.

usage: derive_mpro.py <Module> [--write]
"""
import collections, os, re, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
M = sys.argv[1]
WRITE = '--write' in sys.argv
#  (N-terminal motif, C-terminal motif, expected product length range, where in
#  pp1a the N-terminal site sits as a fraction, from the reference virus)
SITES = {
 'Torovirinae':   (r'F[ST]FQ', r'FKKQ', (250, 340), 0.71),
 'Piscanivirinae': (r'[IVLT][TS][RK]Q', r'[VIL]T[QK]Q', (230, 300), 0.756),
}
nmot, cmot, (lo, hi), frac = SITES[M]


def read_fasta(p):
    h, s = None, []
    for l in open(p):
        if l.startswith('>'):
            if h: yield h, ''.join(s)
            h, s = l[1:].strip(), []
        else: s.append(l.strip())
    if h: yield h, ''.join(s)


rows = list(read_fasta('%s/collections/%s/ORF1A.fasta' % (M, M)))
chem = collections.Counter()
out, refused = [], collections.Counter()
for h, s in rows:
    #  search a window around the expected position, so a motif occurring
    #  elsewhere in a 4,500-residue polyprotein cannot be picked up
    w0, w1 = int(0.55 * len(s)), int(0.95 * len(s))
    nh = [m for m in re.finditer(nmot, s) if w0 <= m.start() <= w1]
    ch = [m for m in re.finditer(cmot, s) if w0 <= m.start() <= w1]
    if not nh: refused['no N-terminal motif in the window'] += 1; continue
    if not ch: refused['no C-terminal motif in the window'] += 1; continue
    n = nh[0].end()                       # first residue of the mature protein
    c = [m for m in ch if m.end() > n]
    if not c: refused['C-terminal motif upstream of the N-terminal one'] += 1; continue
    c = c[0].end()
    prod = s[n:c]
    if not (lo <= len(prod) <= hi):
        refused['derived length %d outside %d-%d' % (len(prod), lo, hi)] += 1; continue
    chem[('P1', s[n - 1])] += 1
    chem[('P1prime', s[n])] += 1
    chem[('P1_C', s[c - 1])] += 1
    chem[('P1prime_C', s[c])] += 1
    out.append((h, prod))
print('%s: %d pp1a sequences, %d derived, %d refused' % (M, len(rows), len(out), sum(refused.values())))
for k, v in refused.most_common(): print('   refused: %-45s %d' % (k, v))
n = max(1, len(out))
for pos in ('P1', 'P1prime', 'P1_C', 'P1prime_C'):
    d = {k[1]: v for k, v in chem.items() if k[0] == pos}
    print('   %-10s %s' % (pos, ', '.join('%s %d%%' % (a, round(100 * b / n))
                                          for a, b in sorted(d.items(), key=lambda x: -x[1]))))
if out:
    lens = sorted(len(s) for _h, s in out)
    print('   derived length %d-%d, median %d' % (lens[0], lens[-1], lens[len(lens) // 2]))
if WRITE and out:
    p = '%s/collections/%s/MPRO.fasta' % (M, M)
    with open(p, 'w') as fh:
        for h, s in out: fh.write('>%s\n%s\n' % (h.split()[0], s))
    print('   wrote', p)
