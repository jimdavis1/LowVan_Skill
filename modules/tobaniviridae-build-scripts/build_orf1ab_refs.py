#!/usr/bin/env python3
"""Nucleotide references for the ORF1ab (pp1ab) transcript_edit feature.

Per in-scope near-complete genome of one module:
  frame of the longest + strand ORF (>=2500 codons) = the ORF1a frame;
  ORF1a start = first ATG after the preceding in-frame stop;
  slip site = the module's literature heptamer in phase 2 of that frame,
  within 150 nt of the ORF1a stop; construct duplicates the heptamer's last
  base (a -1 slip), runs on to the first stop and KEEPS that stop codon.
Rejects: ambiguous bases, ORF1b portion < 1500 codons, internal stop.

usage: build_orf1ab_refs.py <Module> <HEPTAMER>[,<HEPTAMER>...]
"""
import os, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
M, HEPT = sys.argv[1], sys.argv[2].split(',')
CODON = {}
bb, aas = 'TCAG', 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
i = 0
for a in bb:
    for b in bb:
        for c in bb:
            CODON[a + b + c] = aas[i]; i += 1
def tr(s):
    return ''.join(CODON.get(s[j:j + 3], 'X') for j in range(0, len(s) - 2, 3))
meta = {}
for l in open('budget/mod.full.tsv'):
    p = l.rstrip('\n').split('\t'); meta[p[0]] = p
outdir = '%s/Transcript-Editing/%s' % (M, M)
os.makedirs(outdir, exist_ok=True)
fa = open('%s/ORF1AB.fasta' % outdir, 'w')
log = open('%s/ORF1AB.build.tsv' % outdir, 'w')
log.write('genome\tname\tverdict\tstart\theptamer\thept_pos\tpp1a_aa\tpp1ab_aa\n')
kept = 0
for g in sorted(meta):
    p = meta[g]
    if p[3] != M or int(p[6] or 0) < 15000 or not os.path.exists('Contigs/%s.fna' % g):
        continue
    s = ''.join(x.strip() for x in open('Contigs/%s.fna' % g) if not x.startswith('>')).upper()
    best = None
    for f in range(3):
        aa, pos = tr(s[f:]), 0
        for orf in aa.split('*'):
            if best is None or len(orf) > best[2]: best = (f, pos, len(orf))
            pos += len(orf) + 1
    f, aapos, olen = best
    if olen < 2500:
        log.write('%s\t%s\tno_long_orf\t\t\t\t%d\t\n' % (g, p[1], olen)); continue
    orf_nt = f + 3 * aapos
    stop_nt = orf_nt + 3 * olen
    prot = tr(s[orf_nt:stop_nt])
    m = prot.find('M')
    if m < 0:
        log.write('%s\t%s\tno_met\t\t\t\t\t\n' % (g, p[1])); continue
    start = orf_nt + 3 * m
    hits = []
    for k in range(max(start, stop_nt - 150), stop_nt):
        if (k - f) % 3 != 2: continue
        if s[k:k + 7] in HEPT: hits.append(k)
    if not hits:
        log.write('%s\t%s\tno_slip_site\t%d\t\t\t%d\t\n' % (g, p[1], start + 1, (stop_nt - start) // 3)); continue
    k = hits[-1]
    #  Duplicate the LAST base of the slippery heptamer. That is the junction
    #  every RefSeq record of this family uses -- NC_007447 (Breda) heptamer
    #  14184..14190, CDS join(859..14190,14190..21059); NC_008516 (WBV)
    #  14549..14555, join(906..14555,14555..21523); NC_024709 (ball python
    #  nidovirus) 18405..18411, join(1018..18411,18411..25343); KM589359
    #  (bovine nidovirus) join(487..9825,9825..15818) -- and it is also the
    #  convention of the kit's own hand-curated coronavirus references
    #  (HCoV-229E, measured: one extra C after TTTAAAC).
    cons = s[start:k + 7] + s[k + 6:]
    pr = tr(cons)
    end = pr.find('*')
    if end < 0:
        log.write('%s\t%s\tno_stop\t%d\t%s\t%d\t\t\n' % (g, p[1], start + 1, s[k:k + 7], k + 1)); continue
    cons = cons[:3 * (end + 1)]
    body = pr[:end]
    n1b = end - (k + 7 - start) // 3
    if set(cons) - set('ACGT'):
        v = 'ambiguous_base'
    elif n1b < 1500:
        v = 'orf1b_short_%d' % n1b
    elif 'X' in body or '*' in body:
        v = 'internal_stop'
    else:
        v = 'OK'
    log.write('%s\t%s\t%s\t%d\t%s\t%d\t%d\t%d\n'
              % (g, p[1], v, start + 1, s[k:k + 7], k + 1, (stop_nt - start) // 3, end))
    if v == 'OK':
        acc = (p[12] or g).split(',')[0]
        fa.write('>%s|%s/%d-%d %s pp1ab [%s | %s]\n%s\n'
                 % (acc, g, start + 1, start + len(cons), M, p[1], g, cons))
        kept += 1
fa.close(); log.close()
print('%s: %d references written to %s/ORF1AB.fasta' % (M, kept, outdir))
os.system('cut -f3 %s/ORF1AB.build.tsv | tail -n +2 | sort | uniq -c' % outdir)
