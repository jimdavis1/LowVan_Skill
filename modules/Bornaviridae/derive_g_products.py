#!/usr/bin/env python3
"""Cut the G collection at the furin site into G_N and G_C.

Literature anchor: Richt 1998 (PMID 9557754) -- furin cleaves the BoDV-1 ORF-IV
product (p57/gp94) after the arginine at position 249.  Chemistry measured over
this module's own G collection: a single R-X-[K/R]-R site in the 40-62% window,
extended through the basic run, gives P1=Arg 150/150, P1'=Asp 149/150,
P4=Arg 150/150, and puts the BoDV-1 503 aa form's cut at exactly 249.
Members with no site, or more than one, are left out of both products.
"""
import re, collections, os, sys
pat = re.compile(r'(?=(R.[KR]R))')
W = 'collections/Bornaviridae'
seqs = {}; order = []; h = None
for l in open(W + '/G.fasta'):
    if l[0] == '>': h = l.rstrip('\n')[1:]; order.append(h); seqs[h] = ''
    else: seqs[h] += l.strip()
nout = collections.Counter()
with open(W + '/G_N.fasta', 'w') as fn, open(W + '/G_C.fasta', 'w') as fc, \
     open('collections/G_CUT.tsv', 'w') as rep:
    rep.write("id\tlen\tcut\tN_len\tC_len\tP6..P1\tP1'..P5'\n")
    for h in order:
        s = seqs[h]; lo, hi = int(.40 * len(s)), int(.62 * len(s))
        ends = []
        for m in pat.finditer(s):
            e = m.start() + 4
            if not (lo <= e <= hi): continue
            while e < len(s) and s[e] in 'KR': e += 1
            if e not in ends: ends.append(e)
        if len(ends) != 1:
            nout['no single site'] += 1; continue
        c = ends[0]
        if c < 120 or len(s) - c < 120:
            nout['product too short'] += 1; continue
        fn.write('>%s\n%s\n' % (h, s[:c]))
        fc.write('>%s\n%s\n' % (h, s[c:]))
        rep.write('%s\t%d\t%d\t%d\t%d\t%s\t%s\n' % (h.split()[0], len(s), c, c, len(s)-c, s[c-6:c], s[c:c+5]))
        nout['cut'] += 1
print(dict(nout))
