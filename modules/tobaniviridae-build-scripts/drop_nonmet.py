#!/usr/bin/env python3
"""Drop a collection member that does not begin at the protein's start.

check_met_start found one: Torovirinae M holds fig|74501.37.peg.1, a 224 aa
bovine torovirus deposit beginning WPFPDQ. The modal M is 233 aa beginning
MDFLS..., so this deposit is missing its first nine residues -- an unmarked
partial that the length window could not catch, because 224 aa is inside a
220-245 window set from the literature's 233. It had a single-sequence profile
of its own, which is a profile asserting that M starts at Trp.
"""
import glob, json, os, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
DROP = {('Torovirinae', 'M'): ['fig|74501.37.peg.1']}
for (M, K), ids in DROP.items():
    p = '%s/collections/%s/%s.fasta' % (M, M, K)
    rows, h = [], None
    for l in open(p):
        if l.startswith('>'): h = l[1:].strip(); rows.append([h, ''])
        else: rows[-1][1] += l.strip()
    keep = [r for r in rows if r[0].split()[0] not in ids]
    with open(p, 'w') as fh:
        for h, s in keep: fh.write('>%s\n%s\n' % (h, s))
    print('%s %s: %d -> %d sequences' % (M, K, len(rows), len(keep)))
    #  remove the profile and alignment built from it
    for ali in glob.glob('%s/Alignments/%s/%s/reclustered_alis/*.fa' % (M, M, K)):
        txt = open(ali).read()
        if any(i in txt for i in ids):
            n = os.path.basename(ali)[:-3]
            for q in glob.glob('%s/Alignments/%s/%s/pssms/*.%s.pssm' % (M, M, K, n)):
                os.remove(q); print('   removed profile', os.path.basename(q))
            os.remove(ali); print('   removed alignment', os.path.basename(ali))
