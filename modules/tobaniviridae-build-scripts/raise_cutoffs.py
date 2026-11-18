#!/usr/bin/env python3
"""Raise each feature's bit_cutoff to half the lowest ACCEPTED score measured
on its own held-out panel, rounded down to a round number, and never below the
generator's value.

The generator derives a cutoff from the profile, and for a short protein that
lands at 30 bits. Measured: a second "Membrane glycoprotein" was called inside
Berne virus's replicase ORF at 30.34 bits against a cutoff of 30, while the
genuine M on that genome scored 476.8 -- a 15-fold margin. A cutoff is the
only thing standing between a 230-residue profile and a 28,000-nucleotide
genome, and 30 bits is not a threshold.

Half the lowest genuine score keeps every call this build has seen, with a
factor of two in hand for genomes it has not.
"""
import json, os, re, collections, glob
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
scores = collections.defaultdict(list)
for M in ('Torovirinae', 'Piscanivirinae', 'Remotovirinae', 'Serpentovirinae'):
    for f in glob.glob('gto_out/%s/_wd/*/*.stderr.txt' % M):
        for l in open(f, errors='replace'):
            m = re.search(r'Choosing\s+(\S+?)\.(\S+?)\.[\w_]+\.pssm\s+([0-9.]+)', l)
            if m and m.group(1) == M: scores[(M, m.group(2))].append(float(m.group(3)))


def round_down(x):
    for step in (500, 250, 100, 50, 25, 10):
        if x >= step: return int(x // step) * step
    return int(x)


for M in ('Torovirinae', 'Piscanivirinae', 'Remotovirinae', 'Serpentovirinae'):
    p = '%s/%s_Viral_PSSM.json' % (M, M)
    d = json.load(open(p)); feats = d[M]['features']
    for k, v in feats.items():
        obs = scores.get((M, k))
        if not obs or 'bit_cutoff' not in v: continue
        new = max(v['bit_cutoff'], round_down(min(obs) / 2.0))
        if new != v['bit_cutoff']:
            print('%-16s %-8s %4d -> %4d   (lowest accepted call %.1f, n=%d)'
                  % (M, k, v['bit_cutoff'], new, min(obs), len(obs)))
            v['bit_cutoff'] = new
    json.dump(d, open(p, 'w'), indent=3, sort_keys=True)
