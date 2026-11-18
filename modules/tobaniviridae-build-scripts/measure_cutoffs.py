#!/usr/bin/env python3
"""The bit score of every ACCEPTED call on the held-out panels, per feature.

gen_module_json derives a cutoff from the profile, and for a short protein that
can land at 30 bits -- low enough that a chance match inside the replicase ORF
clears it. Measured on Berne virus, a second "Membrane glycoprotein" was called
at 11,322 inside pp1a scoring 30.34 against a cutoff of 30, while the genuine M
on the same genome scored 476.8. This prints the genuine range so a cutoff can
be set with a margin instead of guessed.
"""
import collections, glob, os, re, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
scores = collections.defaultdict(list)
for M in ('Torovirinae', 'Piscanivirinae', 'Remotovirinae', 'Serpentovirinae'):
    for f in glob.glob('gto_out/%s/_wd/*/*.stderr.txt' % M):
        for l in open(f, errors='replace'):
            m = re.search(r'Choosing\s+(\S+?)\.(\S+?)\.[\w_]+\.pssm\s+([0-9.]+)', l)
            if m and m.group(1) == M: scores[(M, m.group(2))].append(float(m.group(3)))
print('%-16s %-8s %5s %9s %9s %9s' % ('module', 'feature', 'n', 'min', 'median', 'max'))
for (M, k), v in sorted(scores.items()):
    v.sort()
    print('%-16s %-8s %5d %9.1f %9.1f %9.1f'
          % (M, k, len(v), v[0], v[len(v) // 2], v[-1]))
