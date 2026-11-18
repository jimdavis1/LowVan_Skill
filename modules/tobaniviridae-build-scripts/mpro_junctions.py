#!/usr/bin/env python3
"""Are the called Mpro boundaries the published cleavage site?

The tiling test in check_cleaved_ends.py asks whether a precursor's products
sum to the precursor; it does not apply to a single INTERNAL mature peptide,
whose two offsets are large by construction. The test that does apply is
chemical: the protease cuts after Gln, so every called Mpro must begin at the
residue after a Gln and end in Gln. This reads the called proteins out of the
Quality GTOs and reports the first and last residue of each.
usage: mpro_junctions.py <Module>"""
import collections, glob, json, os, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
M = sys.argv[1]
first, last, lens = collections.Counter(), collections.Counter(), []
n = 0
for q in sorted(glob.glob('gto_out/%s/*.qual.gto' % M)):
    d = json.load(open(q))
    for f in d.get('features') or []:
        fa = (f.get('family_assignments') or [[None, '']])[0]
        k = fa[1] if len(fa) > 1 else ''
        if not k.startswith('%s.MPRO' % M): continue
        s = f.get('protein_translation') or f.get('aa_sequence') or ''
        if not s: continue
        n += 1; first[s[0]] += 1; last[s[-1]] += 1; lens.append(len(s))
print('%s: %d Mpro call(s) on the panel' % (M, n))
if not n: sys.exit()
lens.sort()
print('   length %d-%d, median %d  (literature: 292 aa EToV, 287 aa PToV, 267 aa WBV)'
      % (lens[0], lens[-1], lens[len(lens) // 2]))
print('   first residue: %s' % ', '.join('%s %d%%' % (a, round(100 * b / n)) for a, b in first.most_common()))
print('   last residue:  %s' % ', '.join('%s %d%%' % (a, round(100 * b / n)) for a, b in last.most_common()))
