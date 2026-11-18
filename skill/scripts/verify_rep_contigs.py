#!/usr/bin/env python3
"""Verify every reference genome by SEQUENCE, never by its deposited name.

Curator rule, 9 Nov 2026: a genome name is a claim exactly like an annotation
string. GenBank and RefSeq are benchmarks, not authorities. So a reference is
only allowed to represent a module if its SEQUENCE agrees.

The test: BLAST every installed reference against every other. A reference
belongs where its nearest neighbours are. If its best match in ANOTHER
module beats its best match inside its own module, the name put it in the
wrong place and it will capture that other module's genomes at routing time.

This is not hypothetical. `Aviadenovirus.4.dna`, deposited "Fowl aviadenovirus
D", is 99.8% identical to FAdV-4 and 72.2% to real FAdV-D. As an
Aviadenovirus reference it took 16 of 27 FAdV-C genomes away from
Aviadenovirus_FowlC, which routed 41%. Removing it alone gave 100%.

Scores use the annotator's own blastn settings, because a check that measures
something other than what routing measures proves nothing.
"""
import os, sys, glob, subprocess, collections, json

#  Where the rep contigs are. Defaults to the installed runtime so this runs on
#  any checkout; it used to name one build's scratch directory, which meant the
#  check could only ever be re-run by the machine that first wrote it.
R = (os.environ.get('LOWVAN_REP_CONTIGS')
     or os.path.join(os.environ.get(
         'LOWVAN_DATA_DIR',
         os.path.expanduser('~/lowvan_runtime/Viral_Annotation')),
         'Viral-Rep-Contigs'))
OUT = os.environ.get('LOWVAN_MISFILED_JSON', 'misfiled_refs.json')
#  The annotator's own blastn settings. A check that measures something other
#  than what routing measures proves nothing -- bare blastn runs megablast
#  (word size 28, reward 1) and under-reports routing about fourfold.
BLAST = ['-evalue', '0.5', '-reward', '2', '-penalty', '-3', '-word_size', '11',
         '-soft_masking', 'false']
#  The adenovirus split this was written for. Pass module names to check others.
MODS = sys.argv[1:] or ['Mastadenovirus', 'Mastadenovirus_TwoFiber',
                        'Mastadenovirus_Primate', 'Aviadenovirus',
                        'Aviadenovirus_FowlA', 'Aviadenovirus_FowlC',
                        'Aviadenovirus_NonFowl']

if not os.path.isdir(R):
    sys.exit('no rep contigs at %s\n'
             'set LOWVAN_DATA_DIR, or LOWVAN_REP_CONTIGS for the directory '
             'itself.' % R)

refs = []
for m in MODS:
    for f in sorted(glob.glob(f'{R}/{m}.*.dna')):
        #  "Aviadenovirus.4.dna" must not be read as module "Aviadenovirus_FowlC"
        base = os.path.basename(f)
        mod = base.rsplit('.', 2)[0]
        if mod == m:
            refs.append((m, base, f))
print('%d reference genome(s) across %d module(s)' % (len(refs), len(MODS)))

best = collections.defaultdict(dict)          # query -> module -> (bits, who)
for qm, qb, qf in refs:
    for sm, sb, sf in refs:
        if qb == sb:
            continue
        out = subprocess.run(['blastn', '-query', qf, '-subject', sf,
                              '-outfmt', '6 bitscore pident'] + BLAST,
                             capture_output=True, text=True).stdout
        top = 0.0, 0.0
        for line in out.splitlines():
            p = line.split('\t')
            if len(p) >= 2 and float(p[0]) > top[0]:
                top = float(p[0]), float(p[1])
        if top[0] > best[qb].get(sm, (0,))[0]:
            best[qb][sm] = (top[0], top[1], sb)

bad = []
print('\n%-30s %-28s %-28s' % ('reference', 'best inside its module', 'best in another module'))
for qm, qb, qf in refs:
    own = best[qb].get(qm, (0, 0, '-'))
    out = max(((v[0], v[1], v[2], k) for k, v in best[qb].items() if k != qm),
              default=(0, 0, '-', '-'))
    #  A flag needs HIGH foreign identity, not merely higher-than-own. Three
    #  shapes produce a false alarm otherwise, and all three appeared on the
    #  first run:
    #    * a module with ONE reference has nothing to compare inside itself,
    #      so "own" is 0 and everything looks foreign (Aviadenovirus_FowlA.1);
    #    * a genuinely divergent genome sits ~70% from everything, and which
    #      side wins by a few bits is noise (the falcon reference, 68.5 vs
    #      67.9, and two NonFowl ducks);
    #    * bitscore favours the longer alignment, so a longer distant match
    #      can outscore a shorter close one.
    #  A real misfile is not subtle: the FAdV-D-labelled FAdV-C genome sits at
    #  99.8% to the other module against 74.1% to its own.
    flag = out[1] >= 95.0 and out[1] > own[1] and own[0] > 0
    if flag:
        bad.append(dict(module=qm, ref=qb, own_bits=own[0], own_pident=own[1],
                        foreign_module=out[3], foreign_bits=out[0],
                        foreign_pident=out[1], foreign_ref=out[2]))
    print('%-30s %8.0f %5.1f%%            %8.0f %5.1f%% %-22s %s'
          % (qb, own[0], own[1], out[0], out[1], out[3], 'MISFILED' if flag else ''))

print('\n%d reference(s) sit closer to another module than to their own.' % len(bad))
for b in bad:
    print('   %s (%s): %.1f%% to %s vs %.1f%% inside %s'
          % (b['ref'], b['module'], b['foreign_pident'], b['foreign_module'],
             b['own_pident'], b['module']))
json.dump(bad, open(OUT, 'w'), indent=1)
print('wrote %s' % os.path.abspath(OUT))
