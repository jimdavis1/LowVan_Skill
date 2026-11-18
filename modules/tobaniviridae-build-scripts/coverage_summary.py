#!/usr/bin/env python3
"""Summarise the whole-taxon run with an explicit completeness cutoff.

"Near-complete" has to be measured against the genome the module expects, not
against a family-wide constant: a torovirus is 28 kb, a bostovirus 20 kb and a
veiled chameleon serpentovirus 36 kb, so one floor either admits partials or
excludes whole genera. Each exemplar is compared with the MEDIAN LENGTH OF ITS
OWN MODULE'S REFERENCE GENOMES, and an exemplar is near-complete when it
reaches 90% of that.
"""
import collections, glob, json, os
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
FRAC = 0.90
MODS = ['Torovirinae', 'Piscanivirinae', 'Remotovirinae', 'Serpentovirinae']
meta = {}
for l in open('tob.full.tsv'):
    p = l.rstrip('\n').split('\t'); meta[p[0]] = (p[1], int(p[6] or 0))
#  expected genome length per module = median of its own reference genomes
expect = {}
name2len = {}
for g, (nm, ln) in meta.items():
    name2len.setdefault(nm, ln)
acc2len = {}
for l in open('Contigs.metadata'):
    p2 = l.rstrip('\n').split('\t')
    if len(p2) >= 7 and p2[5].isdigit():
        for a in p2[6].split(','): acc2len[a] = int(p2[5])
for M in MODS:
    lens = []
    for f, c in json.load(open('%s/Rep-Contigs/close_genomes.json' % M)).items():
        nm = c.get('genome_name', '')
        if nm in name2len: lens.append(name2len[nm])
        elif c.get('genome_ids') in acc2len: lens.append(acc2len[c['genome_ids']])
    lens.sort(); expect[M] = lens[len(lens) // 2] if lens else 0

clu = collections.Counter()
for l in open('cov/work/clu_cluster.tsv'):
    clu[l.split()[0]] += 1
routed, qual, flags, feats = {}, {}, {}, {}
for q in sorted(glob.glob('cov/work/Quality_GTOs/*.qual.gto')):
    g = os.path.basename(q).split('.qual')[0]
    d = json.load(open(q))
    mods, ks = set(), set()
    for f in d.get('features') or []:
        fa = (f.get('family_assignments') or [[None, '']])[0]
        k = fa[1] if len(fa) > 1 else ''
        if k and '.' in k:
            mod, _, rest = k.partition('.')
            head, dot, tail = rest.rpartition('.')
            mods.add(mod); ks.add((mod, head if dot and tail.replace('_', '').isdigit() else rest))
        elif (f.get('function') or '').startswith('ORF1ab'):
            ks.add((sorted(mods)[0] if mods else '?', 'ORF1AB'))
    routed[g] = sorted(mods); qual[g] = (d.get('quality') or {}).get('genome_quality')
    flags[g] = d.get('genome_quality_flags') or []
    feats[g] = ks
allg = sorted(qual)
print('exemplars annotated: %d of %d clusters at 95%% identity' % (len(allg), len(clu)))
print('expected genome length per module (median of its reference genomes): %s'
      % ', '.join('%s %d nt' % (m, expect[m]) for m in MODS))
print('completeness cutoff: >= %d%% of that\n' % (100 * FRAC))
nc, part = [], []
for g in allg:
    mod = routed[g][0] if routed[g] else None
    exp = expect.get(mod)
    (nc if exp and meta.get(g, ('', 0))[1] >= FRAC * exp else part).append(g)
r = sum(1 for g in allg if routed[g])
print('routing: %d/%d exemplars routed (%.1f%%)' % (r, len(allg), 100.0 * r / len(allg)))
good_nc = sum(1 for g in nc if qual[g] == 'Good')
print('near-complete exemplars: %d; good %d (%.1f%%), poor %d'
      % (len(nc), good_nc, 100.0 * good_nc / max(1, len(nc)), len(nc) - good_nc))
good_p = sum(1 for g in part if qual[g] == 'Good')
print('partial exemplars:       %d; good %d, poor %d (%.1f%% correctly flagged incomplete)'
      % (len(part), good_p, len(part) - good_p,
         100.0 * (len(part) - good_p) / max(1, len(part))))
lens = sorted(meta.get(g, ('', 0))[1] for g in part)
if lens:
    print('   partial lengths: median %d nt, range %d-%d' % (lens[len(lens) // 2], lens[0], lens[-1]))
poor_nc = [g for g in nc if qual[g] != 'Good']
if poor_nc:
    print('\nnear-complete exemplars that came out POOR, with the reason:')
    for g in poor_nc:
        print('   %-12s %-42s %6d nt  %s'
              % (g, meta.get(g, ('', 0))[0][:42], meta.get(g, ('', 0))[1],
                 '; '.join(x[:60] for x in flags[g]) or '(no flag recorded)'))
print('\nper-feature calls among NEAR-COMPLETE exemplars, by module:')
bymod = collections.defaultdict(list)
for g in nc:
    for m in routed[g]: bymod[m].append(g)
for m, gs in sorted(bymod.items(), key=lambda x: -len(x[1])):
    f2 = collections.Counter()
    for g in gs:
        for mod, k in feats[g]:
            if mod == m: f2[k] += 1
    print('  %-16s %d genome(s)' % (m, len(gs)))
    for k, v in f2.most_common():
        print('      %-10s %5.1f%% (%d/%d)' % (k, 100.0 * v / len(gs), v, len(gs)))
