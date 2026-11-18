#!/usr/bin/env python3
"""Pick the panel genome for the genome map: every declared feature called
exactly once, and a well-known isolate preferred over a metagenomic contig.
usage: pick_map_genome.py <Module>"""
import collections, glob, json, os, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
M = sys.argv[1]
decl = set(json.load(open('%s/%s_Viral_PSSM.json' % (M, M)))[M]['features'])
meta = {l.split('\t')[0]: l.split('\t')[1] for l in open('Contigs.metadata')}
best = []
for q in sorted(glob.glob('gto_out/%s/*.qual.gto' % M)):
    g = os.path.basename(q)[:-len('.qual.gto')]
    d = json.load(open(q))
    keys = collections.Counter()
    for f in d.get('features') or []:
        fa = (f.get('family_assignments') or [[None, '']])[0]
        k = fa[1] if len(fa) > 1 else ''
        if k.startswith(M + '.'):
            rest = k[len(M) + 1:]
            head, dot, tail = rest.rpartition('.')
            keys[head if dot and tail.replace('_', '').isdigit() else rest] += 1
        elif (f.get('function') or '').endswith('polyprotein') and 'ORF1ab' in str(k or f.get('function')):
            keys['ORF1AB'] += 1
    once = sum(1 for k in decl if keys[k] == 1)
    dup = sum(1 for k in decl if keys[k] > 1)
    best.append((once, -dup, g, meta.get(g, ''), dict(keys)))
best.sort(reverse=True)
for once, ndup, g, name, keys in best[:6]:
    print('%-12s %-44s %d/%d features once, %d duplicated  %s'
          % (g, name[:44], once, len(decl), -ndup, sorted(keys)))
