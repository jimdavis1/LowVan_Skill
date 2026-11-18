#!/usr/bin/env python3
"""The source vocabulary, as BV-BRC holds it: one row per protein feature,
genome_id<TAB>product. This is the "old" side of the rarefaction curve."""
import os, re
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
prod = {}
for l in open('cls/all_products.tsv'):
    a, b = l.rstrip('\n').split('\t'); prod[a] = b
n = 0
with open('pages/bvbrc_products.tsv', 'w') as o:
    for l in open('Tob.id_md5'):
        f, m = l.rstrip('\n').split('\t')
        if not m or f not in prod: continue
        g = re.match(r'fig\|(\d+\.\d+)\.', f).group(1)
        o.write('%s\t%s\n' % (g, prod[f])); n += 1
print(n, 'protein features ->', 'pages/bvbrc_products.tsv')
