#!/usr/bin/env python3
"""Split the family dump into one five-file dump per module.

Module membership comes from cls/genome_module.tsv (sequence-based). Intake
rule: organism names '<taxon> sp.' are out of scope, EXCEPT three ICTV
species-assigned isolates that BV-BRC files under a 'Serpentovirinae sp.'
organism name (see SP_EXCEPTIONS) -- recorded as an open decision.

uniq.* are written by walking one md5 list, so they cannot drift. The
representative feature for each md5 is re-chosen among IN-SCOPE genomes, and
its product is fetched (batched) so the annotation belongs to that feature.
"""
import collections, os, re, sys, time, urllib.parse, requests
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
API = "https://www.bv-brc.org/api"
HDR = {"Content-Type": "application/rqlquery+x-www-form-urlencoded", "Accept": "application/json"}
SP_EXCEPTIONS = {'3071291.3': 'Sertovirus cona', '3071292.3': 'Septovirus foka',
                 '3071293.3': 'Infratovirus latu'}
PREFIX = {'Torovirinae': 'Toro', 'Serpentovirinae': 'Serp', 'Piscanivirinae': 'Pisc',
          'Remotovirinae': 'Remo'}
meta, mod = {}, {}
for l in open('tob.full.tsv'):
    p = l.rstrip('\n').split('\t'); meta[p[0]] = p
for l in open('cls/genome_module.tsv'):
    p = l.rstrip('\n').split('\t'); mod[p[0]] = p[1]
def in_scope(g):
    if g in SP_EXCEPTIONS: return True
    return not re.search(r'\bsp\.($|\s)', meta[g][1])
seq = {l.split('\t')[0]: l.rstrip('\n').split('\t')[1] for l in open('Tob.uniq.seq')}
feats = []
for l in open('Tob.id_md5'):
    f, m = l.rstrip('\n').split('\t')
    feats.append((f, m, re.match(r'fig\|(\d+\.\d+)\.', f).group(1)))
# products for every feature, batched
prod = {}
if os.path.exists('cls/all_products.tsv'):
    for l in open('cls/all_products.tsv'):
        a, b = l.rstrip('\n').split('\t'); prod[a] = b
else:
    ids = [f for f, m, g in feats if m]
    for i in range(0, len(ids), 200):
        q = "in(patric_id,(%s))&select(patric_id,product)&limit(100000)" % ",".join(
            urllib.parse.quote(x, safe="") for x in ids[i:i + 200])
        r = requests.post(API + "/genome_feature/", data=q, headers=HDR, timeout=600).json()
        for x in [y for y in r if isinstance(y, dict)]: prod[x['patric_id']] = (x.get('product') or '').replace('\t', ' ')
    with open('cls/all_products.tsv', 'w') as o:
        for k, v in prod.items(): o.write('%s\t%s\n' % (k, v))
dropped = collections.Counter()
for M, T in PREFIX.items():
    gids = sorted(g for g in meta if mod.get(g) == M and in_scope(g))
    gset = set(gids)
    for g in meta:
        if mod.get(g) == M and g not in gset: dropped[(M, 'sp.')] += 1
    d = 'dumps/' + M
    os.makedirs(d, exist_ok=True)
    with open('%s/%s.id_name_gs' % (d, T), 'w') as o:
        for g in gids: o.write('\t'.join([g, meta[g][1], meta[g][2], meta[g][3] or 'unclassified']) + '\n')
    rep = {}
    with open('%s/%s.id_md5' % (d, T), 'w') as o:
        for f, m, g in feats:
            if g in gset:
                o.write('%s\t%s\n' % (f, m))
                if m and m not in rep: rep[m] = f
    md5s = sorted(rep)
    with open('%s/%s.uniq.md5' % (d, T), 'w') as a, open('%s/%s.uniq.id_ann' % (d, T), 'w') as b, \
         open('%s/%s.uniq.seq' % (d, T), 'w') as c:
        for m in md5s:
            a.write(m + '\n'); b.write('%s\t%s\n' % (rep[m], prod.get(rep[m], ''))); c.write('%s\t%s\n' % (m, seq[m]))
    print('%-16s %4d genomes  %4d unique proteins' % (M, len(gids), len(md5s)))
print('dropped as sp.:', dict(dropped))
