#!/usr/bin/env python3
"""Feature coordinates for every protein feature in the dump (batched)."""
import os, urllib.parse, requests
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
HDR = {"Content-Type": "application/rqlquery+x-www-form-urlencoded", "Accept": "application/json"}
ids = [l.split('\t')[0] for l in open('Tob.id_md5') if l.rstrip('\n').split('\t')[1]]
with open('cls/feature_locs.tsv', 'w') as o:
    for i in range(0, len(ids), 200):
        q = ("in(patric_id,(%s))&select(patric_id,genome_id,accession,start,end,strand,na_length,"
             "aa_length,product,aa_sequence_md5,segments)&limit(100000)") % ",".join(
            urllib.parse.quote(x, safe="") for x in ids[i:i + 200])
        r = requests.post("https://www.bv-brc.org/api/genome_feature/", data=q, headers=HDR, timeout=600).json()
        for x in r:
            if not isinstance(x, dict): continue
            o.write('\t'.join(str(x.get(k, '')) for k in ('patric_id', 'genome_id', 'accession', 'start', 'end',
                    'strand', 'aa_length', 'product', 'aa_sequence_md5', 'segments')) + '\n')
print(sum(1 for _ in open('cls/feature_locs.tsv')))
