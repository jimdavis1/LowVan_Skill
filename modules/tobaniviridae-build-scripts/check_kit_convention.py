#!/usr/bin/env python3
"""Where does the kit's own hand-curated coronavirus ORF1ab reference put the
duplicated base? Align the reference to the genome it came from by finding the
one position where an extra base was inserted."""
import os, re, sys, requests
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
ref, h = {}, None
for l in open('kit/Transcript-Editing/Alphacoronavirus/ORF1ab.fasta'):
    if l.startswith('>'): h = l[1:].strip(); ref[h] = ''
    else: ref[h] += l.strip().upper()
h = next(iter(ref)); s = ref[h]
gid = re.search(r'(\d+\.\d+)', h).group(1)
q = "eq(genome_id,%s)&select(accession,sequence)&limit(50)" % gid
r = requests.post("https://www.bv-brc.org/api/genome_sequence/", data=q,
                  headers={"Content-Type": "application/rqlquery+x-www-form-urlencoded",
                           "Accept": "application/json"}, timeout=300).json()
g = max((x['sequence'] for x in r if isinstance(x, dict)), key=len).upper()
print(h[:70], 'ref %d nt, genome %d nt' % (len(s), len(g)))
# longest common prefix of the reference with the genome from the reference's start
i = g.find(s[:40])
print('reference starts at genome nt', i + 1)
j = 0
while j < len(s) and g[i + j] == s[j]: j += 1
print('first divergence at reference offset %d (genome nt %d)' % (j, i + j + 1))
print('genome   ...%s...' % g[i + j - 12:i + j + 12])
print('ref      ...%s...' % s[j - 12:j + 12])
# confirm the rest of the reference matches the genome shifted by one
k = 0
while j + k < len(s) and g[i + j + k - 1] == s[j + k]: k += 1
print('after the insertion the reference matches the genome shifted -1 for %d nt '
      '(reference tail %d nt)' % (k, len(s) - j))
