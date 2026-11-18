#!/bin/bash
#  Re-run the whole-taxon pass over the NEAR-COMPLETE exemplars only, against
#  the module as it now ships. The first pass ran before the bit cutoffs were
#  raised from the measured score distributions, so its numbers describe a
#  module that is not the one being committed; the near-complete exemplars are
#  the population section 08 of the audit quotes.
set -eu
W=$HOME/lowvan_scratch/tobani
cd "$W"
python3 - <<'PY'
import json, os, shutil, subprocess, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
out = subprocess.run([sys.executable, 'scripts/coverage_summary.py'],
                     capture_output=True, text=True).stdout
ids, seen = [], False
for l in out.splitlines():
    if l.startswith('near-complete exemplars:'): seen = True
    if l.startswith('   ') and seen and len(l.split()) > 2:
        pass
#  the near-complete set is recomputed here rather than scraped
meta = {l.split('\t')[0]: int(l.rstrip('\n').split('\t')[6] or 0) for l in open('tob.full.tsv')}
import glob, collections
expect = {}
name2len = {}
for l in open('tob.full.tsv'):
    p = l.rstrip('\n').split('\t'); name2len.setdefault(p[1], int(p[6] or 0))
for M in ('Torovirinae', 'Piscanivirinae', 'Remotovirinae', 'Serpentovirinae'):
    lens = sorted(name2len[c['genome_name']] for c in
                  json.load(open('%s/Rep-Contigs/close_genomes.json' % M)).values()
                  if c.get('genome_name') in name2len)
    expect[M] = lens[len(lens) // 2]
keep = []
for q in glob.glob('cov/work/Quality_GTOs/*.qual.gto'):
    g = os.path.basename(q).split('.qual')[0]
    d = json.load(open(q)); mod = None
    for f in d.get('features') or []:
        fa = (f.get('family_assignments') or [[None, '']])[0]
        k = fa[1] if len(fa) > 1 else ''
        if k and '.' in k: mod = k.split('.')[0]
    if mod and meta.get(g, 0) >= 0.9 * expect[mod]: keep.append(g)
d = 'cov/nc_contigs'
shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
for g in keep: shutil.copy('Contigs/%s.fna' % g, d)
print('%d near-complete exemplars staged' % len(keep))
PY
./r python3 kit/skill/scripts/evaluate_coverage.py --contigs cov/nc_contigs \
    --repo "$W/VA_tob" --work cov/nc_work --json VA_tob/Viral_PSSM.json \
    --out cov/nc_coverage.json --threads 6 --gto-retries 1
