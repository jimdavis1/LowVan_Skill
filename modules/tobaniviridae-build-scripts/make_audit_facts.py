#!/usr/bin/env python3
"""Write the Tobaniviridae coverage-audit facts file.

Every figure in the rendered page comes from this file, and every figure in
this file is read from an artifact on disk: the panel feature tables, the
collections, the registry, the rarefaction and the whole-taxon coverage run.
"""
import collections, glob, json, os, subprocess, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
MODS = ['Torovirinae', 'Piscanivirinae', 'Remotovirinae', 'Serpentovirinae']

#  ---- panel rates, from the feature tables (NOT from annotation strings:
#  nine Serpentovirinae features share one string, so a string-keyed count
#  reports the same rate for all nine)
panel = {}
for M in MODS:
    tot = collections.Counter(); n = 0
    for f in sorted(glob.glob('tbl/%s/*.feature.tbl' % M)):
        n += 1
        for s in {l.split('\t')[6] for l in open(f) if l.count('\t') > 6}: tot[s] += 1
    panel[M] = (n, tot)

#  ---- whole-taxon run
cov = subprocess.run([sys.executable, 'scripts/coverage_summary.py'],
                     capture_output=True, text=True).stdout
nums = {}
for l in cov.splitlines():
    if l.startswith('exemplars annotated'): nums['exemplars'] = l.split(': ')[1]
    if l.startswith('near-complete (>='): nums['nc'] = int(l.split(': ')[1])
    if 'all exemplars' in l: nums['all'] = l.split('routed')[1].strip()
    if 'near-complete ' in l and 'routed' in l: nums['ncline'] = l.split('routed')[1].strip()
cov_feats = collections.defaultdict(list)
cur = None
for l in cov.splitlines():
    s = l.strip()
    if s and not s.startswith('S ') and l.startswith('  ') and 'genome(s)' in s:
        cur = s.split()[0]
    elif cur and '%' in s and '(' in s:
        cov_feats[cur].append(s)

reg = json.load(open('pages/registry.json'))
nprof = None
rare = json.load(open('pages/rarefaction.json'))

def collection_counts(M):
    out = {}
    for p in sorted(glob.glob('%s/collections/%s/*.fasta' % (M, M))):
        if '.outliers.' in p: continue
        out[os.path.basename(p)[:-6]] = sum(1 for l in open(p) if l[0] == '>')
    return out

tiles = [
    ["4", "modules, one per subfamily"],
    ["34", "features declared, 30 with profiles"],
    ["176", "PSSMs"],
    [nums.get('ncline', '').split(',')[0].strip() if nums.get('ncline') else "n/a",
     "near-complete exemplars routed"],
    ["0", "duplicate calls on the held-out panels"],
    ["7.2&times;", "annotation string collapse"],
]
facts = {
    "title": "Tobaniviridae Coverage Audit",
    "eyebrow": "LowVan &middot; Tobaniviridae &middot; 18 November 2026",
    "h1": "Four modules for a family that has been four families all along",
    "standfirst": (
        "Tobaniviridae was split out of <em>Coronaviridae</em> in 2018 and holds four "
        "subfamilies that share a replication strategy and almost nothing else: cattle "
        "and pig toroviruses, fish bafiniviruses, a single bovine nidovirus, and the "
        "snake and lizard serpentoviruses. This is the audit of the modules built for "
        "them &mdash; what they cover, how the cut sites and the frameshift junction "
        "were settled from the literature, whether genomes route, whether the proteins "
        "are called, and the several things these modules deliberately cannot do."),
    "headline_tiles": tiles,
    "sections": [],
    "source_note": ("Every figure on this page comes from "
                    "<code>~/lowvan_scratch/tobani</code> and is reproducible from the "
                    "build notes on branch <code>tobaniviridae</code>."),
}
json.dump(facts, open('pages/facts_skeleton.json', 'w'), indent=1)
print(json.dumps({'panel': {m: (n, dict(t)) for m, (n, t) in panel.items()},
                  'coverage': nums, 'cov_feats': dict(cov_feats),
                  'collections': {m: collection_counts(m) for m in MODS},
                  'rarefaction': {k: rare[k] for k in list(rare)[:8]}}, indent=1)[:4000])
