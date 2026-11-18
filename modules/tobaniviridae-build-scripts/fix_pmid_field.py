#!/usr/bin/env python3
"""Write PMID_claude_generated as a flat list of ids.

gen_module_json copies the field through from features.json, and features.json
carried it as {"function": [...], "coordinates": [...]} -- which came out of
the generator as a LIST HOLDING ONE STRINGIFIED DICT:

    "PMID_claude_generated": ["{'function': ['9188596', ...]}"]

Nothing errors. check_dlits reads ids out of that field, and a Python repr of
a dict is not an id, so every citation in these modules was invisible to the
check that exists to verify them. Every other module in the kit stores a flat
list (Bornaviridae G: ["9420287", "9557754", ...]), so that is the shape.

The per-key grouping is kept where it belongs -- in the feature's `note`,
which already says which paper establishes the function and which the
coordinates.
"""
import json, os, re
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))


def flatten(v):
    ids = []
    def walk(x):
        if isinstance(x, dict):
            for y in x.values(): walk(y)
        elif isinstance(x, list):
            for y in x: walk(y)
        elif isinstance(x, str):
            if re.fullmatch(r'\d{5,9}', x.strip()): ids.append(x.strip())
            else: ids.extend(re.findall(r'\b\d{7,9}\b', x))
    walk(v)
    seen, out = set(), []
    for i in ids:
        if i not in seen: seen.add(i); out.append(i)
    return out


for M in ('Torovirinae', 'Piscanivirinae', 'Remotovirinae', 'Serpentovirinae'):
    p = '%s/%s_Viral_PSSM.json' % (M, M)
    d = json.load(open(p))
    n = 0
    for k, v in d[M]['features'].items():
        if 'PMID_claude_generated' in v:
            flat = flatten(v['PMID_claude_generated'])
            if flat != v['PMID_claude_generated']:
                v['PMID_claude_generated'] = flat; n += 1
    json.dump(d, open(p, 'w'), indent=3, sort_keys=True)
    print('%-16s %d feature(s) rewritten' % (M, n))
    #  and in features.json, so a regeneration does not reintroduce it
    fp = '%s/features.json' % M
    f = json.load(open(fp))
    for k, v in f.items():
        if 'PMID_claude_generated' in v:
            v['PMID_claude_generated'] = flatten(v['PMID_claude_generated'])
    json.dump(f, open(fp, 'w'), indent=1, sort_keys=True)
