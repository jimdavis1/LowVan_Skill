#!/usr/bin/env python3
"""Restore PMID_claude_generated to its categorised form, as a real dict.

check_dlits.py asks every model proposal to say what it supports -- 'function'
or 'coordinates' -- so the dict form is the right one. The problem was never
the shape in features.json; it was that gen_module_json passed the dict
through str() on its way into the module JSON, producing a list holding one
Python repr. A repr is not an id, so the ids were invisible to the very check
that verifies them.

This writes the dict straight into each module JSON after generation, which is
where the generator's own copy ends up anyway.
"""
import json, os
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
COORD = {'MPRO'}            # the only features whose citation fixes coordinates
for M in ('Torovirinae', 'Piscanivirinae', 'Remotovirinae', 'Serpentovirinae'):
    src = json.load(open('%s/features.json' % M))
    p = '%s/%s_Viral_PSSM.json' % (M, M)
    d = json.load(open(p))
    for k, v in d[M]['features'].items():
        ids = src.get(k, {}).get('PMID_claude_generated')
        if not ids: continue
        if isinstance(ids, dict):
            cat = {kk: list(vv) for kk, vv in ids.items()}
        else:
            cat = {'coordinates' if k in COORD else 'function': list(ids)}
        v['PMID_claude_generated'] = cat
    json.dump(d, open(p, 'w'), indent=3, sort_keys=True)
    print('%-16s %d feature(s) carry categorised proposals'
          % (M, sum(1 for v in d[M]['features'].values() if 'PMID_claude_generated' in v)))
