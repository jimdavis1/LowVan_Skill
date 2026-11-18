#!/usr/bin/env python3
"""Merge the four modules' page inputs into one family-level set.

The four pages are published for the FAMILY, as Rhabdoviridae's were, because
a reader asks "does the Tobaniviridae module work", not "does Remotovirinae".
The registry JSON is keyed by module already, so merging is a dict update; the
string-collapse aggregation is keyed by annotation string, which several
modules share (Spike protein, Membrane glycoprotein ...), so those are summed.
"""
import collections, json, os
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
MODS = ['Torovirinae', 'Piscanivirinae', 'Remotovirinae', 'Serpentovirinae']
reg = {}
for M in MODS: reg.update(json.load(open('pages/registry_%s.json' % M)))
json.dump(reg, open('pages/registry.json', 'w'), indent=1)
print('registry: %d module(s), %d feature(s)'
      % (len(reg), sum(len(v.get('features') or v) for v in reg.values())))


def merge(name, key_is_string=True):
    out = None
    for M in MODS:
        p = '%s/%s' % (M, name)
        if not os.path.exists(p): continue
        d = json.load(open(p))
        if out is None:
            out = d if not isinstance(d, dict) else dict(d)
            continue
        if isinstance(d, dict) and isinstance(out, dict):
            for k, v in d.items():
                if k not in out: out[k] = v
                elif isinstance(v, list) and isinstance(out[k], list): out[k] = out[k] + v
                elif isinstance(v, dict) and isinstance(out[k], dict):
                    for kk, vv in v.items():
                        if kk in out[k] and isinstance(vv, (int, float)) and isinstance(out[k][kk], (int, float)):
                            out[k][kk] += vv
                        elif kk in out[k] and isinstance(vv, list) and isinstance(out[k][kk], list):
                            out[k][kk] = out[k][kk] + vv
                        else: out[k][kk] = vv
                elif isinstance(v, (int, float)) and isinstance(out[k], (int, float)): out[k] += v
        elif isinstance(d, list) and isinstance(out, list): out = out + d
    json.dump(out, open('pages/%s' % name, 'w'), indent=1)
    n = len(out) if out is not None else 0
    print('%-16s %d entries' % (name, n))


for f in ('agg.json', 'unbinned.json', 'typos.json'): merge(f)
