#!/usr/bin/env python3
"""Merge Extracted/ into the collections and refresh the module JSON.

    merge_and_json.py <builddir> <Genus> [patch.json]

Collections: existing  U  projected, deduplicated on sequence. A key that was
split (its Extracted file carries '|split' headers) is REPLACED by the split
set, because the old members are the lumped product.
JSON: for every key whose collection changed, min_len/max_len = 0.8x/1.2x and
bit_cutoff = max(20, 0.45x) of the new median -- the rule every existing
picornavirus feature was built with. patch.json adds or overrides features.
"""
import json, os, sys, glob, collections
B, G = sys.argv[1], sys.argv[2]
patch = json.load(open(sys.argv[3])) if len(sys.argv) > 3 else {}
C = os.path.join(B, "collections", G)
def fasta(p):
    out, h = collections.OrderedDict(), None
    for l in open(p):
        l = l.rstrip("\n")
        if l.startswith(">"): h = l[1:]; out[h] = []
        elif h is not None: out[h].append(l.strip())
    return [(k, "".join(v)) for k, v in out.items()]
jp = os.path.join(B, "%s_Viral_PSSM.json" % G)
J = json.load(open(jp)); F = J[G]["features"]
for k, v in patch.get("features", {}).items():
    F.setdefault(k, {}).update(v)
changed = {}
#  Every key that HAS a collection, not only the ones with an Extracted file.
#  A feature that STOPS being projected -- because its cleavage site turned out
#  not to be conserved across the references -- still needs its window
#  recomputed, or it keeps the one derived from the projected collection it no
#  longer has. Hepatovirus 2B and VPg, Kobuvirus 2A and VP3, and Sapelovirus 2A
#  are all in that position.
keys = sorted({os.path.basename(p)[:-6] for p in glob.glob(os.path.join(C, "*.fasta"))}
              | {os.path.basename(e)[len(G) + 1:-4]
                 for e in glob.glob(os.path.join(B, "Extracted", "%s.*.faa" % G))})
for k in keys:
    e = os.path.join(B, "Extracted", "%s.%s.faa" % (G, k))
    ext = fasta(e) if os.path.exists(e) else []
    cp = os.path.join(C, k + ".fasta")
    old = fasta(cp) if os.path.exists(cp) else []
    if not old and not ext: continue
    split = any("|split" in h for h, _ in ext)
    keep = [] if split else old
    seen = {s for _, s in keep}
    add = [(h, s) for h, s in ext if s not in seen]
    with open(cp, "w") as f:
        for h, s in keep + add: f.write(">%s\n%s\n" % (h, s))
    L = sorted(len(s) for _, s in keep + add); med = L[len(L) // 2]
    changed[k] = (len(old), len(keep + add), med)
    #  POLY is left alone deliberately: its collection never changes and the
    #  rebuild keeps its existing alignments, so its window stays as curated.
    if k in F and k != "POLY" and k not in patch.get("keep_cutoffs", []):
        lo, hi = int(round(0.8 * med)), int(round(1.2 * med))
        bit = max(20, int(round(0.45 * med)))
        #  ADDING TRAINING DATA MUST NOT REMOVE CALLABLE GENOMES. The 0.8x/1.2x
        #  rule is computed from the median, so a projection that brings in
        #  longer sequences drags the whole window up and genomes that used to
        #  sit inside it fall out the bottom: PICO_LDR_VP0 3A went 120-180 ->
        #  134-202 on a median moving 150 -> 168 and lost 13 calls, with more
        #  training data than before. So the window may only ever WIDEN here,
        #  and the bit cutoff may only ever fall.
        #  A SPLIT key is the exception: its collection was replaced, not
        #  extended, and the old window describes the lumped product.
        if not split and "min_len" in F[k]:
            lo, hi = min(lo, F[k]["min_len"]), max(hi, F[k]["max_len"])
            bit = min(bit, F[k].get("bit_cutoff", bit))
        F[k]["min_len"], F[k]["max_len"], F[k]["bit_cutoff"] = lo, hi, bit
def norm(x):
    if isinstance(x, dict): return {k: norm(v) for k, v in x.items()}
    if isinstance(x, list): return [norm(v) for v in x]
    if isinstance(x, float) and x == int(x): return int(x)
    return x
open(jp, "w").write(json.dumps(norm(J), indent=3, separators=(",", " : "), sort_keys=True, ensure_ascii=False) + "\n")
for k, (a, b, m) in sorted(changed.items()):
    f = F.get(k, {})
    print("  %-5s %4d -> %4d seqs  median %4d  window %s-%s  bit %s" % (k, a, b, m, f.get("min_len"), f.get("max_len"), f.get("bit_cutoff")))
