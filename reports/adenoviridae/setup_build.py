#!/usr/bin/env python3
"""Lay out build/<MOD>/collections/<MOD>/ and write features.json per module.

build_features.py expects <workdir>/collections/<module>/<KEY>.fasta, which is
not the layout the triage writes, so the collections are copied into place
rather than the triage being bent to a downstream tool's shape.

PARAMS are scaled to the collection, exactly as the Picornaviridae build did:
a small or short collection cannot support the default -m 5 / -mi 0.8, and the
right response is to lower the threshold rather than lose the feature. Every
departure is recorded in BUILD_PARAMS next to the module.
"""
import json, os, shutil, sys

W = os.path.dirname(os.path.abspath(__file__)) or "."
MODS = ["Mastadenovirus", "Aviadenovirus", "Siadenovirus", "Barthadenovirus"]
for mod in MODS:
    B = os.path.join(W, "build", mod)
    cd = os.path.join(B, "collections", mod)
    os.makedirs(cd, exist_ok=True)
    #  Remove stale copies first -- the same defect that left 143 orphaned
    #  collections in the triage output would otherwise be reproduced here.
    for f in os.listdir(cd):
        if f.endswith(".fasta"):
            os.remove(os.path.join(cd, f))
    src = os.path.join(W, mod, "collections")
    for f in sorted(os.listdir(src)):
        if f.endswith(".fasta"):
            shutil.copy2(os.path.join(src, f), os.path.join(cd, f))
    J = json.load(open(os.path.join(W, "%s_Viral_PSSM.json" % mod)))[mod]["features"]
    spec, tuned = {}, []
    for k, f in sorted(J.items()):
        n, med = f["n_seqs"], f["median_aa"]
        if n < 10 or med < 20:
            p = {"m": 1, "mi": 0.6, "mc": 0.6}; tuned.append((k, n, med, p))
        elif n < 20:
            p = {"m": 2}; tuned.append((k, n, med, p))
        else:
            p = {}
        spec[k] = {"anno": f["anno"], "params": p}
    json.dump(spec, open(os.path.join(B, "features.json"), "w"), indent=1)
    with open(os.path.join(B, "BUILD_PARAMS"), "w") as fh:
        fh.write("# feature\tn_seqs\tmedian_aa\tparams (departure from -m 5 -mi 0.8 -mc 0.8)\n")
        for k, n, med, p in tuned:
            fh.write("%s\t%d\t%d\t%s\n" % (k, n, med, json.dumps(p)))
    print("%-18s %3d features, %d tuned" % (mod, len(spec), len(tuned)))
