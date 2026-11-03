#!/usr/bin/env python3
"""Lay out build/<MOD>/collections/<MOD>/ and write features.json per module.

build_features.py expects <workdir>/collections/<module>/<KEY>.fasta, which is
not the layout the triage writes, so the collections are copied into place
rather than the triage being bent to a downstream tool's shape.

PARAMS: ONLY -m IS TUNED. -m is a minimum CLUSTER SIZE and the reference says
plainly it "is a target, not a rule -- go down to -m 2 rather than lose a
taxon". -mi is a different kind of number: it is the identity at which two
sequences are declared THE SAME PROTEIN, and the curator settled it at 0.80.

An earlier version of this file copied the Picornaviridae pattern and dropped
small collections to mi 0.6, and a retry ladder then walked it lower still.
136 features across the four modules were built that way. That is not tuning,
it is redefining what counts as the same protein in order to make a build
succeed, and the kit's own MI_FLOOR = 0.6 is a hard backstop rather than a
licence to sit on it. Reverted: mi and mc stay at the documented 0.8.

Every departure from -m 5 is still recorded in BUILD_PARAMS.
"""
import json, os, shutil, sys

W = os.path.dirname(os.path.abspath(__file__)) or "."

#  ---------------------------------------------------------------------
#  ENFORCED: -mi and -mc are not tuneable from this script.
#
#  pipeline.md states it plainly -- "-mi is not a tuning dial for coverage;
#  it is the threshold at which you declare two sequences to be the same
#  protein. The documented default is 0.8." That file was READ during this
#  build and the value was lowered anyway, to 0.6, for 136 features across
#  four modules. A note in a doc did not stop it, so this is a check.
#
#  If a module's core genes will not build at 0.8, that is the module asking
#  to be split. Ask the curator; do not reach for the dial.
#  ---------------------------------------------------------------------
SETTLED = {"mi": "0.8", "mc": "0.8"}


def check_params(key, p):
    bad = {k: v for k, v in p.items()
           if k in SETTLED and str(v) != SETTLED[k]}
    if bad:
        raise SystemExit(
            "REFUSED: %s would be built with %s.\n"
            "  Minimum identity/coverage are settled at 0.8 and are not tuned\n"
            "  to make a build succeed. Tune -m instead, or split the module.\n"
            "  If the curator has agreed to a lower value, set it explicitly in\n"
            "  the module JSON and record the reason in BUILD_PARAMS."
            % (key, bad))

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
            p = {"m": 1}; tuned.append((k, n, med, p))
        elif n < 20:
            p = {"m": 2}; tuned.append((k, n, med, p))
        else:
            p = {}
        check_params(k, p)
        spec[k] = {"anno": f["anno"], "params": p}
    json.dump(spec, open(os.path.join(B, "features.json"), "w"), indent=1)
    with open(os.path.join(B, "BUILD_PARAMS"), "w") as fh:
        fh.write("# feature\tn_seqs\tmedian_aa\tparams -- ONLY -m is tuned;\n"
                 "# -mi and -mc stay at the settled 0.8 (identity is not a coverage dial)\n")
        for k, n, med, p in tuned:
            fh.write("%s\t%d\t%d\t%s\n" % (k, n, med, json.dumps(p)))
    print("%-18s %3d features, %d tuned" % (mod, len(spec), len(tuned)))
