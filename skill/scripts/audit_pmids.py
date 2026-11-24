#!/usr/bin/env python3
"""Refuse to let a model-written citation reach the curated `PMID` field.

THE RULE, from Jim, 24 Nov 2026: **Claude never writes `PMID`.** Every citation
a model produces goes in `PMID_claude_generated`, however confident it is and
however directly it was read off the paper. The bare `PMID` field means exactly
one thing -- *a human curator looked this up and confirmed it supports the
claim* -- and that meaning is destroyed the moment a model writes into it. The
two fields are not two formats for the same information; they record who
verified it.

SKILL.md has said this since the kit existed, in three separate places, and it
did not hold. An audit on 24 Nov 2026 found **40** model-written citations
sitting in bare `PMID` across six modules: Alpharhabdovirinae (14),
Betarhabdovirinae (9), Dichorhavirus (6), Novirhabdovirus (6), Togaviridae (4),
Norovirus (1). A rule that is only written down is not an invariant. This is
the check that makes it one.

WHY THIS CANNOT SIMPLY FAIL ON ANY `PMID`. Jim curates: he reads the papers and
promotes citations from `PMID_claude_generated` into `PMID`, and then the module
is shipped again. A check that refused every `PMID` would refuse his work, which
is the opposite of protecting it.

So the test is not "is there a PMID" but "is this PMID NEW". A citation already
committed is the curator's and is left alone; one that appears in the working
tree and not in git HEAD was written by whoever is running now, which on this
path is a model. That is a precise distinction and it needs no annotation,
no honour system, and no memory of who typed what.

    python3 audit_pmids.py                      # every module in the kit
    python3 audit_pmids.py --module Orthohepadnavirus
    python3 audit_pmids.py --strict             # non-zero exit on any finding
    python3 audit_pmids.py --fix                # move new PMIDs into _claude_generated

Also checks the SHAPE of PMID_claude_generated. Two forms are in use and both
are fine -- a bare list of ids, and a {"function": [...], "coordinates": [...]}
object. What is not fine, and was found in both hepadnavirus modules, is a
Python repr inside a JSON string: ["{'function': ['399327']}"]. Nothing reads
that as a citation.
"""
import argparse
import ast
import glob
import json
import os
import subprocess
import sys


def git_show(repo, path):
    """The file as git HEAD has it, or None if it is not committed there."""
    try:
        r = subprocess.run(["git", "-C", repo, "show", "HEAD:" + path],
                           capture_output=True, text=True, timeout=60)
    except Exception:
        return None
    if r.returncode != 0 or not r.stdout.strip():
        return None
    try:
        return json.loads(r.stdout)
    except Exception:
        return None


def features_of(doc):
    for taxon, blk in (doc or {}).items():
        if not isinstance(blk, dict):
            continue
        for name, feat in (blk.get("features") or {}).items():
            if isinstance(feat, dict):
                yield taxon, name, feat


def normalise(value):
    """["{'function': ['1']}"] -> {"function": ["1"]}, or None if already sane."""
    if not (isinstance(value, list) and value):
        return None
    if not all(isinstance(x, str) and x.lstrip().startswith("{") for x in value):
        return None
    out = {}
    for x in value:
        try:
            d = ast.literal_eval(x)
        except Exception:
            return None
        if not isinstance(d, dict):
            return None
        for k, v in d.items():
            out.setdefault(k, []).extend(v if isinstance(v, list) else [v])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kit", default=os.environ.get("LOWVAN_KIT", "."))
    ap.add_argument("--module", help="one module instead of all of them")
    ap.add_argument("--strict", action="store_true",
                    help="exit non-zero if anything is wrong. Use this in a gate.")
    ap.add_argument("--fix", action="store_true",
                    help="move every NEW bare PMID into PMID_claude_generated "
                         "and repair stringified blocks")
    a = ap.parse_args()

    kit = os.path.abspath(a.kit)
    pat = os.path.join(kit, "modules",
                       a.module if a.module else "*", "*_Viral_PSSM.json")
    paths = sorted(glob.glob(pat))
    if not paths:
        raise SystemExit("no module JSON found under %s/modules" % kit)

    n_new = n_curated = n_shape = 0
    for path in paths:
        rel = os.path.relpath(path, kit)
        doc = json.load(open(path))
        head = git_show(kit, rel)
        head_pmids = {}
        for _, name, feat in features_of(head):
            if "PMID" in feat:
                head_pmids[name] = feat["PMID"]

        dirty = False
        for taxon, name, feat in features_of(doc):
            bad = normalise(feat.get("PMID_claude_generated"))
            if bad is not None:
                n_shape += 1
                print("  SHAPE    %-28s %-14s PMID_claude_generated is a Python "
                      "repr in a string" % (rel.split("/")[1], name))
                if a.fix:
                    feat["PMID_claude_generated"] = bad
                    dirty = True

            if "PMID" not in feat:
                continue
            was = head_pmids.get(name)
            if was is not None and was == feat["PMID"]:
                n_curated += 1
                continue
            #  New, or changed, since the last commit. On this path that means
            #  a model wrote it.
            n_new += 1
            print("  NEW PMID %-28s %-14s %s%s"
                  % (rel.split("/")[1], name, feat["PMID"],
                     "   (was %s)" % was if was is not None else ""))
            if a.fix:
                bare = feat.pop("PMID")
                bare = bare if isinstance(bare, list) else [bare]
                tgt = feat.get("PMID_claude_generated")
                if isinstance(tgt, dict):
                    have = tgt.setdefault("function", [])
                    have.extend(x for x in bare if x not in have)
                elif isinstance(tgt, list):
                    tgt.extend(x for x in bare if x not in tgt)
                else:
                    feat["PMID_claude_generated"] = {"function": list(bare)}
                dirty = True

        if dirty:
            json.dump(doc, open(path, "w"), indent=1)
            open(path, "a").write("\n")
            print("  FIXED    %s" % rel)

    print("\n%d module JSON(s) checked" % len(paths))
    print("  %3d citation(s) in a curated PMID field, unchanged since HEAD -- left alone"
          % n_curated)
    print("  %3d NEW or CHANGED bare PMID(s)" % n_new)
    print("  %3d malformed PMID_claude_generated block(s)" % n_shape)

    if n_new or n_shape:
        if not a.fix:
            print("\nA bare PMID means a human curator confirmed the paper supports\n"
                  "the claim. Nothing a model wrote may go there. Re-run with --fix,\n"
                  "or move these into PMID_claude_generated by hand.")
        return 1 if a.strict and not a.fix else 0
    print("\nclean: every citation written here is in PMID_claude_generated.")
    return 0


sys.exit(main())
