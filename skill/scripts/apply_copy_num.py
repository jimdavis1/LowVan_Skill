#!/usr/bin/env python3
"""Add copy_num to the features that earned it, after the module has run.

The safe build order is: ship with no copy_num, measure called/routed on real
genomes, then mark essential only what clears the threshold. copy_num is what
makes a feature essential, so a feature carrying it flags EVERY genome that
lacks it -- three modules scored 0.0% clean that way while their core
features called at 91-100%.

This reads the scored genomes an evaluation produced and rewrites the module
JSON, setting copy_num on features called on >= --threshold of routed
genomes and removing it from those below. It reports both directions and
writes nothing without --write.

    python3 apply_copy_num.py --eval Crinivirus/eval \
            --json Crinivirus/Crinivirus_Viral_PSSM.json --module Crinivirus
    python3 apply_copy_num.py ... --write
"""
import argparse, collections, glob, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from json_canon import write_json


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eval", required=True, help="directory of *.qual.gto")
    ap.add_argument("--json", required=True)
    ap.add_argument("--module", required=True)
    ap.add_argument("--threshold", type=float, default=0.85)
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()

    blk = json.load(open(a.json))
    feats = blk[a.module]["features"]

    called = collections.Counter(); routed = 0
    for q in sorted(glob.glob(os.path.join(a.eval, "*.qual.gto"))):
        try:
            d = json.load(open(q))
        except Exception:
            continue
        if not d.get("features"):
            continue
        routed += 1
        fns = {(f.get("function") or "") for f in d["features"]}
        for k, v in feats.items():
            if v.get("anno") in fns:
                called[k] += 1
    if not routed:
        sys.exit("no scored genome in %s carried any feature" % a.eval)

    add, drop, keep = [], [], []
    for k, v in sorted(feats.items()):
        frac = called[k] / routed
        has = bool(v.get("copy_num"))
        if frac >= a.threshold and not has:
            add.append((k, frac))
        elif frac < a.threshold and has:
            drop.append((k, frac))
        elif has:
            keep.append((k, frac))

    print("  %s: %d routed genome(s), threshold %.0f%%" % (a.module, routed, 100 * a.threshold))
    print("\n  %-12s %7s  %s" % ("feature", "called", "copy_num"))
    for k in sorted(feats, key=lambda x: -called[x]):
        frac = called[k] / routed
        mark = ("ADD" if (k, frac) in add else
                "DROP" if (k, frac) in drop else
                "keep" if feats[k].get("copy_num") else "-")
        print("  %-12s %6.1f%%  %s" % (k, 100 * frac, mark))

    if a.write:
        for k, _f in add:
            feats[k]["copy_num"] = 1
        for k, _f in drop:
            feats[k].pop("copy_num", None)
        write_json(blk, a.json)
        print("\n  wrote %s: +%d copy_num, -%d" % (a.json, len(add), len(drop)))
    else:
        print("\n  dry run -- rerun with --write (+%d, -%d)" % (len(add), len(drop)))


main()
