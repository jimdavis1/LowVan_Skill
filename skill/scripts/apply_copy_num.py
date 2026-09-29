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
    ap.add_argument("--contain", type=float, default=0.80,
                    help="a call at least this contained in a LONGER feature's "
                         "span does not count toward the call rate")
    ap.add_argument("--standalone", action="store_true",
                    help="decide on the STANDALONE rate, discounting calls "
                         "nested in a longer feature. Off by default because "
                         "genuine overlapping-ORF architecture is common; use "
                         "with --expected. The warning prints either way.")
    ap.add_argument("--expected", default="",
                    help="CHILD:PARENT overlaps that are real genome "
                         "architecture, e.g. CP:POLY,MP:POLY. Not discounted.")
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()

    blk = json.load(open(a.json))
    feats = blk[a.module]["features"]

    #  Two rates per feature, per genome:
    #    raw        the feature was called at all
    #    standalone it was called somewhere NOT contained in a longer
    #               feature's span
    #  They differ when a profile fires inside another gene. Closterovirus
    #  P20 is 89.3% raw and 53.3% standalone. But genuine overlapping ORFs
    #  are common -- Capillovirus CP sits inside its own polyprotein and is
    #  called on 100% of genomes with a 2% standalone rate -- so the two
    #  rates are used asymmetrically, see below.
    def span(f):
        l = (f.get("location") or [[None, 0, "+", 0]])[0]
        try:
            c, st0, strand, d0 = l[0], int(l[1]), l[2], int(l[3])
        except Exception:
            return None
        return (c, st0 - d0 + 1, st0) if strand == "-" else (c, st0, st0 + d0 - 1)

    expected = set()
    for pair in a.expected.split(","):
        if ":" in pair:
            c, par = pair.split(":", 1)
            expected.add((c.strip(), par.strip()))

    anno2key = {v["anno"]: k for k, v in feats.items()}
    raw = collections.Counter(); alone = collections.Counter(); routed = 0
    for q in sorted(glob.glob(os.path.join(a.eval, "*.qual.gto"))):
        try:
            d = json.load(open(q))
        except Exception:
            continue
        if not d.get("features"):
            continue
        routed += 1
        fs = [(anno2key.get(f.get("function") or ""), span(f)) for f in d["features"]]
        fs = [x for x in fs if x[0] and x[1]]
        got, got_alone = set(), set()
        for i, (k, sp) in enumerate(fs):
            got.add(k)
            nested = False
            for j, (k2, sp2) in enumerate(fs):
                if i == j or k2 == k or sp[0] != sp2[0]:
                    continue
                if (sp2[2] - sp2[1]) <= (sp[2] - sp[1]):
                    continue
                ov = min(sp[2], sp2[2]) - max(sp[1], sp2[1]) + 1
                mine = sp[2] - sp[1] + 1
                if mine > 0 and ov / mine >= a.contain and (k, k2) not in expected:
                    nested = True
                    break
            if not nested:
                got_alone.add(k)
        for k in got:
            raw[k] += 1
        for k in got_alone:
            alone[k] += 1
    called = raw

    if not routed:
        sys.exit("no scored genome in %s carried any feature" % a.eval)

    #  ASYMMETRIC, deliberately.
    #    ADD    needs BOTH rates over the threshold. A feature whose raw rate
    #           clears it only because the profile fires inside a longer gene
    #           would otherwise earn essentiality on another feature's locus.
    #    DROP   uses the RAW rate only. A low standalone rate is just as often
    #           real overlapping-ORF architecture, and silently demoting a
    #           correct feature is worse than leaving one un-demoted.
    add, drop, keep, blocked = [], [], [], []
    for k, v in sorted(feats.items()):
        frac = raw[k] / routed
        sa = alone[k] / routed
        has = bool(v.get("copy_num"))
        if frac >= a.threshold and not has:
            (add if sa >= a.threshold else blocked).append((k, frac, sa))
        elif frac < a.threshold and has:
            drop.append((k, frac))
        elif has:
            keep.append((k, frac))

    print("  %s: %d routed genome(s), threshold %.0f%%%s"
          % (a.module, routed, 100 * a.threshold,
             ""))
    #  Warn either way. A feature whose raw rate clears the threshold on
    #  another feature's locus is the failure this exists to prevent --
    #  Closterovirus P20 is 89.3% raw and 53.3% standalone -- but genuine
    #  overlapping ORFs are common (Capillovirus CP sits inside its own
    #  polyprotein), so the default does NOT silently demote them.
    print("\n  %-12s %7s  %s" % ("feature", "called", "copy_num"))
    for k in sorted(feats, key=lambda x: -called[x]):
        frac = called[k] / routed
        sa = alone[k] / routed
        mark = ("ADD" if any(x[0] == k for x in add) else
                "DROP" if (k, frac) in drop else
                "BLOCKED" if any(x[0] == k for x in blocked) else
                "keep" if feats[k].get("copy_num") else "-")
        print("  %-12s %6.1f%%  %s" % (k, 100 * frac, mark))

    if blocked:
        print("\n  NOT marked essential -- raw rate clears the threshold but the")
        print("  standalone rate does not, so the calls are on another feature's")
        print("  locus. Resolve with qc_duplicate_calls.py, or pass --expected")
        print("  CHILD:PARENT if the overlap is real genome architecture:")
        for k, f, sa in blocked:
            print("    %-11s %5.1f%% raw  ->  %5.1f%% standalone" % (k, 100 * f, 100 * sa))

    if a.write:
        for k, _f, _s in add:
            feats[k]["copy_num"] = 1
        for k, _f in drop:
            feats[k].pop("copy_num", None)
        write_json(blk, a.json)
        print("\n  wrote %s: +%d copy_num, -%d, %d blocked"
              % (a.json, len(add), len(drop), len(blocked)))
    else:
        print("\n  dry run -- rerun with --write (+%d, -%d, %d blocked)"
              % (len(add), len(drop), len(blocked)))


main()
