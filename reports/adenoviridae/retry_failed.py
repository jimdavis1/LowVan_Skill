#!/usr/bin/env python3
"""Rebuild features that produced no PSSM, lowering the clustering identity.

WHY THEY FAILED, measured. Barthadenovirus E4_1 has 13 sequences and ALL 13
went to Leftover_Seqs.aa: at -mi 0.8 -mc 0.8 no two of them are 80% identical
over 80% coverage, so no cluster reached min_seqs and nothing was aligned. The
genus spans ruminant, avian and reptilian hosts, and 0.8 is simply the wrong
identity floor for it. 17 of 43 Barthadenovirus features failed this way,
including core genes -- DBP, IVa2, 52K, IIIa, protease, pVI, pX, UXP.

This is the documented Betarhabdovirinae case in references/pipeline.md, where
dropping -mi took P from 20% to 71% and M from 30% to 82%.

-mi is NOT a coverage dial. It is the identity at which two sequences are
declared to be the same protein, so it is lowered deliberately, one step at a
time, and the step that worked is recorded in BUILD_PARAMS next to the module.
The escalation stops at the FIRST setting that produces a profile rather than
going straight to the lowest, so no feature is clustered looser than it needs.
"""
import glob, json, os, subprocess, sys

W = os.path.dirname(os.path.abspath(__file__)) or "."
KIT = os.environ["LOWVAN_KIT"]
LADDER = [{"m": 2, "mi": 0.6, "mc": 0.7},
          {"m": 2, "mi": 0.5, "mc": 0.6},
          {"m": 2, "mi": 0.4, "mc": 0.5},
          {"m": 1, "mi": 0.3, "mc": 0.4}]


def built(mod, key):
    return glob.glob("%s/build/%s/Alignments/%s/%s/pssms/*.pssm" % (W, mod, mod, key))


def main():
    mod = sys.argv[1]
    B = "%s/build/%s" % (W, mod)
    J = json.load(open("%s/%s_Viral_PSSM.json" % (W, mod)))[mod]["features"]
    miss = [k for k in sorted(J) if not built(mod, k)]
    print("%s: %d feature(s) with no profile" % (mod, len(miss)))
    spec = json.load(open("%s/features.json" % B))
    rescued, failed = [], []
    for k in miss:
        for p in LADDER:
            spec[k] = {"anno": J[k]["anno"], "params": p}
            json.dump(spec, open("%s/features.json" % B, "w"), indent=1)
            subprocess.run([sys.executable,
                            "%s/skill/scripts/build_features.py" % KIT,
                            "--workdir", ".", "--module", mod,
                            "--features", "features.json", "--only", k],
                           cwd=B, capture_output=True, text=True)
            if built(mod, k):
                n = len(built(mod, k))
                rescued.append((k, p, n))
                print("  %-14s rescued at mi=%.1f mc=%.1f m=%d -> %d PSSM(s)"
                      % (k, p["mi"], p["mc"], p["m"], n))
                break
        else:
            failed.append(k)
            print("  %-14s STILL FAILED at every setting (n=%d, median %d aa)"
                  % (k, J[k]["n_seqs"], J[k]["median_aa"]))
    with open("%s/BUILD_PARAMS" % B, "a") as fh:
        fh.write("\n# rescued by lowering the clustering identity floor\n")
        for k, p, n in rescued:
            fh.write("%s\t%s\t%d PSSM(s)\n" % (k, json.dumps(p), n))
        for k in failed:
            fh.write("%s\tNO PROFILE AT ANY SETTING\n" % k)
    print("  %d rescued, %d still without a profile" % (len(rescued), len(failed)))


main()
