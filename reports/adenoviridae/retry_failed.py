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

#  Two dials, and they are not equivalent. -mi is the identity at which two
#  sequences are declared THE SAME PROTEIN, so it is lowered reluctantly.
#  -m is only the minimum cluster size, so lowering it costs much less --
#  but it costs something, because a one- or two-member profile is a lookup
#  table rather than a generalisation.
#
#  So at each identity level, try -m 2 first and -m 1 second, and only then
#  drop the identity. THE FIRST VERSION OF THIS LADDER lowered -m only on the
#  bottom rung, paired with mi 0.3, where everything over-merges -- and
#  Barthadenovirus LH1 "failed at every setting" while -m 1 at mi 0.5 gives it
#  5 clusters and 5 profiles. LH1 is 11 sequences in five same-species pairs
#  at 96-99% identity with under 50% between pairs: the pairs are real
#  clusters and only the size floor was rejecting them.
#
#  -mc MATTERS INDEPENDENTLY, and the second version of this ladder still
#  missed it. Barthadenovirus RH4 "failed at every setting" while mi 0.5 with
#  mc 0.5 gives it 2 profiles -- the ladder had only ever paired mi 0.5 with
#  mc 0.6. RH4's eight sequences are duplicate pairs across same-species
#  genomes: each tight cluster holds two IDENTICAL sequences which collapse to
#  one, and psiblast -in_msa cannot build a profile from a single sequence.
#  Loosening coverage merges the pairs into alignments that carry real
#  diversity. So the ladder walks coverage as well as identity.
#  ONLY -m MOVES. The curator settled minimum identity at 0.80, and -mi is
#  the identity at which two sequences are declared THE SAME PROTEIN -- not a
#  dial for making a build succeed. An earlier ladder here walked mi down to
#  0.6/0.5/0.4 and reported rescues; the kit silently clamped every value
#  below MI_FLOOR = 0.6 back up, so the log overstated what had been tried,
#  and the rescues that were real had redefined "same protein" to get there.
#
#  -m is a different kind of number: a minimum cluster size, which the
#  reference says outright is "a target, not a rule". Lowering it costs
#  generalisation -- a one- or two-member profile is closer to a lookup table
#  -- but it does not change what the profile claims to be.
#
#  A feature that cannot be profiled at mi 0.8 is REPORTED as unprofiled. If a
#  whole module needs a lower identity to build its CORE genes, that is not a
#  parameter problem; it is the module telling you it should be split.
LADDER = [{"m": 2}, {"m": 1}]


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
            check_params(k, p)
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
                print("  %-14s rescued at m=%d (mi stays 0.8) -> %d PSSM(s)"
                      % (k, p["m"], n))
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
