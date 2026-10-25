#!/usr/bin/env python3
"""Report every installed feature whose alignments and profiles disagree.

A profile with no alignment cannot be rebuilt or audited -- `rebuild_pssms.py`
refreshes an existing profile from its alignment and cannot recreate one that
has none. An alignment with no profile is either a leftover from a build that
renumbered its cluster, or a curated file whose profile was dropped.

`install_module.py` copies both `corrected_alis/` and `reclustered_alis/` into
the flat `PSSM-Alignments/<M>/<FEAT>/`, but modules installed before that was
fixed are missing their n-terminal splits. This finds them. 25 Oct 2026: eight
modules, 207 alignments.

Run it before any kit sync. Exit status is non-zero if anything disagrees, so
it can gate a loop.

    python3 audit_alignment_parity.py --repo <VA checkout> [--module M ...]

`Leftover_Seqs.fa` is ignored: it holds the sequences that did not cluster, no
profile derives from it, and it is provenance rather than an alignment.
"""
import argparse
import json
import os

LEFTOVER = "Leftover_Seqs.fa"


def stems_fa(d):
    if not os.path.isdir(d):
        return set()
    return {x[:-3] for x in os.listdir(d) if x.endswith(".fa") and x != LEFTOVER}


def stems_pssm(d):
    if not os.path.isdir(d):
        return set()
    return {x[:-5].split(".")[-1] for x in os.listdir(d) if x.endswith(".pssm")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--module", action="append",
                    help="limit to these; default is every installed module")
    a = ap.parse_args()

    master = json.load(open(os.path.join(a.repo, "Viral_PSSM.json")))
    mods = a.module or sorted(master)

    bad_mods = 0
    tot_noali = tot_nopssm = 0
    for mod in mods:
        if mod not in master:
            print("%-26s NOT INSTALLED" % mod)
            bad_mods += 1
            continue
        feats = master[mod].get("features", {})
        rows = []
        for feat in sorted(feats):
            if feats[feat].get("special") in ("splice", "transcript_edit"):
                continue
            A = stems_fa(os.path.join(a.repo, "PSSM-Alignments", mod, feat))
            P = stems_pssm(os.path.join(a.repo, "Viral-PSSMs", "%s.pssms" % mod, feat))
            if A != P:
                rows.append((feat, sorted(P - A), sorted(A - P)))
        if not rows:
            continue
        bad_mods += 1
        print("=== %s ===" % mod)
        for feat, noali, nopssm in rows:
            tot_noali += len(noali)
            tot_nopssm += len(nopssm)
            print("  %-14s profiles with no alignment: %-3d %s" %
                  (feat, len(noali), noali[:8] if noali else "-"))
            if nopssm:
                print("  %-14s alignments with no profile: %-3d %s" %
                      ("", len(nopssm), nopssm[:8]))

    print("\n%d module(s) checked, %d with a disagreement; "
          "%d profile(s) with no alignment, %d alignment(s) with no profile"
          % (len(mods), bad_mods, tot_noali, tot_nopssm))
    if tot_noali:
        print("\nProfiles with no alignment are repairable from a build workdir:")
        print("  python3 install_missing_alignments.py --repo <repo> "
              "--search <workdir root> --module <M> [--apply]")
    return 1 if bad_mods else 0


raise SystemExit(main())
