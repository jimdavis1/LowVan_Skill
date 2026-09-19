#!/usr/bin/env python3
"""Rename leftover profiles from loN to plain integers.

build_leftover_pssms.py used to name second-pass clusters lo1, lo2 ... so they
could not collide with the main pass. Collision avoidance never needed a
prefix -- it only needed not restarting at 1 -- and the prefix put a
build-pipeline detail into a permanent identifier: the id is written into
every annotated genome as

    family_assignments -> ["LOWVAN", "<Module>.<FEAT>.<id>", ...]

so a second id shape existed only in modules this kit built. The generator now
emits integers; this migrates what is already on disk.

Renaming is per feature. Each loN takes the next integer above the highest one
the feature already uses, in loN order, so the leftovers stay contiguous and
keep their relative order. Ids like 2_1 or 3.1 (N-terminal splits) contribute
their leading integer to the maximum but are otherwise untouched.

Both trees move together, because a profile and its alignment are matched by
id and install_module.py now prunes any alignment that has no profile:

    <repo>/Viral-PSSMs/<Module>.pssms/<FEAT>/<Module>.<FEAT>.<id>.pssm
    <repo>/PSSM-Alignments/<Module>/<FEAT>/<id>.fa
    <workdir>/Alignments/<Module>/<FEAT>/{corrected_alis,reclustered_alis}/<id>.fa

The .pssm files are ASN.1 and carry no self-reference, so the filename is the
only thing that identifies them and a rename is sufficient.

Run --check first; it prints the full mapping and touches nothing.
"""
import argparse
import glob
import json
import os
import re
import sys


def lead(x):
    m = re.match(r"(\d+)", x)
    return int(m.group(1)) if m else 0


def plan_feature(ids):
    """Return {old_id: new_id} for the loN ids in this feature."""
    lo = sorted((i for i in ids if re.fullmatch(r"lo\d+", i)),
                key=lambda x: int(x[2:]))
    if not lo:
        return {}
    fixed = [i for i in ids if i not in lo]
    taken = set(fixed)
    n = max([lead(i) for i in fixed] or [0])
    out = {}
    for old in lo:
        n += 1
        while str(n) in taken:
            n += 1
        out[old] = str(n)
        taken.add(str(n))
    return out


def feature_ids(pssm_dir, module, feat):
    pre = "%s.%s." % (module, feat)
    ids = []
    for f in os.listdir(pssm_dir):
        if f.endswith(".pssm") and f.startswith(pre):
            ids.append(f[len(pre):-len(".pssm")])
    return ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True, help="Viral_Annotation checkout")
    ap.add_argument("--workdir", action="append", default=[],
                    help="a module working directory; repeatable")
    ap.add_argument("--module", action="append", default=[],
                    help="limit to these modules (default: every module with loN)")
    ap.add_argument("--check", action="store_true", help="print the plan, change nothing")
    a = ap.parse_args()

    pb = os.path.join(a.repo, "Viral-PSSMs")
    ab = os.path.join(a.repo, "PSSM-Alignments")
    if not os.path.isdir(pb):
        sys.exit("not a Viral_Annotation checkout: %s" % a.repo)

    #  workdir lookup by module name, so a module without one still migrates
    #  in the repo and simply reports the workdir as absent.
    wd = {}
    for w in a.workdir:
        for d in glob.glob(os.path.join(w, "Alignments", "*")):
            if os.path.isdir(d):
                wd[os.path.basename(d)] = d

    mods = sorted(os.path.basename(d)[:-len(".pssms")]
                  for d in glob.glob(os.path.join(pb, "*.pssms")))
    if a.module:
        mods = [m for m in mods if m in set(a.module)]

    total = 0
    renames = []          # (kind, src, dst)
    report = {}
    for m in mods:
        mp = os.path.join(pb, m + ".pssms")
        feats = sorted(f for f in os.listdir(mp) if os.path.isdir(os.path.join(mp, f)))
        for feat in feats:
            fp = os.path.join(mp, feat)
            ids = feature_ids(fp, m, feat)
            mapping = plan_feature(ids)
            if not mapping:
                continue
            report.setdefault(m, {})[feat] = mapping
            total += len(mapping)
            for old, new in mapping.items():
                renames.append(("pssm",
                                os.path.join(fp, "%s.%s.%s.pssm" % (m, feat, old)),
                                os.path.join(fp, "%s.%s.%s.pssm" % (m, feat, new))))
                ad = os.path.join(ab, m, feat)
                if os.path.isdir(ad):
                    src = os.path.join(ad, old + ".fa")
                    if os.path.exists(src):
                        renames.append(("align", src, os.path.join(ad, new + ".fa")))
                base = wd.get(m)
                if base:
                    for sub in ("corrected_alis", "reclustered_alis"):
                        src = os.path.join(base, feat, sub, old + ".fa")
                        if os.path.exists(src):
                            renames.append(("workdir",
                                            src, os.path.join(base, feat, sub, new + ".fa")))

    for m in sorted(report):
        print("  %s" % m)
        for feat in sorted(report[m]):
            mm = report[m][feat]
            print("    %-14s %s" % (feat, ", ".join(
                "%s->%s" % (k, mm[k]) for k in sorted(mm, key=lambda x: int(x[2:])))))
    print("\n  %d profile(s) to renumber, %d file(s) to move" % (total, len(renames)))
    if not a.workdir:
        print("  NOTE: no --workdir given; workdir alignments are left untouched")

    if a.check:
        print("  --check: nothing written")
        return

    #  Collisions are impossible by construction, but a rename that lands on an
    #  existing file would destroy it silently, so refuse rather than clobber.
    for _k, src, dst in renames:
        if os.path.exists(dst):
            sys.exit("refusing: %s already exists" % dst)

    #  Two phases via a temporary suffix. A direct rename could still collide
    #  if a feature ever mixed the two id shapes in a way the planner did not
    #  foresee; staging removes the ordering question entirely.
    for _k, src, dst in renames:
        os.rename(src, src + ".renumber_tmp")
    for _k, src, dst in renames:
        os.rename(src + ".renumber_tmp", dst)

    with open(os.path.join(a.repo, "LEFTOVER_RENUMBERING.json"), "w") as fh:
        json.dump(report, fh, indent=2, sort_keys=True)
    print("  done; mapping written to %s/LEFTOVER_RENUMBERING.json" % a.repo)


if __name__ == "__main__":
    main()
