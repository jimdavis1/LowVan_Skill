#!/usr/bin/env python3
"""Install alignments that exist in a build workdir but never reached the repo.

The symptom this repairs: a feature's `Viral-PSSMs/<M>.pssms/<FEAT>/` holds a
profile named `<N>_<K>.pssm` for which `PSSM-Alignments/<M>/<FEAT>/<N>_<K>.fa`
does not exist. Those `_K` names come from the n-terminal reclustering pass,
which writes its split alignments to `reclustered_alis/` beside the feature's
`corrected_alis/`. The installer only ever copied `corrected_alis/`, so the
split alignments were dropped in transit while their profiles installed fine.

A module in that state cannot be rebuilt from its own alignments, which is the
premise `modules/` rests on -- the division of a cluster survives only inside
the profile. Verified across the rhabdovirus workdirs: for every feature,
`corrected_alis union reclustered_alis` equals the profile set exactly, so the
union is the alignment set that should have been installed.

This script is ADDITIVE. It never deletes. An installed alignment with no
profile may be a leftover or may be a curated file whose workdir copy was
cleaned; telling those apart needs evidence this script does not have, so it
reports them and leaves them alone.

    python3 install_missing_alignments.py --repo <VA checkout> \
            --search <dir to scan for workdirs> --module M [--module ...] [--apply]

Without --apply it reports and writes nothing.
"""
import argparse
import os
import shutil
import json

CAND_SUBDIRS = ("reclustered_alis", "corrected_alis", "alis")


def copy_verified(src, dst):
    """Box returns short and empty reads while exiting 0. Check the bytes."""
    want = os.path.getsize(src)
    if want == 0:
        raise IOError("source is empty, which on Box means unreadable: %s" % src)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copyfile(src, dst)
    got = os.path.getsize(dst)
    if got != want:
        raise IOError("short copy: %s is %d bytes, %s is %d" % (src, want, dst, got))
    return want


def index_workdirs(search_root, modules):
    """Map (module, feature, stem) -> source path, preferring reclustered_alis.

    A workdir looks like <anything>/<Module>/<FEAT>/<subdir>/<stem>.fa. Walking
    is done once because each stat over Box is slow.
    """
    want = set(modules)
    index = {}
    for dirpath, dirnames, filenames in os.walk(search_root):
        base = os.path.basename(dirpath)
        if base not in CAND_SUBDIRS:
            continue
        feat = os.path.basename(os.path.dirname(dirpath))
        mod = os.path.basename(os.path.dirname(os.path.dirname(dirpath)))
        if mod not in want:
            continue
        rank = CAND_SUBDIRS.index(base)
        for fn in filenames:
            if not fn.endswith(".fa"):
                continue
            key = (mod, feat, fn[:-3])
            prev = index.get(key)
            if prev is None or rank < prev[0]:
                index[key] = (rank, os.path.join(dirpath, fn))
    return {k: v[1] for k, v in index.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--search", required=True)
    ap.add_argument("--module", action="append", required=True)
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()

    master = json.load(open(os.path.join(a.repo, "Viral_PSSM.json")))
    print("indexing workdirs under %s ..." % a.search)
    index = index_workdirs(a.search, a.module)
    print("  %d candidate alignment files found\n" % len(index))

    installed = unresolved = orphan_ali = 0
    for mod in a.module:
        if mod not in master:
            print("%s: NOT INSTALLED" % mod)
            continue
        print("=== %s ===" % mod)
        feats = master[mod].get("features", {})
        for feat in sorted(feats):
            if feats[feat].get("special") in ("splice", "transcript_edit"):
                continue
            ad = os.path.join(a.repo, "PSSM-Alignments", mod, feat)
            pd = os.path.join(a.repo, "Viral-PSSMs", "%s.pssms" % mod, feat)
            alis = {x[:-3] for x in os.listdir(ad) if x.endswith(".fa")} if os.path.isdir(ad) else set()
            ps = {x[:-5].split(".")[-1] for x in os.listdir(pd)
                  if x.endswith(".pssm")} if os.path.isdir(pd) else set()
            missing = sorted(ps - alis)
            extra = sorted(alis - ps)
            if extra:
                orphan_ali += len(extra)
                print("  %-14s %d installed alignment(s) with no profile (left alone): %s"
                      % (feat, len(extra), extra[:8]))
            if not missing:
                continue
            got, lost = [], []
            for stem in missing:
                src = index.get((mod, feat, stem))
                if src is None:
                    lost.append(stem)
                    continue
                if a.apply:
                    copy_verified(src, os.path.join(ad, "%s.fa" % stem))
                got.append(stem)
            installed += len(got)
            unresolved += len(lost)
            print("  %-14s missing=%-3d resolved=%-3d UNRESOLVED=%s"
                  % (feat, len(missing), len(got), lost or "-"))
        print()

    print("%s %d alignment(s); %d unresolved; %d installed alignments without a profile"
          % ("installed" if a.apply else "WOULD install", installed, unresolved, orphan_ali))
    return 0


raise SystemExit(main())
