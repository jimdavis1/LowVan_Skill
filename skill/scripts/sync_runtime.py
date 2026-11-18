#!/usr/bin/env python3
"""Keep a local copy of the Viral_Annotation runtime, and top it up.

The runtime lives in Box, and Box serves it through placeholder files that
download on first read. When that stops working the failure is quiet and
misleading: directory listings return instantly, so everything looks
present, while every attempt to READ a .pssm times out. tblastn then gets a
short file and reports "Unsupported format for PSSM", which reads like a
corrupt profile. It is not -- the profile never arrived. Two modules were
written off as broken on that evidence before the pattern was spotted.

So the annotator should run against a local copy, and this builds one.

Two sources, because neither is complete on its own:

  - github.com/CEPI-dxkb/Viral_Annotation, which is authoritative but
    behind: it tracks 34 of the 110 module directories Box holds.
  - Box itself, for everything the repo does not carry.

Both are merged into one tree. The run is resumable and never overwrites a
module that is already complete, so it can be re-run whenever Box starts
cooperating and it will only fetch what is still short. A module is counted
complete only when its local file count matches Box's, which is the check
that the quiet failure defeats.

    sync_runtime.py --status          what is present, what is short
    sync_runtime.py --git             clone or update from the repo
    sync_runtime.py --from-box        try Box for whatever is still short
"""
import argparse, os, shutil, subprocess, sys

REPO = "https://github.com/CEPI-dxkb/Viral_Annotation.git"
DEFAULT_ROOT = os.path.expanduser("~/lowvan_runtime/Viral_Annotation")
DEFAULT_BOX = os.path.expanduser(
    "~/Library/CloudStorage/Box-Box/1_Projects/CEPI/LowVan_Skill_Development"
    "/_infrastructure/Viral_Annotation")


def nfiles(p):
    n = 0
    for _, _, fs in os.walk(p):
        n += len(fs)
    return n


def survey(root, box):
    """(complete, short, absent) module names, judged on file count."""
    src = os.path.join(box, "Viral-PSSMs")
    dst = os.path.join(root, "Viral-PSSMs")
    complete, short, absent = [], [], []
    if not os.path.isdir(src):
        return complete, short, absent
    for m in sorted(os.listdir(src)):
        if not m.endswith(".pssms"):
            continue
        try:
            want = nfiles(os.path.join(src, m))
        except OSError:
            want = -1
        have = nfiles(os.path.join(dst, m)) if os.path.isdir(
            os.path.join(dst, m)) else 0
        if want > 0 and have >= want:
            complete.append(m)
        elif have > 0:
            short.append((m, have, want))
        else:
            absent.append((m, want))
    return complete, short, absent


def do_git(root):
    tmp = root + ".git-tmp"
    shutil.rmtree(tmp, ignore_errors=True)
    print("cloning %s ..." % REPO)
    r = subprocess.run(["git", "clone", "--depth", "1", REPO, tmp])
    if r.returncode:
        sys.exit("clone failed")
    os.makedirs(root, exist_ok=True)
    moved = 0
    for dirpath, _, files in os.walk(tmp):
        if ".git" in dirpath.split(os.sep):
            continue
        rel = os.path.relpath(dirpath, tmp)
        out = os.path.join(root, rel) if rel != "." else root
        os.makedirs(out, exist_ok=True)
        for f in files:
            d = os.path.join(out, f)
            if not os.path.exists(d):          # never clobber a Box-sourced file
                shutil.copy2(os.path.join(dirpath, f), d)
                moved += 1
    shutil.rmtree(tmp, ignore_errors=True)
    print("added %d file(s) from the repo" % moved)


def do_box(root, box):
    """Copy only what is still short. Box failures are per-file and expected."""
    complete, short, absent = survey(root, box)
    todo = [m for m, _, _ in short] + [m for m, _ in absent]
    if not todo:
        print("nothing short; local copy already matches Box")
        return
    print("%d module(s) to fetch from Box" % len(todo))
    ok = fail = 0
    for m in todo:
        s = os.path.join(box, "Viral-PSSMs", m)
        d = os.path.join(root, "Viral-PSSMs", m)
        got = err = 0
        for dirpath, _, files in os.walk(s):
            rel = os.path.relpath(dirpath, s)
            out = os.path.join(d, rel) if rel != "." else d
            os.makedirs(out, exist_ok=True)
            for f in files:
                t = os.path.join(out, f)
                if os.path.exists(t):
                    continue
                try:
                    shutil.copy2(os.path.join(dirpath, f), t)
                    got += 1
                except OSError as e:           # the placeholder did not hydrate
                    err += 1
                    if os.path.exists(t):
                        os.remove(t)           # never leave a short file behind
        if err:
            fail += 1
            print("  %-34s %d copied, %d would not read" % (m, got, err))
        else:
            ok += 1
    print("\n%d module(s) completed, %d still short" % (ok, fail))
    if fail:
        print("Box is not serving those files. Re-run later; this resumes.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=DEFAULT_ROOT)
    ap.add_argument("--box", default=DEFAULT_BOX)
    ap.add_argument("--git", action="store_true")
    ap.add_argument("--from-box", dest="boxsync", action="store_true")
    ap.add_argument("--status", action="store_true")
    a = ap.parse_args()

    if a.git:
        do_git(a.root)
    if a.boxsync:
        do_box(a.root, a.box)

    complete, short, absent = survey(a.root, a.box)
    print("\nlocal runtime: %s" % a.root)
    print("  complete modules : %d" % len(complete))
    print("  partial          : %d" % len(short))
    print("  absent           : %d" % len(absent))
    if a.status:
        for m, have, want in short:
            print("     partial  %-32s %d/%d" % (m, have, want))
        for m, want in absent:
            print("     absent   %-32s 0/%d" % (m, want))
    print("\nPoint the annotator at it with:")
    print("  export LOWVAN_DATA_DIR=%s" % a.root)


main()
