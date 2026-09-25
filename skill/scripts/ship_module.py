#!/usr/bin/env python3
"""Put a finished module everywhere it has to go, and say what is still missing.

A module is not built when its PSSMs exist in a temp directory. It is built
when it is in the kit, in Box, in the runtime repo, and committed. Twice in
this project a module was measured, reported as finished, and found days later
to exist only under /tmp -- so this makes the last mile mechanical.

    python3 ship_module.py --module Carlavirus --workdir /tmp/lowvan_work/Carlavirus
    python3 ship_module.py --module Carlavirus --workdir ... --check

--check reports without writing. Exit status is non-zero while anything is
missing, so it can gate a loop.
"""
import argparse, glob, json, os, re, shutil, subprocess, sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--module", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--panel", help="panel directory, for measurements/")
    a = ap.parse_args()
    M, W = a.module, os.path.abspath(a.workdir)
    KIT, PROJ = os.environ["LOWVAN_KIT"], os.environ["LOWVAN_PROJ"]
    RT = os.environ["LOWVAN_DATA_DIR"]
    K, B = os.path.join(KIT, "modules", M), os.path.join(PROJ, "work", M)

    def have_kit():
        return len(glob.glob(os.path.join(K, "PSSM-Alignments", "*", "*", "*.fa")))
    def have_rt():
        try: return M in json.load(open(os.path.join(RT, "Viral_PSSM.json")))
        except Exception: return False
    #  Deliberately NOT detecting git state here. git invoked from a Python
    #  subprocess on macOS can die with "xcrun: error: unable to load libxcrun"
    #  and return 1 whatever the tree looks like, and routing through a login
    #  shell brought its own problems. The three placements below are what this
    #  script can verify honestly; the commit is left to the caller and printed
    #  as a reminder, because a check that lies is worse than no check.
    state = {"kit alignments": have_kit(), "Box work": os.path.isdir(B),
             "Box runtime JSON": have_rt(), }
    npssm = len(glob.glob(os.path.join(W, "Alignments", M, "*", "pssms", "*.pssm")))
    print("%s: %d PSSMs in the workdir" % (M, npssm))
    for k, v in state.items():
        print("   %-18s %s" % (k, v if not isinstance(v, bool) else ("yes" if v else "NO")))
    missing = [k for k, v in state.items() if not v]
    if a.check:
        if missing:
            print("\nMISSING: %s" % ", ".join(missing))
        else:
            print("\nplaced everywhere; confirm with:  git -C $LOWVAN_KIT status "
                  "--porcelain modules/%s" % M)
        return 1 if missing else 0

    SUB = ("corrected_alis", "reclustered_alis", "truncated_alis", "alis")
    src = os.path.join(W, "Alignments", M)
    n = 0
    for feat in sorted(os.listdir(src)):
        fd = os.path.join(src, feat)
        if not os.path.isdir(fd) or feat.startswith("."): continue
        dst = os.path.join(K, "PSSM-Alignments", feat, "corrected_alis")
        os.makedirs(dst, exist_ok=True)
        for p in sorted(glob.glob(os.path.join(fd, "pssms", "*.pssm"))):
            num = re.sub(r"\.pssm$", "", os.path.basename(p)).split(".")[-1]
            s = next((c for sd in SUB for c in [os.path.join(fd, sd, num + ".fa")]
                      if os.path.exists(c)), None)
            if s: shutil.copyfile(s, os.path.join(dst, num + ".fa")); n += 1
    os.makedirs(os.path.join(K, "Rep-Contigs"), exist_ok=True)
    for f in glob.glob(os.path.join(W, "Rep-Contigs", "*")):
        shutil.copy(f, os.path.join(K, "Rep-Contigs"))
    for f in glob.glob(os.path.join(W, "%s_Viral_PSSM.json" % M)):
        shutil.copy(f, K)
    os.makedirs(os.path.join(B, "measurements"), exist_ok=True)
    for sub in ("Alignments", "Rep-Contigs", "collections"):
        if os.path.isdir(os.path.join(W, sub)):
            subprocess.run(["rsync", "-a", os.path.join(W, sub) + "/",
                            os.path.join(B, sub) + "/"])
    for pat in ("*.py", "*.json", "*.sh", "*.log", "*.tsv"):
        for f in glob.glob(os.path.join(W, pat)): shutil.copy(f, B)
    if a.panel:
        for f in glob.glob(os.path.join(a.panel, "quality*.tsv")):
            shutil.copy(f, os.path.join(B, "measurements"))
    subprocess.run([sys.executable, os.path.join(KIT, "skill/scripts/install_module.py"),
                    "--workdir", W, "--repo", RT])
    print("\nshipped %d alignments to the kit; now commit and push modules/%s" % (n, M))
    return 0

sys.exit(main())
