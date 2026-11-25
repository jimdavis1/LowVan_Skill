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
    #  THE CURATED PMID FIELD IS NOT OURS TO WRITE. A bare `PMID` means a human
    #  read the paper and confirmed it supports the claim; everything a model
    #  produced goes in `PMID_claude_generated`. SKILL.md has said so in three
    #  places since the kit existed and 40 model-written citations still ended
    #  up in `PMID` across six modules. So it is checked here, where a module
    #  cannot be shipped past it, rather than only written down. The audit
    #  compares against git HEAD, so a citation Jim curated and committed is
    #  left alone and only a NEW one fails.
    pm = subprocess.run([sys.executable,
                         os.path.join(KIT, "skill/scripts/audit_pmids.py"),
                         "--module", M, "--strict"],
                        capture_output=True, text=True)
    if pm.returncode != 0:
        print(pm.stdout)
        raise SystemExit(
            "REFUSED: %s carries a citation in the curated `PMID` field that is "
            "not in git HEAD.\nMove it with:  python3 $SK/audit_pmids.py "
            "--module %s --fix" % (M, M))

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
    #  Only the module's OWN features ship. A feature built from several
    #  collections keeps its sources on disk as the build record -- the
    #  Orthohepadnavirus p17 genotype split is merged into P17 by
    #  merge_feature_profiles.py while P17_A and P17_NONA stay beside it --
    #  and those directories are not features. install_module.py applies the
    #  same filter; without it here the kit carries alignments for profiles
    #  no shipped feature references.
    try:
        keys = set(json.load(open(os.path.join(
            W, "%s_Viral_PSSM.json" % M)))[M].get("features", {}))
    except Exception:
        keys = set()
    for feat in sorted(os.listdir(src)):
        fd = os.path.join(src, feat)
        if not os.path.isdir(fd) or feat.startswith("."): continue
        if keys and feat not in keys:
            print("   %-18s %s (not a feature in the module JSON)"
                  % ("not shipped", feat))
            continue
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
    #  features_all.json GOES TO THE KIT TOO, because it is the only place the
    #  curation record lives. Its _note fields carry why each copy_num was set
    #  or dropped, which non-callers were checked and what tblastn returned on
    #  them, and why a feature is spliced rather than profiled. The 7b9c6f08
    #  commit message says the Spumaretrovirinae POL and TAS misses are
    #  "logged open ... in the feature notes" -- and they are, in this file,
    #  which until now was copied only into the Box work directory and was
    #  tracked by nothing. The repo held the profiles and the pages and not
    #  one line of the reasoning behind them.
    for f in glob.glob(os.path.join(W, "features_all.json")):
        shutil.copy(f, K)
    #  THE AUDITS AND THE PROTECTED PAIRS GO TO THE KIT, for the same reason
    #  features_all.json does: they are the evidence the module is sound, and
    #  nothing else in the repo holds them.
    #
    #  build_features.py refuses to build until LENGTH_AUDIT and BIN_AUDIT
    #  both read `Unresolved: 0`, so every shipped module HAS a clean pair of
    #  reports -- and until now they stayed in the work directory, which is a
    #  temp dir on one machine. A reviewer could see the profiles and not the
    #  measurement that licensed them. Orthohepadnavirus LENGTH_AUDIT.md was
    #  in the repo only because it was copied by hand.
    #
    #  protected_pairs.json matters more than it looks. A protected pair tells
    #  the bin audit NOT to judge two collections that are near-identical by
    #  design -- the nested isoform families, L/M/S and preC/C here. That
    #  suppresses findings, so it is a curation decision of exactly the kind
    #  the repo should carry, and without it the audit is not reproducible:
    #  re-running it elsewhere reports defects that were deliberately waived.
    for cands in ([os.path.join(W, "LENGTH_AUDIT.%s.md" % M),
                   os.path.join(W, M, "LENGTH_AUDIT.md"),
                   os.path.join(W, "LENGTH_AUDIT.md")],
                  [os.path.join(W, "BIN_AUDIT.%s.md" % M),
                   os.path.join(W, M, "BIN_AUDIT.md"),
                   os.path.join(W, "BIN_AUDIT.md")],
                  [os.path.join(W, "collections", M, "protected_pairs.json"),
                   os.path.join(W, M, "collections", "protected_pairs.json")]):
        src = next((p for p in cands if os.path.exists(p)), None)
        if src:
            dst = os.path.join(K, os.path.basename(cands[-1]
                               if src.endswith(".json") else cands[1]))
            shutil.copy2(src, dst)
            print("   %-18s %s" % ("curation record", os.path.basename(dst)))
    #  SPLICE AND TRANSCRIPT-EDIT REFERENCES GO TO THE REPO TOO.
    #  install_module.py puts them in the runtime and nothing put them in the
    #  kit, so a module declaring a special feature could not be reinstalled
    #  from the repository alone -- the profiles and the module JSON were
    #  there and the references the feature is actually called from were not.
    #  Betaretrovirus REM.fasta and Spumaretrovirinae BET.fasta are in the
    #  repo only because they were placed by hand. They are curated artifacts,
    #  not build output: Rem's were mapped from NCBI because the dump has no
    #  Rem at all, and Deltaretrovirus REX took five lineage queries and a
    #  de-duplication pass to assemble.
    for kind in ("Splice-Variants", "Transcript-Editing"):
        for cand in (os.path.join(W, kind, M), os.path.join(W, kind)):
            fs = glob.glob(os.path.join(cand, "*.fasta"))
            if not fs:
                continue
            dd = os.path.join(KIT, kind, M)
            os.makedirs(dd, exist_ok=True)
            for f in sorted(fs):
                shutil.copy2(f, os.path.join(dd, os.path.basename(f)))
            print("   %-18s %d file(s) -> %s/%s/" % (kind, len(fs), kind, M))
            break
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
