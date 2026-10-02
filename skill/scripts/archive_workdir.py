#!/usr/bin/env python3
"""Preserve a module's build tree out of scratch, and check that it arrived.

Builds run in scratch because Box sync turns many-small-file I/O into load
average with the job itself at single-digit CPU. That is the right call for
speed and the wrong one for memory: scratch is `/private/tmp`, a reboot takes
it, and the build tree is the only record of HOW a module was made --
`clusters/` (what grouped with what), `alis/` before correction and
`corrected_alis/` after, `Curation_Report`, `Truncation_Report`,
`FLAGGED_NTERM`, `BUILD_PARAMS`, `Leftover_Seqs.aa`. None of it is
reconstructible from the shipped module, and a curator cannot check the work
without it.

**You say where it goes.** There is no default: `--dest` is required, because
the archive is the user's storage decision, not this script's.

    python3 archive_workdir.py --workdir $LOWVAN_SCRATCH --module Enterovirus \\
            --dest "$LOWVAN_PROJ/work"

That writes <dest>/<Module>/ holding the whole tree, matching how every other
taxon stores it:

    <Taxon>/work/<Module>/Alignments/<Module>/<FEATURE>/{alis, clusters,
        corrected_alis, reclustered_alis, Curation_Report, Truncation_Report,
        FLAGGED_NTERM, BUILD_PARAMS, Leftover_Seqs.aa, pssms, ...}

**The archive is local, not a git artifact.** The repo carries the module --
JSON, rep contigs, and the alignments behind the shipped profiles. The build
tree stays wherever --dest puts it; `.gitignore` excludes `Alignments/` and
`collections/` so it cannot be committed by accident.

With --repo it also checks the one thing that is easy to get wrong: that the
alignments preserved for each feature account for every PSSM that was
installed. `corrected_alis/` alone does NOT -- the leftover pass and the
N-terminal re-clustering write to `reclustered_alis/`, and on Enterovirus that
was 67 of 701 profiles, 32 of them in 3D alone, with nothing reporting the
shortfall.
"""
import argparse, os, subprocess, sys, glob

ap = argparse.ArgumentParser()
ap.add_argument("--workdir", required=True, help="the scratch build directory")
ap.add_argument("--module", required=True)
ap.add_argument("--dest", required=True,
                help="where the build tree is kept. Required on purpose: this "
                     "is your storage decision, not a default")
ap.add_argument("--repo", default=None,
                help="Viral_Annotation checkout; enables the installed-PSSM "
                     "coverage check")
ap.add_argument("--dry-run", action="store_true")
a = ap.parse_args()

if not os.path.isdir(a.workdir): sys.exit("no such workdir: %s" % a.workdir)
dest = os.path.join(a.dest, a.module)

#  tmp/ is scratch within scratch; everything else is evidence
EXCLUDE = ["tmp/", "*.blastdb", "__pycache__/"]
cmd = ["rsync", "-a"] + ["--exclude=%s" % e for e in EXCLUDE]
if a.dry_run: cmd.append("--dry-run")
cmd += [a.workdir.rstrip("/") + "/", dest + "/"]

if not a.dry_run: os.makedirs(dest, exist_ok=True)
print("archiving %s -> %s" % (a.workdir, dest))
r = subprocess.run(cmd)
if r.returncode: sys.exit("rsync failed (%d)" % r.returncode)
if a.dry_run: sys.exit(0)

adir = os.path.join(dest, "Alignments", a.module)
if not os.path.isdir(adir):
    sys.exit("ARCHIVED, BUT no Alignments/%s in the result -- check --workdir" % a.module)

feats = sorted(d for d in os.listdir(adir) if os.path.isdir(os.path.join(adir, d)))
print("\n  %-14s %8s %8s %8s %8s" % ("feature", "alis", "corrected", "reclust", "reports"))
problems = []
for f in feats:
    d = os.path.join(adir, f)
    n = lambda sub: len(glob.glob(os.path.join(d, sub, "*.fa")))
    reports = [r for r in ("Curation_Report", "Truncation_Report", "BUILD_PARAMS",
                           "FLAGGED_NTERM", "Leftover_Seqs.aa")
               if os.path.exists(os.path.join(d, r))]
    print("  %-14s %8d %8d %8d %8d" % (f, n("alis"), n("corrected_alis"),
                                       n("reclustered_alis"), len(reports)))
    if not n("corrected_alis"): problems.append("%s: no corrected_alis" % f)
    if not os.path.exists(os.path.join(d, "BUILD_PARAMS")):
        problems.append("%s: no BUILD_PARAMS -- the departures from the "
                        "documented defaults are unrecorded" % f)

if a.repo:
    print("\n  checking that every installed PSSM has its alignment preserved:")
    short = 0
    for f in feats:
        inst = glob.glob(os.path.join(a.repo, "Viral-PSSMs", a.module + ".pssms",
                                      f, "*.pssm"))
        have = set()
        for sub in ("corrected_alis", "reclustered_alis"):
            have |= {os.path.basename(p)[:-3]
                     for p in glob.glob(os.path.join(adir, f, sub, "*.fa"))}
        #  installed profiles are named <Module>.<FEAT>.<n>.pssm; the
        #  alignment they came from is <n>.fa. Compare the cluster id, not
        #  the filename.
        pre = "%s.%s." % (a.module, f)
        ids = set()
        for p in inst:
            b = os.path.basename(p)[:-5]
            ids.add(b[len(pre):] if b.startswith(pre) else b)
        miss = ids - have
        if miss:
            short += len(miss)
            print("    %-14s %d of %d installed profile(s) have no preserved "
                  "alignment: %s" % (f, len(miss), len(inst),
                                     ", ".join(sorted(miss)[:6])))
            problems.append("%s: %d profile(s) with no alignment" % (f, len(miss)))
    print("    %s" % ("every installed profile has its alignment" if not short
                      else "%d profile(s) unaccounted for" % short))

print()
if problems:
    print("ARCHIVED WITH %d PROBLEM(S):" % len(problems))
    for p in problems: print("   %s" % p)
    sys.exit(1)
print("archived cleanly: %d feature(s) under %s" % (len(feats), adir))
