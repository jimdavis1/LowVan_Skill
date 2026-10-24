#!/usr/bin/env python3
"""Stage an INSTALLED module out of a Viral_Annotation checkout into the kit's modules/.

`ship_module.py` stages from the build workdir, which is wrong twice over.

1. A workdir accumulates profiles the module no longer declares. HIV1 carries
   885 `.pssm` files there against the 763 that are installed, so staging from
   it ships 122 alignments for profiles that do not exist, and the kit then
   rebuilds a module the annotator never ran.
2. It never copied `Splice-Variants/` or `Transcript-Editing/`. Those were put
   in by hand, which worked for Merhavirus and Togaviridae and would silently
   not happen for the next module whose author did not know to do it.

The installed repo is the pruned, authoritative state -- it is what the
annotator actually loads -- so that is what this reads. Verified against
Norovirus: repo `PSSM-Alignments/<M>/<FEAT>/N.fa` is byte-identical to kit
`modules/<M>/PSSM-Alignments/<FEAT>/corrected_alis/N.fa`.

Layout written, which is what `install_module.py` and the modules/ README
expect:

    modules/<M>/<M>_Viral_PSSM.json              {"<M>": block}, json indent=1
    modules/<M>/Rep-Contigs/<M>.N.dna            routing references
    modules/<M>/Rep-Contigs/close_genomes.json   lifted out of the block
    modules/<M>/PSSM-Alignments/<FEAT>/corrected_alis/N.fa
    modules/<M>/Splice-Variants/<FEAT>.fasta     if any feature is special:splice
    modules/<M>/Transcript-Editing/<FEAT>.fasta  if special:transcript_edit

Every copy is size-verified. Box Drive returns short or empty reads on
un-cached content and exits 0 while doing it, so a copy that is not checked is
not a copy.

    python3 stage_module.py --repo <VA checkout> --kit <LowVan_Skill> \
            --module HIV1 [--module ...] [--check]

--check reports and writes nothing. Exit status is non-zero if anything is
missing or inconsistent, so it can gate a loop.
"""
import argparse
import json
import os
import shutil
import sys

PAGE_SUBDIR = "corrected_alis"


def copy_verified(src, dst):
    """Copy and confirm the bytes arrived. Returns the size written."""
    want = os.path.getsize(src)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copyfile(src, dst)
    got = os.path.getsize(dst)
    if got != want:
        raise IOError("short copy: %s is %d bytes, %s is %d" % (src, want, dst, got))
    if want == 0:
        raise IOError("source is empty, which on Box means unreadable: %s" % src)
    return got


def stage(repo, kit, module, check=False):
    master = json.load(open(os.path.join(repo, "Viral_PSSM.json")))
    if module not in master:
        print("  %-24s NOT INSTALLED in %s" % (module, repo))
        return False
    block = master[module]
    feats = block.get("features", {})
    dest = os.path.join(kit, "modules", module)

    repo_alis = os.path.join(repo, "PSSM-Alignments", module)
    repo_pssms = os.path.join(repo, "Viral-PSSMs", "%s.pssms" % module)

    #  Which features carry profiles, and which are called by other machinery.
    #  A splice or transcript_edit feature has reference fastas instead, and a
    #  declared feature with neither is the defect install_module.py refuses.
    special = {f: v.get("special") for f, v in feats.items()}
    pssm_feats = [f for f, s in special.items() if s not in ("splice", "transcript_edit")]

    problems, plan = [], []

    n_ali = n_pssm = 0
    for feat in sorted(pssm_feats):
        fd = os.path.join(repo_alis, feat)
        if not os.path.isdir(fd):
            problems.append("%s: no alignments in the repo" % feat)
            continue
        alis = sorted(f for f in os.listdir(fd) if f.endswith(".fa"))
        pd = os.path.join(repo_pssms, feat)
        npr = len([f for f in os.listdir(pd) if f.endswith(".pssm")]) if os.path.isdir(pd) else 0
        if not alis:
            problems.append("%s: declared but has no alignment" % feat)
        if npr != len(alis):
            problems.append("%s: %d alignments against %d profiles" % (feat, len(alis), npr))
        n_ali += len(alis)
        n_pssm += npr
        for a in alis:
            plan.append((os.path.join(fd, a),
                         os.path.join(dest, "PSSM-Alignments", feat, PAGE_SUBDIR, a)))

    #  Routing references.
    reps = sorted(f for f in os.listdir(os.path.join(repo, "Viral-Rep-Contigs"))
                  if f.startswith(module + ".") and f.endswith(".dna"))
    if not reps:
        problems.append("no rep contigs")
    for r in reps:
        plan.append((os.path.join(repo, "Viral-Rep-Contigs", r),
                     os.path.join(dest, "Rep-Contigs", r)))

    #  Reference fastas for features no PSSM can call.
    for kind, sub in (("splice", "Splice-Variants"), ("transcript_edit", "Transcript-Editing")):
        want = [f for f, s in special.items() if s == kind]
        if not want:
            continue
        srcdir = os.path.join(repo, sub, module)
        for feat in sorted(want):
            cand = os.path.join(srcdir, "%s.fasta" % feat)
            if os.path.exists(cand):
                plan.append((cand, os.path.join(dest, sub, "%s.fasta" % feat)))
            else:
                problems.append("%s: declared special:%s with no %s/%s.fasta"
                                % (feat, kind, sub, feat))

    print("  %-24s %2d features, %d alignments / %d profiles, %d refs%s"
          % (module, len(feats), n_ali, n_pssm, len(reps),
             "" if not problems else "   <-- %d PROBLEM(S)" % len(problems)))
    for p in problems:
        print("       ! %s" % p)
    if check or problems:
        return not problems

    #  Write. Replace the alignment and reference trees wholesale so a feature
    #  that was dropped since the last staging does not survive as a leftover --
    #  the stale-alignment defect install_module.py had.
    for sub in ("PSSM-Alignments", "Rep-Contigs", "Splice-Variants", "Transcript-Editing"):
        d = os.path.join(dest, sub)
        if os.path.isdir(d):
            shutil.rmtree(d)
    total = 0
    for src, dst in plan:
        total += copy_verified(src, dst)

    #  close_genomes is also written beside the contigs it describes, because
    #  that is where a curator looks to see what a reference actually is. It
    #  stays in the block as well -- install_module.py reads it from there, and
    #  every module already staged carries it twice.
    cg = block.get("close_genomes")
    if cg is not None:
        with open(os.path.join(dest, "Rep-Contigs", "close_genomes.json"), "w") as fh:
            fh.write(json.dumps(cg, indent=1, sort_keys=True))
    with open(os.path.join(dest, "%s_Viral_PSSM.json" % module), "w") as fh:
        fh.write(json.dumps({module: block}, indent=1, sort_keys=True))

    print("       staged %d files, %.1f MB" % (len(plan), total / 1e6))
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True, help="Viral_Annotation checkout to read")
    ap.add_argument("--kit", required=True, help="LowVan_Skill checkout to write")
    ap.add_argument("--module", action="append", required=True)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    print("%s %d module(s) from %s"
          % ("checking" if a.check else "staging", len(a.module), a.repo))
    ok = [stage(a.repo, a.kit, m, a.check) for m in a.module]
    bad = len([x for x in ok if not x])
    print("\n%d ok, %d with problems" % (len(ok) - bad, bad))
    return 1 if bad else 0


sys.exit(main())
