#!/usr/bin/env python3
"""Collect profiles built from several collections into ONE feature.

    python3 merge_feature_profiles.py --workdir . --module Orthohepadnavirus \
        --into P17 --from P17_NONA --from P17_A

WHY A FEATURE WOULD EVER BE BUILT IN PIECES.  Normally one collection is one
feature.  The exception is a protein whose ENDPOINT is decided by something
the clustering cannot see, so that a single alignment of all its members has a
ragged terminus and no profile can place the cut.

Orthohepadnavirus p17 (HBeAg), 24 Nov 2026.  p17 is p22 minus the C-terminal
arginine-rich domain, removed by furin after an RXXR motif -- and WHICH motif
is genotype-dependent.  Non-A genotypes cut after the first, core
`151-RRGR-154`, giving 164 aa.  Genotype A carries a two-residue insertion
that makes the first motif `151-RRDR-154`, a poor furin substrate, so the cut
shifts to the second, `156-RSPR-159`, giving 169 aa.

Clustered together at 97% identity those two populations merge, the alignment
frays over the last five columns, and the profile calls an endpoint that is
wrong for one genotype or the other.  Split by motif before clustering, they
come out clean -- measured here as 2,986 sequences at exactly 164 aa and 419
at exactly 169, with the alignments ending in `VRRRGR` and `RGRSPR`.

WHY NOT SHIP THEM AS TWO FEATURES.  Because at ~97% identity and five residues
apart they cross-call: each profile scores well on the other's genomes, and
the annotator emits both, or the wrong one.  That is the HBwx lesson in a
different costume -- a PSSM cannot veto on a handful of columns.

As two profiles under ONE feature the ambiguity disappears: the annotator
scores both and keeps the better, which is that genome's genotype and
therefore that genome's correct C-terminus.  One feature, one call, one
endpoint, and the genotype decided by score instead of by a rule.

WHAT THIS DOES.  Copies each source's profiles and their corrected alignments
into <into>/, renumbered from 1 in source order, and writes MERGED_FROM
recording where each number came from.  The sources are left alone: they stay
on disk as the build record, and their collections stay in collections/ so
both audits still see them.  Re-running replaces the merged directory.

The merged feature must appear in features.json and the sources must NOT,
or they ship as features in their own right and cross-call exactly as above.
"""
import argparse
import glob
import json
import os
import re
import shutil
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--module", required=True)
    ap.add_argument("--into", required=True,
                    help="the feature name that ships, e.g. P17")
    ap.add_argument("--from", dest="sources", action="append", required=True,
                    help="a built feature directory to take profiles from; "
                         "repeat, and the order sets the numbering")
    ap.add_argument("--features", help="features.json, checked for the "
                    "merged name present and the sources absent")
    a = ap.parse_args()

    base = os.path.join(a.workdir, "Alignments", a.module)
    if not os.path.isdir(base):
        raise SystemExit("no Alignments/%s under %s" % (a.module, a.workdir))

    #  Refuse before writing anything if a source was never built, rather
    #  than producing a merged feature that is quietly missing a genotype.
    pairs = []
    for src in a.sources:
        sd = os.path.join(base, src)
        pss = sorted(glob.glob(os.path.join(sd, "pssms", "*.pssm")),
                     key=lambda p: int(re.sub(r"\.pssm$", "", p).split(".")[-1]))
        if not pss:
            raise SystemExit("%s has no profiles in %s/pssms -- build it first"
                             % (src, sd))
        for p in pss:
            num = re.sub(r"\.pssm$", "", os.path.basename(p)).split(".")[-1]
            fa = os.path.join(sd, "corrected_alis", num + ".fa")
            if not os.path.exists(fa):
                fa = os.path.join(sd, "alis", num + ".fa")
            if not os.path.exists(fa):
                raise SystemExit("%s profile %s has no alignment beside it"
                                 % (src, num))
            pairs.append((src, num, p, fa))

    dst = os.path.join(base, a.into)
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    os.makedirs(os.path.join(dst, "pssms"))
    os.makedirs(os.path.join(dst, "corrected_alis"))

    record = []
    for i, (src, num, p, fa) in enumerate(pairs, 1):
        shutil.copy2(p, os.path.join(dst, "pssms", "%s.%s.%d.pssm"
                                     % (a.module, a.into, i)))
        shutil.copy2(fa, os.path.join(dst, "corrected_alis", "%d.fa" % i))
        n = sum(1 for line in open(fa) if line.startswith(">"))
        record.append({"profile": i, "from": src, "was": num, "seqs": n})
        print("  %-3d <- %s/%s  (%d sequences)" % (i, src, num, n))

    with open(os.path.join(dst, "MERGED_FROM"), "w") as o:
        o.write("# %s is built from several collections; see "
                "merge_feature_profiles.py\n" % a.into)
        json.dump(record, o, indent=1)
        o.write("\n")

    print("\n%s: %d profile(s) from %s"
          % (a.into, len(pairs), ", ".join(a.sources)))

    if a.features and os.path.exists(a.features):
        f = json.load(open(a.features))
        f = f.get("features", f)
        bad = [s for s in a.sources if s in f]
        if a.into not in f:
            print("\nWARNING: %s is not in %s -- it will not ship"
                  % (a.into, a.features), file=sys.stderr)
        if bad:
            print("\nWARNING: %s still in %s. The sources must NOT be "
                  "features or they cross-call the merged one."
                  % (", ".join(bad), a.features), file=sys.stderr)
    return 0


sys.exit(main())
