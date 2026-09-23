#!/usr/bin/env python3
"""Annotate GTOs and score their quality -- the good-vs-poor half of step 11.

Feature tables answer "which proteins were called". They cannot answer "is this
genome well annotated", because they cannot represent special features and they
miscount fragments: a genome whose L is split across two contigs shows a short
L rather than a fragmented one, and the quality check then reports "Feature is
too short: 776" for a module that is behaving correctly.

So the quality questions run on GTOs:

    make_gto.py            fasta -> GTO, ids issued by the ID server
    run_gto_eval.py        annotate, then score, then tabulate

    python3 run_gto_eval.py --gto-dir gto --repo $LOWVAN_DATA_DIR \\
            --out gto_out --jobs 8

Reports the good/poor split and, for the poor ones, which flags fired and how
often -- which is the useful output, because it separates "the module missed
something" from "this record is a 900-base partial submission".
"""
import argparse, collections, glob, json, os, shutil, subprocess, sys

SPECIAL_DIR = {"te": "Transcript-Editing", "sv": "Splice-Variants"}
from concurrent.futures import ThreadPoolExecutor


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gto-dir", required=True)
    ap.add_argument("--repo", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--report", help="write the per-genome table here")
    ap.add_argument("--tbl-dir",
                    help="also save each genome's feature table here as "
                         "<genome>.feature.tbl. The annotator already writes one "
                         "-- annotate_by_viral_pssm-GTO.pl passes -p and -tbl, so "
                         "the 20-column table lands in _wd/<genome>/<genome>.stdout.txt "
                         "-- and it was simply being thrown away with _wd. "
                         "annotation_rarefaction.py and analyze_feature_gap.py both "
                         "need it, and hand-synthesising one is how a gap analysis "
                         "came back '0 genomes' from an empty gene-symbol column.")
    ap.add_argument("--skip-special", action="store_true",
                    help="do not call transcript-edited or spliced features. "
                         "The default is to call them: a coverage audit that "
                         "omits a module's special features understates it "
                         "without saying so.")
    ap.add_argument("--perl5lib", help="extra PERL5LIB entries, colon-separated")
    ap.add_argument("--perl", default="perl",
                    help="perl for the GTO wrapper. GenomeTypeObject needs a UUID "
                         "module (Data::UUID or UUID); the BV-BRC runtime perl has "
                         "one and a conda perl usually does not. The inner "
                         "annotate_by_viral_pssm.pl still resolves via PATH, so it "
                         "keeps the interpreter its own deps are installed under.")
    args = ap.parse_args()

    R = os.path.abspath(args.repo)
    #  Everything below joins onto args.out and args.report, and the
    #  subprocesses run with cwd=<out>/_wd, so relative values resolve
    #  against the wrong directory. --repo was already absolutised; these
    #  were not, giving 'Could not open output file gto_out/x.ann.gto' on
    #  every genome while the run still exited 0.
    args.out = os.path.abspath(args.out)
    if args.report: args.report = os.path.abspath(args.report)
    os.makedirs(args.out, exist_ok=True)
    #  the GTO script needs GenomeTypeObject.pm, and it shells out to
    #  annotate_by_viral_pssm.pl by bare name, so the repo must be on PATH.
    p5 = [os.path.expanduser("~/perl5/lib/perl5")]
    if args.perl5lib:
        p5 = args.perl5lib.split(":") + p5
    if os.environ.get("PERL5LIB"):
        p5.append(os.environ["PERL5LIB"])
    env = dict(os.environ, LOWVAN_DATA_DIR=R,
               PATH=R + ":" + os.environ["PATH"],
               PERL5LIB=":".join(p5))

    #  The annotator runs with cwd=<out>/_wd, so a relative --gto-dir
    #  resolves against the wrong directory and every genome fails with
    #  "Could not open input file gto/x.gto". --repo is already
    #  absolutised above; this was not. 669/669 failures on Hepeviridae.
    if args.tbl_dir:
        os.makedirs(args.tbl_dir, exist_ok=True)
    gtos = sorted(glob.glob(os.path.join(os.path.abspath(args.gto_dir), "*.gto")))
    print("  %d GTO(s)" % len(gtos))

    def one(g):
        b = os.path.splitext(os.path.basename(g))[0]
        ann = os.path.join(args.out, b + ".ann.gto")
        qual = os.path.join(args.out, b + ".qual.gto")
        if os.path.exists(qual) and os.path.getsize(qual) > 0:
            return b
        wd = os.path.join(args.out, "_wd", b)
        os.makedirs(wd, exist_ok=True)
        if not (os.path.exists(ann) and os.path.getsize(ann) > 0):
            r = subprocess.run([args.perl, os.path.join(R, "annotate_by_viral_pssm-GTO.pl"),
                                "-i", g, "-o", ann, "-t", "1", "-x", b],
                               capture_output=True, text=True, env=env, cwd=wd)
            if not os.path.exists(ann) or os.path.getsize(ann) == 0:
                with open(os.path.join(args.out, b + ".FAIL"), "w") as fh:
                    fh.write((r.stderr or r.stdout or "")[-3000:])
                return b
        #  Special features BEFORE quality scoring, so the quality score sees
        #  them. This used not to happen at all, and the consequence was
        #  quiet: every coverage audit in the report series measured only the
        #  PSSM features of its module, with no indication that anything was
        #  missing. Togaviridae TF, the coronavirus ORF1ab and NSP12, the
        #  paramyxovirus V/W pairs and Orthoflavivirus NS1' were all absent
        #  from their own numbers. A module is not measured until its
        #  declared features are measured, whichever program calls them.
        #
        #  Each pass detects its module from the GTO and is a no-op when that
        #  module declares nothing of its kind, so running both costs little
        #  on modules without special features. --skip-special opts out.
        src = ann
        if not args.skip_special:
            for prog, kind in (("get_transcript_edited_features.pl", "te"),
                               ("get_splice_variant_features.pl", "sv")):
                if not os.path.isdir(os.path.join(R, SPECIAL_DIR[kind])):
                    continue
                dst = os.path.join(wd, b + "." + kind + ".gto")
                rs = subprocess.run([args.perl, os.path.join(R, prog),
                                     "-i", src, "-o", dst,
                                     "-j", os.path.join(R, "Viral_PSSM.json"),
                                     "-d", os.path.join(R, SPECIAL_DIR[kind]),
                                     "-a", str(max(1, args.jobs // 2))],
                                    capture_output=True, text=True, env=env, cwd=wd)
                if os.path.exists(dst) and os.path.getsize(dst) > 0:
                    src = dst
                else:
                    with open(os.path.join(args.out, b + ".SFAIL." + kind), "w") as fh:
                        fh.write((rs.stderr or rs.stdout or "")[-3000:])
            #  the .ann.gto downstream tools read should carry everything
            if src != ann:
                shutil.copyfile(src, ann)
        if args.tbl_dir:
            src_tbl = os.path.join(wd, b + ".stdout.txt")
            if os.path.exists(src_tbl) and os.path.getsize(src_tbl) > 0:
                shutil.copyfile(src_tbl, os.path.join(args.tbl_dir, b + ".feature.tbl"))
        r = subprocess.run([args.perl, os.path.join(R, "viral_genome_quality.pl"),
                            "-i", ann, "-o", qual, "-p", b],
                           capture_output=True, text=True, env=env, cwd=wd)
        if not os.path.exists(qual) or os.path.getsize(qual) == 0:
            with open(os.path.join(args.out, b + ".QFAIL"), "w") as fh:
                fh.write((r.stderr or r.stdout or "")[-3000:])
        return b

    done = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for _ in pool.map(one, gtos):
            done += 1
            if done % 100 == 0:
                print("    %d/%d" % (done, len(gtos)), flush=True)

    #  tabulate
    good = poor = unann = featflagged = 0
    flags = collections.Counter()
    ffl = collections.Counter()
    rows = []
    for q in sorted(glob.glob(os.path.join(args.out, "*.qual.gto"))):
        b = os.path.basename(q)[:-9]
        try:
            d = json.load(open(q))
        except Exception:
            continue
        gf = d.get("genome_quality_flags") or []
        cf = [f for c in d.get("contigs", []) for f in (c.get("contig_quality_flags") or [])]
        feat = d.get("features") or []
        nf = len(feat)
        for f in feat:
            for x in (f.get("feature_quality_flags") or []):
                ffl[x.split(":")[0].strip()] += 1
        allf = gf + cf
        for x in allf:
            flags[x.split(";")[0].split(":")[0].strip()] += 1
        hasff = any(f.get("feature_quality_flags") for f in feat)
        if nf == 0:
            unann += 1
        elif allf:
            poor += 1
        else:
            good += 1
            if hasff: featflagged += 1
        #  The feature-flag count belongs in the table. "Genuinely clean" --
        #  the figure every coverage audit leads with -- is no genome, contig
        #  OR feature flag, and without this column the table cannot reproduce
        #  it: reading quality.tsv gives the "good" number instead, which is
        #  several points higher and answers a different question.
        nff = sum(len(f.get("feature_quality_flags") or []) for f in feat)
        rows.append((b, nf, len(gf), len(cf), nff, "; ".join(allf)[:200]))

    tot = good + poor + unann
    print("\n  scored %d genome(s)" % tot)
    if tot:
        #  The verdict is genome+contig flags only; feature-level flags are
        #  reported but do not disqualify a genome. Say so in the label. A
        #  Trivirinae genome with its replicase called 2,254 aa
        #  against a max of 2,195 and its movement protein 563 against 498
        #  counts as "good", and "good (no flags)" made that read as a genome
        #  with nothing wrong with it.
        print("    good (no genome/contig flags) %5d  %5.1f%%" % (good, 100.0*good/tot))
        print("    poor (genome/contig flagged)  %5d  %5.1f%%" % (poor, 100.0*poor/tot))
        if featflagged:
            print("      of which %d carry feature-level flags -- counted below,"
                  " but not against the verdict" % featflagged)
            print("    genuinely clean               %5d  %5.1f%%"
                  % (good - featflagged, 100.0 * (good - featflagged) / tot))
        print("    no features called    %5d  %5.1f%%" % (unann, 100.0*unann/tot))
    print("\n  genome/contig flags, most common first:")
    for k, v in flags.most_common(15):
        print("    %5d  %s" % (v, k))
    if ffl:
        print("\n  feature-level flags:")
        for k, v in ffl.most_common(10):
            print("    %5d  %s" % (v, k))
    fails = len(glob.glob(os.path.join(args.out, "*.FAIL")))
    qfails = len(glob.glob(os.path.join(args.out, "*.QFAIL")))
    if fails or qfails:
        print("\n  %d annotation failure(s), %d quality failure(s)" % (fails, qfails))

    if args.report:
        with open(args.report, "w") as fh:
            fh.write("genome\tn_features\tn_genome_flags\tn_contig_flags"
                     "\tn_feature_flags\tflags\n")
            for r in sorted(rows):
                fh.write("\t".join(str(x) for x in r) + "\n")
        print("\n  wrote %s" % args.report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
