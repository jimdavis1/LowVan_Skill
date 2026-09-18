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
import argparse, collections, glob, json, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gto-dir", required=True)
    ap.add_argument("--repo", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--report", help="write the per-genome table here")
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
        rows.append((b, nf, len(gf), len(cf), "; ".join(allf)[:200]))

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
            fh.write("genome\tn_features\tn_genome_flags\tn_contig_flags\tflags\n")
            for r in sorted(rows):
                fh.write("\t".join(str(x) for x in r) + "\n")
        print("\n  wrote %s" % args.report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
