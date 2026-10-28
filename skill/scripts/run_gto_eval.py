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


def _special_functions(repo):
    """Every `anno` belonging to a feature a special pass calls, any module."""
    try:
        j = json.load(open(os.path.join(repo, "Viral_PSSM.json")))
    except Exception:
        return set()
    out = set()
    for blk in j.values():
        if not isinstance(blk, dict):
            continue
        for f in (blk.get("features") or {}).values():
            if isinstance(f, dict) and f.get("special") and f.get("anno"):
                out.add(f["anno"])
    return out


def _has_special(ann_gto, repo):
    sp = _special_functions(repo)
    if not sp:
        return False
    try:
        g = json.load(open(ann_gto))
    except Exception:
        return False
    return any(f.get("function") in sp for f in g.get("features", []))

def _append_special_rows(tbl, ann_gto):
    """Add the spliced / transcript-edited features to the feature table.

    The table is annotate_by_viral_pssm-GTO.pl's own stdout, and that program
    runs BEFORE get_splice_variant_features.pl -- so the table it writes can
    never contain a spliced feature, while the .ann.gto beside it does. The
    two then disagree about the same genome, and every table-driven analysis
    silently loses the special features: measured on HIV-1, 26 coverage
    tables carried zero Tat and zero Rev rows while their GTOs carried 22 and
    24. A vocabulary rarefaction built from those tables would have reported
    the module's controlled vocabulary as 7 strings instead of 9, and nothing
    in the output would have looked wrong.

    Reconciles by FUNCTION COUNT rather than by feature id: the table and the
    GTO number their features independently, but a function the annotator
    could not produce is missing from the table by exactly the number of
    times the special pass added it. Columns follow the 20-column layout, with
    the function in column 18 and the translation in column 20, which is what
    annotation_rarefaction.py and analyze_feature_gap.py read.
    """
    import collections as _c
    try:
        g = json.load(open(ann_gto))
    except Exception:
        return
    want = _c.Counter(f.get("function", "") for f in g.get("features", []))
    want.pop("", None)
    have = _c.Counter()
    rows = []
    for line in open(tbl, errors="replace"):
        rows.append(line)
        c = line.rstrip("\n").split("\t")
        if len(c) >= 18 and c[17]:
            have[c[17]] += 1
    missing = {k: want[k] - have.get(k, 0) for k in want if want[k] > have.get(k, 0)}
    if not missing:
        return
    #  Identity columns come from the rows the annotator already wrote, not
    #  from the GTO. The two disagree: on the HIV-1 coverage set a table
    #  named for genome 11676.1004069 held a GTO whose id was 11676.1248678
    #  and annotator rows keyed 11676.2c8645 -- three identities for one
    #  genome. Appending the GTO's id would leave two different genome ids in
    #  column 1 of a single table and silently split any consumer that groups
    #  by it. Fall back to the GTO only for a table with no rows to copy.
    gid = g.get("id", "")
    sci = g.get("scientific_name", "")
    for line in rows:
        c = line.rstrip("\n").split("\t")
        if len(c) >= 18 and c[17] and c[0]:
            gid = c[0]
            if len(c) >= 16 and c[15]: sci = c[15]
            break
    contig = ""
    for f in g.get("features", []):
        loc = f.get("location") or []
        if loc:
            contig = loc[0][0]
            break
    out = []
    for f in g.get("features", []):
        fn = f.get("function", "")
        if missing.get(fn, 0) <= 0:
            continue
        missing[fn] -= 1
        loc = f.get("location") or [[contig, "", "+", ""]]
        c = [""] * 20
        c[0] = gid
        c[1] = gid
        c[2] = loc[0][0]
        c[3] = "LV"
        c[4] = "CDS"
        c[5] = f.get("id", "")
        c[6] = (f.get("gene") or fn.split()[0].lower())
        c[7] = str(loc[0][1])
        c[9] = str(loc[0][2])
        c[15] = sci
        c[17] = fn
        c[19] = f.get("protein_translation", "")
        out.append("\t".join(c))
    if out:
        with open(tbl, "a") as fh:
            if rows and not rows[-1].endswith("\n"):
                fh.write("\n")
            fh.write("\n".join(out) + "\n")

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

    #  Refuse to start on a box that is already annotating. Two concurrent
    #  evaluations on ten cores drove the load average to 148 and throughput
    #  to nearly zero: each annotator spawns ~10 short-lived BLAST children,
    #  so 12 jobs is ~120 runnable processes. This was written down as a
    #  lesson and then re-tripped three times in one session, so it is a
    #  check now rather than a paragraph.
    try:
        import subprocess as _sp
        _n = _sp.run(["pgrep", "-fc", "annotate_by_viral_pssm.pl"],
                     capture_output=True, text=True).stdout.strip()
        _n = int(_n) if _n.isdigit() else 0
    except Exception:
        _n = 0
    if _n and not os.environ.get("LOWVAN_ALLOW_CONCURRENT_EVAL"):
        sys.exit("refusing to start: %d annotate_by_viral_pssm.pl already running.\n"
                 "Two evaluations at once thrash the box -- finish or stop the other\n"
                 "first, or set LOWVAN_ALLOW_CONCURRENT_EVAL=1 if you are certain." % _n)


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
            #  A SKIPPED genome still gets its table reconciled against its
            #  GTO. The skip is what makes a resumed run cheap, but it also
            #  freezes whatever table the earlier run left behind -- and a
            #  table written before _append_special_rows existed carries no
            #  spliced rows at all. Measured on the HIV-1 coverage set after
            #  the 28 Oct resume: 269 of 500 tables held zero Tat and zero
            #  Rev rows while their GTOs held 421 of them, and the 231 tables
            #  this run wrote agreed with their GTOs exactly. Nothing warned.
            #  The appender reconciles by function COUNT, so it is a no-op on
            #  a table that is already complete and cannot double-append.
            if args.tbl_dir:
                t = os.path.join(args.tbl_dir, b + ".feature.tbl")
                if os.path.exists(t) and os.path.getsize(t) > 0 \
                   and os.path.exists(ann) and os.path.getsize(ann) > 0:
                    _append_special_rows(t, ann)
            return b
        wd = os.path.join(args.out, "_wd", b)
        os.makedirs(wd, exist_ok=True)
        reused = os.path.exists(ann) and os.path.getsize(ann) > 0
        if not reused:
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
        #  A REUSED .ann.gto has already been through the special passes, and
        #  those passes APPEND -- they do not replace. Re-running them on it
        #  calls every spliced feature a second time on identical
        #  coordinates. Measured on Spumaretrovirinae: deleting only the
        #  *.qual.gto to re-score after a copy_num change (the documented
        #  cheap path, since copy_num changes scoring and not annotation)
        #  left 18 of 45 genomes carrying BET twice, same start, same end,
        #  same 487-residue translation, under two feature ids. The call RATE
        #  was unaffected, so nothing in the coverage report looked wrong --
        #  it surfaced only because qc_duplicate_calls.py was run afterwards.
        #  Skipped only when the reused file already carries one of the
        #  module's special functions: a genome where none was called still
        #  gets its chance, and --skip-special on the earlier run is not
        #  silently inherited by this one.
        if reused and _has_special(ann, R):
            pass
        elif not args.skip_special:
            for prog, kind in (("get_transcript_edited_features.pl", "te"),
                               ("get_splice_variant_features.pl", "sv")):
                if not os.path.isdir(os.path.join(R, SPECIAL_DIR[kind])):
                    continue
                dst = os.path.join(wd, b + "." + kind + ".gto")
                #  annotate_by_viral_pssm.pl is run with -tmp and removes its
                #  own working directory on the way out, so by the time the
                #  special-feature passes run, cwd=wd may no longer exist and
                #  subprocess dies with FileNotFoundError on the CWD rather
                #  than on the program -- which reads like a missing script.
                #  Intermittent: 21 genomes succeeded before the first hit.
                os.makedirs(wd, exist_ok=True)
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
                dst_tbl = os.path.join(args.tbl_dir, b + ".feature.tbl")
                shutil.copyfile(src_tbl, dst_tbl)
                _append_special_rows(dst_tbl, ann)
        os.makedirs(wd, exist_ok=True)   # same reason as above
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
        #  The flags column used to be truncated at 200 characters. Anything
        #  reading it back -- gen_facts.py's tally, any per-species breakdown --
        #  then splits a half-written flag into its own category: a Rotavirus
        #  panel produced "Genome is missing essent" and "Mature C-" as if they
        #  were distinct findings. The column is the audit trail; it is written
        #  in full.
        rows.append((b, nf, len(gf), len(cf), nff, "; ".join(allf)))

    #  Essentiality audit. copy_num marks a feature essential, and an
    #  essential feature that most genomes do not carry flags all of them --
    #  three modules have scored 0% this way (Quinvirinae, Closterovirus,
    #  Ampelovirus) while their core features called at 91-100%. This is the
    #  only place that can see the answer, because the test is the fraction of
    #  GENOMES carrying the feature and that is what a panel measures.
    try:
        blk = json.load(open(os.path.join(R, "Viral_PSSM.json")))
        fam = None
        for q in sorted(glob.glob(os.path.join(args.out, "*.qual.gto")))[:1]:
            fam = (json.load(open(q)).get("viral_family") or [None])[0] \
                  if isinstance(json.load(open(q)).get("viral_family"), list) \
                  else json.load(open(q)).get("viral_family")
        feats = (blk.get(fam) or {}).get("features") or {}
        if feats:
            called = collections.Counter(); nrt = 0
            for q in sorted(glob.glob(os.path.join(args.out, "*.qual.gto"))):
                d = json.load(open(q))
                if not d.get("features"): continue
                nrt += 1
                fns = {(f.get("function") or "") for f in d["features"]}
                for k, v in feats.items():
                    if v.get("anno") in fns: called[k] += 1
            thin = [(k, called[k] / nrt) for k, v in feats.items()
                    if v.get("copy_num") and nrt and called[k] / nrt < 0.85]
            if thin:
                print("\n  ESSENTIAL BUT RARE -- these carry copy_num and so flag every")
                print("  genome that lacks them. Drop copy_num unless the feature is")
                print("  genuinely present and merely uncallable:")
                for k, f in sorted(thin, key=lambda x: x[1]):
                    print("    %-8s called on %5.1f%% of routed genomes" % (k, 100.0 * f))

            #  The complement, and the reason this block exists at all. The
            #  safe build order is to ship with NO copy_num, measure here,
            #  and then add it to whatever earned it -- a missing copy_num
            #  costs nothing, a wrong one fails every genome in the taxon.
            #  Reporting only the first direction left the second half of
            #  that workflow with nothing to read.
            fat = [(k, called[k] / nrt) for k, v in feats.items()
                   if not v.get("copy_num") and nrt and called[k] / nrt >= 0.85]
            if fat:
                print("\n  EARNS copy_num -- called on >= 85%% of routed genomes and")
                print("  not currently marked essential. Add copy_num: 1 to these:")
                for k, f in sorted(fat, key=lambda x: -x[1]):
                    print("    %-8s called on %5.1f%% of routed genomes" % (k, 100.0 * f))

            print("\n  call rate, every declared feature (%d routed genome(s)):" % nrt)
            for k in sorted(feats, key=lambda x: -called[x]):
                print("    %-11s %5.1f%%  %s" % (k, 100.0 * called[k] / nrt if nrt else 0.0,
                                                 "copy_num" if feats[k].get("copy_num") else ""))
    except Exception:
        pass

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
