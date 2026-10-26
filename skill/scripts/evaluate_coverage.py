#!/usr/bin/env python3
"""Run the finished module over the taxon's whole BV-BRC holding and report.

run_gto_eval.py scores a held-out panel of a few dozen genomes. That is a
correctness check, not a coverage one: the panel is small, curated, and usually
the same genomes the profiles were built from. This runs the module over
everything BV-BRC has for the taxon and answers four questions the small panel
cannot:

  1. How many genomes does the BLASTn routing step miss entirely? A genome no
     rep contig claims is never annotated at all, whatever the profiles can do.
     Misses that cluster in one genus mean a missing rep contig.
  2. How many genomes come out good versus poor quality?
  3. When poor, which flags fire? This is the audit of min_len/max_len and
     copy_num -- if "Feature is too short" dominates, the length bounds are
     wrong, not the genomes.
  4. What is the distribution of called proteins? In complete genomes the core
     proteins should appear in near-equal numbers. A core protein lagging the
     others is a coverage hole in that feature's profiles.

    # 1. download (uses the taxon's min/max from the JSON if it has segments)
    perl Other_Scripts/New-annotate-viral-taxon.pl -i Rhabdoviridae \
         -download-only -d Contigs -metaout Rhabdoviridae.metadata \
         -min 1000 -max 50000

    # 2. cluster, pick exemplars, annotate, report
    python3 evaluate_coverage.py --contigs Contigs --repo ~/Viral_Annotation \
            --min-complete 10000 --threads 6 --out coverage.json

Why cluster first: BV-BRC holds thousands of near-identical genomes for the
well-sequenced species, and annotating all of them measures sequencing effort
rather than module coverage. One exemplar per 95% nucleotide cluster gives each
distinct genome type one vote.

Why --min-complete: most BV-BRC records for a virus family are single-gene
partial submissions. Those are flagged poor for missing essential features,
which is correct and tells you nothing about the module. Restrict to
near-complete genomes and question 4 becomes diagnostic.
"""

import argparse
import collections
import glob
import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor


def read_contigs(path):
    seqs, cur = [], []
    for line in open(path, errors="replace"):
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
            cur = []
        else:
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


def stage(contigs, work):
    """Genomes -> genomes.fna, plus every contig of each genome.

    The longest contig is what goes into genomes.fna, because that is only a
    fingerprint for the redundancy clustering. ANNOTATION must see the whole
    genome: a segmented taxon keeps its genes on separate molecules, and
    staging one contig means only the longest segment is ever annotated.

    That is not hypothetical. Rotavirus has 11 segments and its polymerase sits
    on the longest; evaluated this way the module called RDRP on 82% of routed
    genomes and every other feature on under 11%, because the other ten
    segments were never shown to the annotator. 356 of 391 routed genomes
    scored exactly one feature.
    """
    os.makedirs(work, exist_ok=True)
    fa = os.path.join(work, "genomes.fna")
    n = 0
    seqs = {}
    with open(fa, "w") as fh:
        for f in sorted(glob.glob(os.path.join(contigs, "*"))):
            if os.path.isdir(f):
                continue
            s = read_contigs(f)
            if not s:
                continue
            gid = re.sub(r"\.(contigs|fna|fa|fasta)$", "", os.path.basename(f))
            best = max(s, key=len)
            seqs[gid] = list(s)          # every contig, for annotation
            fh.write(">%s\n%s\n" % (gid, best))   # longest only, for clustering
            n += 1
    print("  staged %d genome(s) -> %s" % (n, fa))
    return fa, seqs


def cluster(fa, work, pid, cov, threads):
    pre = os.path.join(work, "clu")
    tsv = pre + "_cluster.tsv"
    if not os.path.exists(tsv):
        subprocess.run(["mmseqs", "easy-cluster", fa, pre, os.path.join(work, "tmp_clu"),
                        "--min-seq-id", str(pid), "-c", str(cov), "--cov-mode", "1",
                        "--threads", str(threads)], capture_output=True)
    members = collections.defaultdict(list)
    for line in open(tsv):
        rep, mem = line.split()[:2]
        members[rep].append(mem)
    print("  %d cluster(s) at %.0f%% nucleotide identity" % (len(members), 100 * pid))
    return members


#  QUALITY MUST COME FROM THE GTO, NOT THE FEATURE TABLE.
#
#  annotate_by_viral_pssm.pl writes a flat feature table, and it is tempting to
#  assess quality straight from it: compare each row's length to min_len and
#  count copies. That reimplements viral_genome_quality.pl badly and it will
#  mislead you in two ways.
#
#  First, the table is one row per called region, not one row per feature. A
#  polymerase broken by a frameshift comes back as several rows with separate
#  feature ids -- L was fragmented in 18% of Rhabdoviridae genomes -- and
#  comparing each fragment's length to the full-length floor reports a flood of
#  "Feature is too short" that is an artefact of the comparison.
#
#  Second, and worse, the flat table has no representation of special features.
#  Rhabdoviridae happens to contain no splice variants or transcript-edited
#  products in the modules that ship, so the shortcut appeared to work. A taxon
#  with a spliced or edited feature would be scored as if those features did
#  not exist.
#
#  So run the real pipeline: New-annotate-viral-taxon.pl in contig-dir mode
#  builds Contig_GTOs, Anno_GTOs and Quality_GTOs, and the Quality GTOs carry
#  genome_quality, feature_quality and feature_quality_flags as the annotator
#  itself computed them.
#
#  This needs the BV-BRC dev kit (GenomeTypeObject.pm, IDclient.pm,
#  rast-create-genome). Those are not on PyPI or CPAN and pulling individual
#  .pm files from GitHub does not work -- the dependency chain keeps going.
#  Build the dev container, or run this step on a machine that has it.


def annotate_gto(contigs_dir, meta, work, repo, taxon, threads):
    """Annotate via the GTO pipeline and return the Quality_GTOs directory."""
    qdir = os.path.join(work, "Quality_GTOs")
    #  The output directories MUST be relative and the process MUST run in
    #  `work`. New-annotate-viral-taxon.pl builds every output path as
    #  "$base/$dir/..." where $base is getcwd at startup, so an absolute -c
    #  lands at /work//work/Contig_GTOs. It does not fail on that: the inner
    #  redirect fails per genome, the pipeline runs to the end, and the script
    #  exits 0 with three empty directories. The only symptom upstream was
    #  "0 routed, 2 UNROUTED (100.0%)".
    os.makedirs(work, exist_ok=True)
    #  -contigs is prefixed with $base too (line 292: "--contigs
    #  $base/$contigDir/$item"), so it must be relative to `work` as well.
    #  -meta and -j are used as given and stay absolute.
    cmd = ["perl", os.path.join(repo, "Other_Scripts", "New-annotate-viral-taxon.pl"),
           "-contigs", os.path.relpath(os.path.abspath(contigs_dir), work),
           "-meta", os.path.abspath(meta), "-t", str(threads),
           "-c", "Contig_GTOs", "-a", "Anno_GTOs", "-q", "Quality_GTOs",
           "-j", os.path.join(repo, "Viral_PSSM.json")]
    if taxon:
        cmd += ["-i", taxon]
    print("  (cwd=%s) %s" % (work, " ".join(cmd)))
    r = subprocess.run(cmd, cwd=work, env=dict(os.environ, LOWVAN_DATA_DIR=repo))
    #  An empty Quality_GTOs is the failure above wearing a success exit code,
    #  so emptiness is checked, not just existence.
    if r.returncode != 0 or not os.path.isdir(qdir) or not os.listdir(qdir):
        raise SystemExit(
            "the GTO pipeline did not produce %s.\n"
            "It needs the BV-BRC dev kit on PERL5LIB/PATH: GenomeTypeObject.pm,\n"
            "IDclient.pm and rast-create-genome. Run this step where that is\n"
            "installed -- do not fall back to parsing the feature table, which\n"
            "cannot see special features and miscounts fragmented ones." % qdir)
    return qdir


def read_quality(qdir):
    """genome-level and flag-level tallies straight out of the Quality GTOs."""
    good = poor = 0
    flags = collections.Counter()
    per_feature = collections.defaultdict(collections.Counter)
    for f in sorted(glob.glob(os.path.join(qdir, "*"))):
        try:
            g = json.load(open(f))
        except Exception:
            continue
        q = (g.get("genome_quality") or g.get("quality", {}).get("genome_quality"))
        if str(q).lower().startswith("good"):
            good += 1
        else:
            poor += 1
        for feat in g.get("features", []):
            fam = feat.get("function") or feat.get("product") or "?"
            per_feature[fam][feat.get("feature_quality", "?")] += 1
            for fl in feat.get("feature_quality_flags", []) or []:
                flags[fl] += 1
    return good, poor, flags, per_feature


def annotate(exemplars, seqs, work, repo, threads):
    exd = os.path.join(work, "ex"); os.makedirs(exd, exist_ok=True)
    outd = os.path.join(work, "ann"); os.makedirs(outd, exist_ok=True)
    env = dict(os.environ, LOWVAN_DATA_DIR=repo)
    jobs = []
    for g in exemplars:
        fa = os.path.join(exd, g + ".fna")
        if not os.path.exists(fa):
            with open(fa, "w") as fh:
                for i, c in enumerate(seqs[g], 1):
                    fh.write(">%s_%d\n%s\n" % (g, i, c))
        jobs.append((g, fa))

    def one(job):
        g, fa = job
        pre = os.path.join(outd, g)
        if os.path.exists(pre + ".feature.tbl") or os.path.exists(pre + ".NOCALL"):
            return
        r = subprocess.run(["perl", os.path.join(repo, "annotate_by_viral_pssm.pl"),
                            "-i", fa, "-p", pre, "-threads", "1"],
                           capture_output=True, text=True, env=env)
        if not os.path.exists(pre + ".feature.tbl"):
            with open(pre + ".NOCALL", "w") as fh:
                fh.write((r.stderr or r.stdout or "")[-3000:])

    with ThreadPoolExecutor(max_workers=threads) as pool:
        for i, _ in enumerate(pool.map(one, jobs), 1):
            if i % 100 == 0:
                print("    %d/%d" % (i, len(jobs)), flush=True)
    return outd



def gto_ok(work, g):
    """True when this exemplar actually produced a non-empty annotation."""
    for d in ("Anno_GTOs", "Quality_GTOs"):
        for pat in (g + ".anno.gto", g + ".qual.gto"):
            p = os.path.join(work, d, pat)
            if os.path.exists(p) and os.path.getsize(p) > 0:
                return True
    return False


def stage_exemplars(ex, seqs, work, sub="ex_gto", meta_name="exemplars.metadata"):
    """One FASTA per exemplar plus the metadata New-annotate needs.

    Contig-dir mode wants `file_name TAB genome_name TAB taxon_id`. The taxon
    id is the part of a BV-BRC genome id before the dot; anything else gets
    10239, the Viruses node, which is only used to stamp the GTO.
    """
    exd = os.path.join(work, sub)
    os.makedirs(exd, exist_ok=True)
    meta = os.path.join(work, meta_name)
    with open(meta, "w") as mh:
        for g in ex:
            fa = os.path.join(exd, g + ".fna")
            if not os.path.exists(fa) or os.path.getsize(fa) == 0:
                with open(fa, "w") as fh:
                    for i, c in enumerate(seqs[g], 1):
                        fh.write(">%s_%d\n%s\n" % (g, i, c))
            head = g.split(".")[0]
            tax = head if head.isdigit() else "10239"
            mh.write("%s\t%s\t%s\n" % (os.path.basename(fa), g, tax))
    print("  staged %d exemplar(s) -> %s" % (len(ex), exd))
    return exd, meta


def report_gto(ex, seqs, work, out):
    """Same report as report(), read from the GTOs so specials are included.

    Module and feature come from family_assignments, which both the PSSM pass
    and the special passes fill in: 'HIV2.NEF.9' from a profile, 'HIV2.REV'
    from the splice caller. The fourth element names the program, so a feature
    can be attributed to the pass that called it.

    The flat path read its 'flags' out of column 13, which is the reference
    contig, not a flag -- every previous coverage.json carries a flags block
    that is really a tally of which rep contig routed how many genomes. Here
    the flags are the real feature_quality_flags from the Quality GTOs.
    """
    annod = os.path.join(work, "Anno_GTOs")
    quald = os.path.join(work, "Quality_GTOs")
    unrouted, per_mod = [], collections.Counter()
    feats, flags = collections.Counter(), collections.Counter()
    per_mod_feat = collections.defaultdict(collections.Counter)
    by_caller = collections.Counter()
    rows = []

    #  Exact names only. A glob of "<id>*" is wrong twice here: Anno_GTOs also
    #  holds <id>.stdout.txt and <id>.stderr.txt, which parse as nothing, and
    #  genome ids are prefixes of each other -- 11665.150* also matches
    #  11665.1505. Both failures look identical downstream: "UNROUTED".
    def gto_for(d, g):
        for pat in (g + ".anno.gto", g + ".qual.gto", g + ".gto", g + ".ann.gto"):
            p = os.path.join(d, pat)
            if os.path.exists(p) and os.path.getsize(p) > 0:
                return p
        return None

    for g in ex:
        p = gto_for(annod, g) or gto_for(quald, g)
        if not p:
            unrouted.append(g); continue
        try:
            d = json.load(open(p))
        except Exception:
            unrouted.append(g); continue
        mod = None; calls = []
        for f in d.get("features", []) or []:
            fam = (f.get("family_assignments") or [[None, "", "", ""]])[0]
            key = (fam[1] or "") if len(fam) > 1 else ""
            if not key:
                continue
            bits = key.split(".")
            mod = bits[0] or mod
            gene = dict(f.get("alias_pairs") or []).get("gene")
            calls.append(gene or (bits[1].lower() if len(bits) > 1 else key))
            if len(fam) > 3 and fam[3]:
                by_caller[fam[3]] += 1
        if not calls:
            unrouted.append(g); continue
        per_mod[mod or "(none)"] += 1
        for c in calls:
            feats[c] += 1
            per_mod_feat[mod or "(none)"][c] += 1
        rows.append({"genome": g, "len": sum(map(len, seqs[g])), "module": mod,
                     "n_calls": len(calls), "calls": sorted(set(calls))})

    good = poor = 0
    for q in sorted(glob.glob(os.path.join(quald, "*"))):
        if os.path.isdir(q):
            continue
        try:
            d = json.load(open(q))
        except Exception:
            continue
        qa = d.get("genome_quality") or (d.get("quality") or {}).get("genome_quality")
        if str(qa).lower().startswith("good"):
            good += 1
        else:
            poor += 1
        for f in d.get("features", []) or []:
            for fl in f.get("feature_quality_flags") or []:
                flags[fl] += 1

    res = {"exemplars": len(ex), "unrouted": unrouted, "per_module": dict(per_mod),
           "features": dict(feats),
           "per_module_features": {k: dict(v) for k, v in per_mod_feat.items()},
           "flags": dict(flags), "rows": rows,
           "genome_quality": {"good": good, "poor": poor},
           "calls_by_program": dict(by_caller),
           "source": "GTO pipeline -- special features included"}
    json.dump(res, open(out, "w"), indent=1)
    print("\n  %d exemplars | %d routed | %d UNROUTED (%.1f%%)"
          % (len(ex), len(ex) - len(unrouted), len(unrouted),
             100.0 * len(unrouted) / max(len(ex), 1)))
    if by_caller:
        print("  calls by program:")
        for k, v in sorted(by_caller.items(), key=lambda x: -x[1]):
            print("    %-40s %d" % (k, v))
    print("  genome quality: %d good, %d poor" % (good, poor))
    print("  wrote %s" % out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs", required=True, help="directory of per-genome contig files")
    ap.add_argument("--repo", required=True)
    ap.add_argument("--work", default="coverage_eval")
    ap.add_argument("--pident", type=float, default=0.95)
    ap.add_argument("--cov", type=float, default=0.80)
    ap.add_argument("--min-complete", type=int, default=0,
                    help="completeness floor in bases. Default 0 means take it "
                         "per module from the JSON's declared segment min_len, "
                         "which is the only correct choice for a taxon with "
                         "segmented members -- see --json.")
    ap.add_argument("--json", default=None,
                    help="module JSON, for per-module completeness thresholds")
    ap.add_argument("--threads", type=int, default=6)
    ap.add_argument("--out", default="coverage.json")
    ap.add_argument("--taxon", default=None,
                    help="taxon name, passed to New-annotate-viral-taxon.pl")
    ap.add_argument("--gto-retries", type=int, default=3,
                    help="rounds of sequential retry for exemplars whose "
                         "Contig GTO came back empty (ID-server 403).")
    ap.add_argument("--no-special", action="store_true",
                    help="force the flat feature-table path even when the JSON "
                         "declares spliced or transcript-edited features. The "
                         "report then CANNOT contain them -- only use this when "
                         "the BV-BRC dev kit is unavailable.")
    args = ap.parse_args()

    work = os.path.abspath(args.work)
    fa, seqs = stage(args.contigs, work)
    members = cluster(fa, work, args.pident, args.cov, args.threads)

    #  A single completeness floor is wrong the moment the taxon contains a
    #  segmented member. Rhabdoviridae looked like it had no complete
    #  Dichorhavirus genome anywhere in BV-BRC -- 23 clusters, all "partial".
    #  They are complete: Dichorhavirus is bisegmented and its segments are
    #  5310-7371 bases, so a 10 kb floor excluded the entire genus by
    #  construction and the evaluation would have silently said nothing about
    #  a module that ships 26 profiles.
    #
    #  So take the floor from each module's own declared segment min_len.
    floors = {}
    specials = {}
    if args.json and os.path.exists(args.json):
        j = json.load(open(args.json))
        for m, b in j.items():
            segs = b.get("segments") or {}
            if segs:
                floors[m] = min(v["min_len"] for v in segs.values())
            sp = sorted(k for k, e in (b.get("features") or {}).items()
                        if e.get("special"))
            if sp:
                specials[m] = sp
    #  This function annotates with annotate(), the flat-feature-table path,
    #  which cannot represent a special feature at all -- see the long comment
    #  above annotate_gto(). If the JSON declares any, say so loudly rather
    #  than silently reporting 0% for them. Togaviridae's TF is called on 293
    #  of 381 genomes by the full pipeline and appears nowhere in this report.
    #
    #  THIS USED TO BE A WARNING AND NOTHING ELSE.
    #
    #  It printed the list and then ran the flat path anyway, so every coverage
    #  audit of a module with a special feature reported that feature at zero
    #  without ever saying zero -- it simply was not a row. The lentivirus
    #  audit of 25 Oct measured HIV2 and EIAV with REV and TAT missing, and the
    #  warning scrolled past in a log. A warning that does not change what the
    #  program does is a comment addressed to no one.
    #
    #  So the declaration now decides the route: if any module in scope
    #  declares a special feature, annotation goes through the GTO pipeline,
    #  which pipes annotate_by_viral_pssm-GTO.pl into
    #  get_transcript_edited_features.pl and get_splice_variant_features.pl.
    #  --no-special restores the old behaviour and says loudly what it costs.
    use_gto = bool(specials) and not args.no_special
    if specials:
        print()
        for m, sp in sorted(specials.items()):
            print("    %-24s special: %s" % (m, ", ".join(sp)))
        if use_gto:
            print("  -> annotating via the GTO pipeline so these are measured.")
        else:
            print("  WARNING: --no-special given. These features are called by an")
            print("           external program and CANNOT appear in this report.")
            print("           Their coverage here is not low -- it is absent.")
        print()
    default_floor = args.min_complete or (min(floors.values()) if floors else 10000)
    if floors:
        print("  completeness floor per module: %s"
              % ", ".join("%s=%d" % kv for kv in sorted(floors.items())))
    print("  floor applied to genomes of unknown module: %d" % default_floor)

    ex, partial_only = [], 0
    for _rep, mems in members.items():
        cand = [m for m in mems if sum(map(len, seqs.get(m, []))) >= default_floor]
        if not cand:
            partial_only += 1
            continue
        ex.append(max(cand, key=lambda m: sum(map(len, seqs[m]))))
    ex.sort()
    print("  %d exemplar(s); %d cluster(s) below every floor, skipped"
          % (len(ex), partial_only))

    repo = os.path.abspath(args.repo)
    if use_gto:
        exd, meta = stage_exemplars(ex, seqs, work)
        annotate_gto(exd, meta, work, repo, args.taxon, args.threads)
        #  RETRY THE ID-SERVER 403s.
        #
        #  rast-create-genome registers a genome id, and the service rate
        #  limits on concurrency by answering **403 Forbidden**, not 429 --
        #  so the failure reads like an expired token and is not one.
        #  New-annotate-viral-taxon.pl has no retry: it leaves a zero-byte
        #  Contig GTO, the annotation pipeline runs on nothing, and the
        #  genome surfaces downstream as "UNROUTED". The genus half of the
        #  26 Oct lentivirus audit lost 23 of 101 exemplars that way, with a
        #  token good for another 3,315 hours; re-running one of them by hand
        #  failed once and then succeeded twice.
        #
        #  So retry the failures sequentially. --threads 1 is the point: the
        #  limit is on concurrency, and make_gto.py found that even 4 was
        #  enough to lose 89 of 120 genomes on Trivirinae.
        for rnd in range(1, args.gto_retries + 1):
            bad = [g for g in ex if not gto_ok(work, g)]
            if not bad:
                break
            print("  retry %d/%d: %d exemplar(s) have no annotation"
                  % (rnd, args.gto_retries, len(bad)), flush=True)
            time.sleep(min(30 * rnd, 120))
            rexd, rmeta = stage_exemplars(bad, seqs, work,
                                          sub="ex_retry%d" % rnd,
                                          meta_name="retry%d.metadata" % rnd)
            try:
                annotate_gto(rexd, rmeta, work, repo, args.taxon, 1)
            except SystemExit as e:
                print("  retry %d did not complete: %s" % (rnd, e))
        still = [g for g in ex if not gto_ok(work, g)]
        if still:
            print("  %d exemplar(s) still unannotated after %d retries: %s"
                  % (len(still), args.gto_retries, ", ".join(still[:10])))
        report_gto(ex, seqs, work, args.out)
    else:
        outd = annotate(ex, seqs, work, repo, args.threads)
        report(ex, seqs, outd, args.out)


def report(ex, seqs, outd, out):
    unrouted, per_mod, feats, flags = [], collections.Counter(), collections.Counter(), collections.Counter()
    per_mod_feat = collections.defaultdict(collections.Counter)
    quality = collections.Counter()
    rows = []
    for g in ex:
        tbl = os.path.join(outd, g + ".feature.tbl")
        if not os.path.exists(tbl):
            unrouted.append(g); continue
        mod = None; calls = []
        for line in open(tbl, errors="replace"):
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            calls.append(f[6])
            mod = f[11] or mod
            if len(f) > 12 and f[12].strip():
                for fl in f[12].split(";"):
                    if fl.strip():
                        flags[fl.strip()] += 1
        if not calls:
            unrouted.append(g); continue
        per_mod[mod or "(none)"] += 1
        for c in calls:
            feats[c] += 1
            per_mod_feat[mod or "(none)"][c] += 1
        rows.append({"genome": g, "len": sum(map(len, seqs[g])), "module": mod, "n_calls": len(calls),
                     "calls": sorted(set(calls))})
    res = {"exemplars": len(ex), "unrouted": unrouted, "per_module": dict(per_mod),
           "features": dict(feats),
           "per_module_features": {k: dict(v) for k, v in per_mod_feat.items()},
           "flags": dict(flags), "rows": rows}
    json.dump(res, open(out, "w"), indent=1)
    print("\n  %d exemplars | %d routed | %d UNROUTED (%.1f%%)"
          % (len(ex), len(ex) - len(unrouted), len(unrouted),
             100.0 * len(unrouted) / max(len(ex), 1)))
    print("  wrote %s" % out)


if __name__ == "__main__":
    sys.exit(main())
