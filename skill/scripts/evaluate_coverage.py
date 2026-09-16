#!/usr/bin/env python3
"""Run the finished module over the taxon's whole BV-BRC holding and report.

evaluate_module.py scores a handful of hand-picked reference genomes. That is a
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
    """One representative sequence per genome file -> genomes.fna."""
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
            seqs[gid] = best
            fh.write(">%s\n%s\n" % (gid, best))
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
    cmd = ["perl", os.path.join(repo, "Other_Scripts", "New-annotate-viral-taxon.pl"),
           "-contigs", contigs_dir, "-meta", meta, "-t", str(threads),
           "-c", os.path.join(work, "Contig_GTOs"),
           "-a", os.path.join(work, "Anno_GTOs"),
           "-q", qdir, "-j", os.path.join(repo, "Viral_PSSM.json")]
    if taxon:
        cmd += ["-i", taxon]
    print("  " + " ".join(cmd))
    r = subprocess.run(cmd, env=dict(os.environ, LOWVAN_DATA_DIR=repo))
    if r.returncode != 0 or not os.path.isdir(qdir):
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
                fh.write(">%s\n%s\n" % (g, seqs[g]))
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
    if specials:
        print()
        print("  WARNING: the following features are called by an external program and")
        print("           CANNOT appear in this evaluation. Their coverage here is not")
        print("           low -- it is absent. Run New-annotate-viral-taxon.pl to see them.")
        for m, sp in sorted(specials.items()):
            print("             %-24s %s" % (m, ", ".join(sp)))
        print()
    default_floor = args.min_complete or (min(floors.values()) if floors else 10000)
    if floors:
        print("  completeness floor per module: %s"
              % ", ".join("%s=%d" % kv for kv in sorted(floors.items())))
    print("  floor applied to genomes of unknown module: %d" % default_floor)

    ex, partial_only = [], 0
    for _rep, mems in members.items():
        cand = [m for m in mems if len(seqs.get(m, "")) >= default_floor]
        if not cand:
            partial_only += 1
            continue
        ex.append(max(cand, key=lambda m: len(seqs[m])))
    ex.sort()
    print("  %d exemplar(s); %d cluster(s) below every floor, skipped"
          % (len(ex), partial_only))

    outd = annotate(ex, seqs, work, os.path.abspath(args.repo), args.threads)
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
        rows.append({"genome": g, "len": len(seqs[g]), "module": mod, "n_calls": len(calls),
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
