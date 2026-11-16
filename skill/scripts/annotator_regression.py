#!/usr/bin/env python3
"""Run two versions of annotate_by_viral_pssm.pl over a fixed genome set and
report every call that differs.

Why this exists. The annotator is shared by every LowVan module, so a change
made to fix one family is a change to all of them. On 12 Nov 2026 a ranking fix
for adenovirus (prefer a profile whose HSP can be delivered without a stop)
was first written so bluntly that it treated a CDS's own terminal stop as a
defect: E4ORF5 dropped from 137 to 79 aa, K100 from 695 to 616, and dozens of
features lost their last residue -- on the very family it was meant to fix.
The curator's rule since: any annotator change is re-run against this set,
old against new, and every difference is explained before it is committed.

The set (skill/assets/annotator_regression/genomes.tsv) is one complete genome
per installed non-adenovirus module -- all segments for segmented taxa -- plus
the 11 adenovirus genomes that the 12 Nov change altered or had to leave
alone. The contigs are stored with it so the set never drifts with a data
release.

    python3 annotator_regression.py \
        --old  $LOWVAN_DATA_DIR/annotate_by_viral_pssm.pl.pre_change \
        --new  $LOWVAN_KIT/annotate/annotate_by_viral_pssm.pl \
        --out  regression_run/  [--threads 4] [--jobs 4] [--only ID,ID]

Both versions run against the SAME LOWVAN_DATA_DIR, so the only variable is
the code. Outputs land in <out>/old.<id>.* and <out>/new.<id>.*; existing
non-empty outputs are reused, so an interrupted run resumes. The report
(<out>/REPORT.md, also printed) lists per genome: the module it was routed
to, the feature count, and every feature whose calls differ with the old and
new (start, end, strand, length aa, profile).

A difference is not a failure. It is a finding to be explained, one by one,
in the commit message or the review that ships the change.
"""
import argparse, collections, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ASSETS = os.path.join(HERE, "..", "assets", "annotator_regression")


def read_set(path):
    rows = []
    for l in open(path):
        if l.startswith("#") or not l.strip():
            continue
        c = l.rstrip("\n").split("\t")
        rows.append({"id": c[0], "module": c[1], "fna": os.path.join(ASSETS, c[2])})
    return rows


def run(script, fna, prefix, threads, maxlen):
    if os.path.exists(prefix + ".feature.tbl") and os.path.getsize(prefix + ".feature.tbl") > 0:
        return 0
    with open(prefix + ".err", "w") as err:
        return subprocess.run(["perl", script, "-i", fna, "-p", prefix, "-max", str(maxlen),
                               "-threads", str(threads)], stdout=subprocess.DEVNULL,
                              stderr=err).returncode


def calls(tbl):
    """feature.tbl -> {annotation: {(start, end, strand, aa, profile)}}; None if absent."""
    if not os.path.exists(tbl):
        return None
    d = collections.defaultdict(set)
    for l in open(tbl):
        c = l.rstrip("\n").split("\t")
        if len(c) < 17:
            continue
        prof = c[16].split(".", 1)[1] if "." in c[16] else c[16]
        d[c[6]].add((int(c[7]), int(c[8]), c[9], int(c[10]) // 3, prof))
    return d


def routed(err):
    """The module the annotator routed to, from its own log, if it says."""
    if not os.path.exists(err):
        return "?"
    for l in open(err, errors="replace"):
        if "Best" in l and "match" in l.lower():
            return l.strip()[:80]
    return "?"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--old", required=True, help="annotate_by_viral_pssm.pl before the change")
    ap.add_argument("--new", required=True, help="annotate_by_viral_pssm.pl after the change")
    ap.add_argument("--out", required=True)
    ap.add_argument("--set", default=os.path.join(ASSETS, "genomes.tsv"))
    ap.add_argument("--threads", type=int, default=2, help="blast threads per annotation")
    ap.add_argument("--jobs", type=int, default=4, help="annotations run at once")
    ap.add_argument("--max", type=int, default=60000, help="annotator -max (longest genome, nt)")
    ap.add_argument("--only", default="", help="comma-separated ids")
    a = ap.parse_args()
    if not os.environ.get("LOWVAN_DATA_DIR"):
        sys.exit("set LOWVAN_DATA_DIR: both versions must read the same data")
    os.makedirs(a.out, exist_ok=True)
    rows = read_set(a.set)
    if a.only:
        keep = set(a.only.split(","))
        rows = [r for r in rows if r["id"] in keep]
    jobs = []
    for r in rows:
        for tag, script in (("old", a.old), ("new", a.new)):
            jobs.append((script, r["fna"], os.path.join(a.out, "%s.%s" % (tag, r["id"]))))
    with ThreadPoolExecutor(a.jobs) as ex:
        rcs = list(ex.map(lambda j: run(*j, a.threads, a.max), jobs))
    failed = [j[2] for j, rc in zip(jobs, rcs) if rc]

    out = ["# Annotator regression", "", "old: `%s`" % a.old, "new: `%s`" % a.new,
           "data: `%s`" % os.environ["LOWVAN_DATA_DIR"], ""]
    tot = collections.Counter()
    for r in rows:
        o = calls(os.path.join(a.out, "old.%s.feature.tbl" % r["id"]))
        n = calls(os.path.join(a.out, "new.%s.feature.tbl" % r["id"]))
        if o is None or n is None:
            out.append("## %s (%s): NO OUTPUT from %s" % (r["id"], r["module"],
                       "old" if o is None else "new"))
            tot["no_output"] += 1
            continue
        diff = sorted(k for k in set(o) | set(n) if o.get(k, set()) != n.get(k, set()))
        tot["genomes"] += 1
        tot["features"] += len(set(o) | set(n))
        tot["differing"] += len(diff)
        out.append("## %s (%s): %d features, %d differ" % (r["id"], r["module"],
                   len(set(o) | set(n)), len(diff)))
        for k in diff:
            out.append("- **%s**" % k)
            out.append("  - old: %s" % sorted(o.get(k, set())))
            out.append("  - new: %s" % sorted(n.get(k, set())))
    out += ["", "**Total:** %d genomes, %d features, %d differ; %d genomes with no output."
            % (tot["genomes"], tot["features"], tot["differing"], tot["no_output"])]
    if failed:
        out += ["", "Non-zero exits (read the .err files): " + ", ".join(sorted(failed))]
    rep = "\n".join(out) + "\n"
    open(os.path.join(a.out, "REPORT.md"), "w").write(rep)
    print(rep)


if __name__ == "__main__":
    main()
