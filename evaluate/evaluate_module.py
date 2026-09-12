#!/usr/bin/env python3
"""Run an installed LowVan module over held-out genomes and score the calls.

For each accession it fetches the nucleotide sequence and the GenBank CDS
coordinates, runs annotate_by_viral_pssm.pl, and compares. Reference CDS come
from the submitter, so they are a yardstick, not truth -- the point is to see
which genes the module misses and where it doubles up, not to match GenBank.

Two failures matter more than the recall number:

  DUPLICATE   two features called on the same coordinates. Always a bug: some
              cluster is binned under the wrong feature and out-scores the
              right PSSM at that locus. Run qc_cross_feature.py.
  MISSED      a reference CDS with no call. Either a feature with no PSSM, or
              a bit_cutoff set too high.

    python3 evaluate_module.py --repo ~/Viral_Annotation --acc NC_001542,EF206707
    python3 evaluate_module.py --repo ~/Viral_Annotation --acc-file accs.txt --out eval.tsv
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def fetch(acc, rettype):
    url = "%s/efetch.fcgi?%s" % (EUTILS, urllib.parse.urlencode(
        {"db": "nuccore", "id": acc, "rettype": rettype, "retmode": "text"}))
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=90) as fh:
                return fh.read().decode("utf8", "replace")
        except Exception:
            if attempt == 2:
                raise
            time.sleep(3)


def genbank_cds(gb):
    """[(start, end, product)] from a GenBank flatfile, 1-based inclusive."""
    out = []
    cur = None
    for line in gb.split("\n"):
        m = re.match(r"^     (CDS|mat_peptide)\s+(\S.*)$", line)
        if m:
            if cur:
                out.append(cur)
            loc = m.group(2)
            span = re.findall(r"(\d+)\.\.>?(\d+)", loc)
            if span:
                cur = [int(span[0][0]), int(span[-1][1]), "", m.group(1)]
            else:
                cur = None
            continue
        if cur is not None:
            p = re.match(r'^\s+/product="([^"]*)"?', line)
            if p:
                cur[2] = p.group(1)
    if cur:
        out.append(cur)
    return [c for c in out if c[3] == "CDS"]


def overlap(a1, a2, b1, b2):
    lo, hi = max(a1, b1), min(a2, b2)
    return max(0, hi - lo + 1)


def run_one(repo, acc, tmp, threads):
    fa = os.path.join(tmp, acc + ".fa")
    with open(fa, "w") as fh:
        fh.write(fetch(acc, "fasta"))
    gb = fetch(acc, "gbwithparts")
    ref = genbank_cds(gb)

    pre = os.path.join(tmp, acc)
    env = dict(os.environ, LOWVAN_DATA_DIR=repo)
    r = subprocess.run(["perl", os.path.join(repo, "annotate_by_viral_pssm.pl"),
                        "-i", fa, "-p", pre, "-threads", str(threads)],
                       capture_output=True, text=True, env=env)

    tbl = pre + ".feature.tbl"
    #  A crashed annotator and a module that calls nothing look identical in
    #  the score table: both report 0 calls and 0% recall. That is the worst
    #  possible failure mode for an evaluation tool, because 0% reads as a
    #  finding about the module. It is usually the environment -- PERL5LIB
    #  without File::Slurp or gjoseqlib is enough to abort before a single
    #  PSSM runs. Never let a non-run be scored as a result.
    if r.returncode != 0 or not os.path.exists(tbl):
        err = (r.stderr or r.stdout or "").strip().splitlines()
        raise SystemExit(
            "annotate_by_viral_pssm.pl failed on %s (exit %d) and produced no\n"
            "feature table, so there is nothing to score. Last output:\n\n  %s\n\n"
            "Most often PERL5LIB: the annotator needs File::Slurp, JSON::XS,\n"
            "gjoseqlib.pm and BlastInterface.pm. Run one genome by hand first."
            % (acc, r.returncode, "\n  ".join(err[-6:]) or "(no output)"))
    calls = []
    if os.path.exists(tbl):
        for line in open(tbl):
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            calls.append({"type": f[4], "key": f[6], "start": int(f[7]),
                          "end": int(f[8]), "module": f[11]})
    return ref, calls


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    ap.add_argument("--acc", help="comma-separated accessions")
    ap.add_argument("--acc-file", help="file with one accession per line")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", help="write a per-feature TSV here")
    args = ap.parse_args()

    accs = []
    if args.acc:
        accs += [a.strip() for a in args.acc.split(",") if a.strip()]
    if args.acc_file:
        accs += [l.split("#")[0].strip() for l in open(args.acc_file)
                 if l.split("#")[0].strip()]
    if not accs:
        sys.exit("give --acc or --acc-file")

    repo = os.path.abspath(os.path.expanduser(args.repo))
    tmp = tempfile.mkdtemp(prefix="lveval.")
    rows, totals = [], {"ref": 0, "called": 0, "hit": 0, "dup": 0, "extra": 0}

    print("%-12s %-20s %5s %5s %5s %5s %5s   %s"
          % ("accession", "module", "ref", "call", "hit", "miss", "dup", "missed products"))
    for acc in accs:
        try:
            ref, calls = run_one(repo, acc, tmp, args.threads)
        except Exception as exc:
            print("%-12s ERROR %s" % (acc, exc))
            continue
        cds = [c for c in calls if c["type"] == "CDS"]
        module = cds[0]["module"] if cds else (calls[0]["module"] if calls else "-")

        # duplicates: same coordinates called under two keys
        seen, dups = {}, []
        for c in cds:
            k = (c["start"], c["end"])
            if k in seen:
                dups.append((seen[k], c["key"], k))
            else:
                seen[k] = c["key"]

        matched, missed = set(), []
        for rs, re_, prod, _t in ref:
            best = None
            for i, c in enumerate(cds):
                ov = overlap(rs, re_, c["start"], c["end"])
                if ov > 0.5 * (re_ - rs + 1):
                    best = i
                    break
            if best is None:
                missed.append(prod or "?")
            else:
                matched.add(best)
        extra = len(cds) - len(matched)

        totals["ref"] += len(ref); totals["called"] += len(cds)
        totals["hit"] += len(matched); totals["dup"] += len(dups)
        totals["extra"] += max(0, extra)
        print("%-12s %-20s %5d %5d %5d %5d %5d   %s"
              % (acc, module[:20], len(ref), len(cds), len(matched),
                 len(missed), len(dups), "; ".join(m[:26] for m in missed[:3])))
        for a, b, k in dups:
            print("%-12s   DUPLICATE %s and %s both at %d-%d" % ("", a, b, k[0], k[1]))
        for c in cds:
            rows.append((acc, module, c["key"], c["start"], c["end"]))

    print("\n  reference CDS %d | called %d | matched %d (%.0f%%) | duplicates %d | unmatched calls %d"
          % (totals["ref"], totals["called"], totals["hit"],
             100.0 * totals["hit"] / max(totals["ref"], 1), totals["dup"], totals["extra"]))
    if totals["dup"]:
        print("  duplicates are always a defect -- run qc_cross_feature.py")

    if args.out:
        with open(args.out, "w") as fh:
            fh.write("accession\tmodule\tfeature\tstart\tend\n")
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")
        print("  wrote %s" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
