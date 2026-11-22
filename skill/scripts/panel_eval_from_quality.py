#!/usr/bin/env python3
"""Build gen_facts.py's panel_eval.tsv from Quality_GTOs already on disk.

`run_gto_eval.py` produces this table, but it re-annotates to do it, and the
coverage run has already annotated and scored exactly this panel.  Re-running
the annotator to recover six columns per genome costs hours and can only
reproduce what `Quality_GTOs/` already holds.

Columns, in the order gen_facts.py reads them:
  0 genome
  1 n_features        features in the annotated GTO
  2 contig flags      non-empty Exceptions in <g>.contig_quality
  3 genome flags      1 when genome_quality is not Good
  4 feature flags     non-empty Exceptions in <g>.feature_quality
  5 flag text         "; "-joined, for the tally in section 06

Usage: panel_eval_from_quality.py <coverage work dir> <out.tsv>
"""
import collections
import json
import os
import sys

work, out = sys.argv[1], sys.argv[2]
qd = os.path.join(work, "Quality_GTOs")
ad = os.path.join(work, "Anno_GTOs")


def exceptions(path, col):
    """Non-empty Exceptions cells in a *_quality table, as a list."""
    if not os.path.exists(path):
        return []
    hits = []
    with open(path, errors="replace") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        try:
            i = head.index("Exceptions")
        except ValueError:
            return []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > i and f[i].strip():
                hits.append(f[i].strip())
    return hits


rows = []
for fn in sorted(os.listdir(qd)):
    if not fn.endswith(".qual.gto"):
        continue
    g = fn[:-9]
    try:
        gq = json.load(open(os.path.join(qd, fn))).get("quality", {}).get("genome_quality")
    except Exception:
        gq = None
    cflag = exceptions(os.path.join(qd, g + ".contig_quality"), "Exceptions")
    fflag = exceptions(os.path.join(qd, g + ".feature_quality"), "Exceptions")
    nf = 0
    ap = os.path.join(ad, g + ".anno.gto")
    if os.path.exists(ap):
        try:
            nf = len(json.load(open(ap)).get("features") or [])
        except Exception:
            nf = 0
    rows.append([g, str(nf), str(len(cflag)),
                 "0" if str(gq).lower() == "good" else "1",
                 str(len(fflag)), "; ".join(sorted(set(cflag + fflag)))])

with open(out, "w") as fh:
    fh.write("genome\tn_features\tcontig_flags\tgenome_flags\tfeature_flags\tflags\n")
    for r in rows:
        fh.write("\t".join(r) + "\n")

t = collections.Counter()
for r in rows:
    for f in r[5].split("; "):
        if f.strip():
            t[f.strip()] += 1
print("  %d genome(s) -> %s" % (len(rows), out))
print("  %d with nothing called, %d with no flag at all"
      % (sum(1 for r in rows if r[1] == "0"),
         sum(1 for r in rows if r[1] != "0" and r[2] == "0" and r[3] == "0" and r[4] == "0")))
for k, v in t.most_common(10):
    print("    %-50s %d" % (k[:50], v))
