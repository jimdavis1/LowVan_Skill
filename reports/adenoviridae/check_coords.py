#!/usr/bin/env python3
"""Compare each collection's length distribution against the REFERENCE genome.

The curator's question is whether the features have proper coordinates. A
collection cannot be checked against itself: a profile built from a mixture
of two genes has a perfectly self-consistent length distribution. What it can
be checked against is the annotated length of that gene in a genome somebody
curated -- so this reads the reference lengths from a table written by hand
from the GenBank record and flags every collection whose median is far from
it.

A big gap means one of three things, and all three have happened in this
family: the collection is mixed (two genes in one), the collection is really
a DIFFERENT gene than its name says (Barthadenovirus LH3 vs "E1B 55K"), or
the gene is genuinely a different size in other species of the genus, which
is a length WINDOW to widen rather than a defect.
"""
import json, os, statistics, sys

TOL = 0.25   # fraction of the reference length

def lens(p):
    if not os.path.exists(p):
        return []
    out, cur = [], 0
    for l in open(p):
        if l.startswith(">"):
            if cur:
                out.append(cur)
            cur = 0
        else:
            cur += len(l.strip())
    if cur:
        out.append(cur)
    return out

def main():
    ref = json.load(open(sys.argv[1]))
    for mod, feats in ref.items():
        print("=" * 70)
        print("%s   reference %s" % (mod, feats.pop("_reference", "?")))
        print("%-12s %7s %7s %7s %7s  %s" %
              ("feature", "ref_aa", "n", "median", "range", "verdict"))
        for k in sorted(feats):
            r = feats[k]
            L = sorted(lens(os.path.join(mod, "collections", "%s.fasta" % k)))
            if not L:
                print("%-12s %7d %7d %7s %7s  EMPTY" % (k, r, 0, "-", "-"))
                continue
            m = statistics.median(L)
            d = abs(m - r) / float(r)
            v = "ok" if d <= TOL else "OFF by %.0f%%" % (100 * d)
            print("%-12s %7d %7d %7.0f %4d-%-4d  %s"
                  % (k, r, len(L), m, L[0], L[-1], v))

main()
