#!/usr/bin/env python3
"""Measure how fast the vocabulary of annotation strings grows.

A controlled vocabulary and an uncontrolled one both label every protein. The
difference only becomes visible when you count DISTINCT strings as genomes
accumulate: a controlled vocabulary saturates, because the Nth genome reuses
names the first N-1 already established, while free text keeps growing roughly
linearly, because every submitter spells the same protein a new way.

Shuffle the genomes into random order, walk the order accumulating distinct
annotation strings, and average over replicates -- the rarefaction curve from
ecology, where the same question is how fast new species appear.

    python3 annotation_rarefaction.py --new coverage_eval/ann \
            --old bvbrc_annotations.tsv --out rarefaction.json --replicates 100

Reads the new annotations from the module's own feature tables and the old ones
from a TSV of genome id and product. Writes the two mean curves plus the final
ratio, ready to plot.

Report the ratio, not just the picture: strings per 100 genomes at the end of
each curve says how much collapse the vocabulary bought.
"""
import argparse, collections, glob, json, os, random, sys


def curve(sets, reps, rng):
    n = len(sets)
    acc = [0.0] * n
    for _ in range(reps):
        order = list(range(n)); rng.shuffle(order)
        seen = set()
        for i, j in enumerate(order):
            seen |= sets[j]; acc[i] += len(seen)
    return [a / reps for a in acc]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--new", required=True, help="directory of *.feature.tbl")
    ap.add_argument("--old", required=True, help="TSV: genome_id<TAB>product")
    ap.add_argument("--out", required=True)
    ap.add_argument("--replicates", type=int, default=100)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    new = {}
    for t in sorted(glob.glob(os.path.join(args.new, "*.feature.tbl"))):
        if os.path.getsize(t) == 0:
            continue
        g = os.path.basename(t)[:-12]; s = set()
        for line in open(t, errors="replace"):
            c = line.rstrip("\n").split("\t")
            #  feature.tbl columns: 6 symbol, 7 start, 8 stop, 9 strand,
            #  11 module, 16 pssm, 17 annotation. The annotation is 17 -- 9 is
            #  the strand, which yields a "vocabulary" of exactly two.
            if len(c) >= 18 and c[17].strip():
                s.add(c[17].strip())
        if s:
            new[g] = s

    old = collections.defaultdict(set)
    for line in open(args.old, errors="replace"):
        c = line.rstrip("\n").split("\t")
        if len(c) >= 2 and c[1].strip() and not c[0].startswith("genome"):
            old[c[0].strip()].add(c[1].strip())

    shared = sorted(set(new) & set(old))
    print("  %d genome(s) annotated by both" % len(shared))
    if len(shared) < 20:
        print("  ERROR: too few shared genomes to rarefy", file=sys.stderr)
        return 1
    ns = [new[g] for g in shared]
    os_ = [old[g] for g in shared]
    print("  distinct strings, all genomes pooled:")
    print("    source vocabulary : %d" % len(set().union(*os_)))
    print("    module vocabulary : %d" % len(set().union(*ns)))

    rng = random.Random(args.seed)
    cn = curve(ns, args.replicates, rng)
    rng = random.Random(args.seed)
    co = curve(os_, args.replicates, rng)

    #  where each curve reaches 95% of its own final value
    def knee(c):
        t = 0.95 * c[-1]
        return next(i + 1 for i, v in enumerate(c) if v >= t)

    n = len(shared)
    print("\n  genomes  source  module")
    for f in (0.01, 0.05, 0.1, 0.25, 0.5, 1.0):
        i = max(0, int(f * n) - 1)
        print("  %7d  %6.1f  %6.1f" % (i + 1, co[i], cn[i]))
    print("\n  95%% of final vocabulary reached at:")
    print("    source : %d genome(s) (%.0f%% of the set)" % (knee(co), 100.0 * knee(co) / n))
    print("    module : %d genome(s) (%.0f%% of the set)" % (knee(cn), 100.0 * knee(cn) / n))
    print("\n  strings per 100 genomes at the end:")
    print("    source : %.1f" % (100.0 * co[-1] / n))
    print("    module : %.1f" % (100.0 * cn[-1] / n))
    print("    %.1f-fold collapse" % (co[-1] / cn[-1]))

    json.dump({"n_genomes": n, "replicates": args.replicates,
               "source": co, "module": cn,
               "source_total": len(set().union(*os_)),
               "module_total": len(set().union(*ns)),
               "source_knee": knee(co), "module_knee": knee(cn)},
              open(args.out, "w"))
    print("\n  wrote %s" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
