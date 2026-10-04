#!/usr/bin/env python3
"""Measure how fast the vocabulary of annotation strings grows.

A controlled vocabulary and an uncontrolled one both label every protein. The
difference only becomes visible when you count DISTINCT strings as proteins
accumulate: a controlled vocabulary saturates, because the Nth protein reuses
a name the first N-1 already established, while free text keeps growing,
because every submitter spells the same protein a new way.

The experiment is a rarefaction over PROTEINS. Pool every annotated protein --
one record per called feature in the module's tables, one per product in the
source -- shuffle the pool, draw the proteins one at a time, and count how
many distinct strings have been seen after each draw. Average over replicates.
Both curves start at (0, 0): no string can have been seen before the first
protein is drawn, and k proteins can carry at most k strings.

An earlier version shuffled GENOMES and added a whole genome's set of strings
per step, so the curve started at that genome's vocabulary size (13 or 15
rather than 0) and its x-axis counted genomes. That is not a rarefaction of
annotations, and every page built from it was wrong.

    python3 annotation_rarefaction.py --new coverage_eval/ann \\
            --old bvbrc_annotations.tsv --out rarefaction.json --replicates 100

Only genomes present in both sets are used, so the two pools describe the
same genomes. The pools differ in size -- the source often records one
polyprotein where the module calls a dozen mature peptides -- so the ratio is
reported at the depth both curves reach, as well as at each curve's end.
"""
import argparse, collections, glob, json, os, random, sys

NPTS = 400   # points kept per curve; the curves are smooth at this density


def curve(pool, reps, rng):
    """Mean distinct-string count after k draws, k = 0..len(pool)."""
    n = len(pool)
    acc = [0.0] * (n + 1)
    idx = {s: i for i, s in enumerate(sorted(set(pool)))}
    codes = [idx[s] for s in pool]
    for _ in range(reps):
        rng.shuffle(codes)
        seen = bytearray(len(idx)); d = 0
        for k, c in enumerate(codes, 1):
            if not seen[c]:
                seen[c] = 1; d += 1
            acc[k] += d
    return [a / reps for a in acc]


def thin(c):
    n = len(c) - 1
    xs = sorted({round(n * i / (NPTS - 1)) for i in range(NPTS)} | {0, n})
    return {"n": n, "x": xs, "y": [round(c[x], 3) for x in xs]}


def knee(c):
    t = 0.95 * c[-1]
    return next(k for k, v in enumerate(c) if v >= t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--new", required=True, help="directory of *.feature.tbl")
    ap.add_argument("--old", required=True, help="TSV: genome_id<TAB>product, one row per protein")
    ap.add_argument("--out", required=True)
    ap.add_argument("--replicates", type=int, default=100)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    new = collections.defaultdict(list)
    for t in sorted(glob.glob(os.path.join(args.new, "*.feature.tbl"))):
        g = os.path.basename(t)[:-12]
        for line in open(t, errors="replace"):
            c = line.rstrip("\n").split("\t")
            #  feature.tbl columns: 6 symbol, 7 start, 8 stop, 9 strand,
            #  11 module, 16 pssm, 17 annotation. The annotation is 17.
            if len(c) >= 18 and c[17].strip():
                new[g].append(c[17].strip())

    old = collections.defaultdict(list)
    for line in open(args.old, errors="replace"):
        c = line.rstrip("\n").split("\t")
        if len(c) >= 2 and c[1].strip() and not c[0].startswith("genome"):
            old[c[0].strip()].append(c[1].strip())

    shared = sorted(set(new) & set(old))
    print("  %d genome(s) annotated by both" % len(shared))
    if len(shared) < 20:
        print("  ERROR: too few shared genomes to rarefy", file=sys.stderr)
        return 1
    npool = [s for g in shared for s in new[g]]
    opool = [s for g in shared for s in old[g]]
    print("  proteins: source %d, module %d" % (len(opool), len(npool)))

    co = curve(opool, args.replicates, random.Random(args.seed))
    cn = curve(npool, args.replicates, random.Random(args.seed))
    m = min(len(opool), len(npool))

    print("\n  proteins  source  module")
    for f in (0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0):
        k = max(1, int(f * m))
        print("  %8d  %6.1f  %6.1f" % (k, co[k], cn[k]))
    print("\n  distinct strings, all proteins pooled: source %d, module %d"
          % (len(set(opool)), len(set(npool))))
    print("  95%% of final vocabulary reached after: source %d protein(s), module %d"
          % (knee(co), knee(cn)))
    print("  at the depth both reach (%d proteins): source %.1f, module %.1f -> %.1f-fold"
          % (m, co[m], cn[m], co[m] / cn[m]))

    json.dump({"unit": "proteins", "n_genomes": len(shared), "replicates": args.replicates,
               "source": thin(co), "module": thin(cn),
               "source_total": len(set(opool)), "module_total": len(set(npool)),
               "source_knee": knee(co), "module_knee": knee(cn),
               "common_n": m, "source_at_common": round(co[m], 3),
               "module_at_common": round(cn[m], 3)},
              open(args.out, "w"))
    print("\n  wrote %s" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
