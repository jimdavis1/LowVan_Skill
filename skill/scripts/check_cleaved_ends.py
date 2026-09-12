#!/usr/bin/env python3
"""Check that cleaved feature ends are crisp and have not crept.

Where a feature declares upstream_ext 0 or downstream_ext 0, the annotator does
not scan to a start or stop codon: the end is wherever the tblastn alignment
stopped. That is correct for a protease cut site, whose position is not marked
by a codon -- and silently wrong for any terminus that is NOT a cut site, where
the boundary then drifts with the alignment.

The cost is easy to miss and compounds. Alpharhabdovirinae G_MAT carried
downstream_ext 0, implying its C-terminus was a cut site. It is not: signal
peptidase removes the leader and nothing else, so the mature chain runs to the
precursor's own stop codon. Measured over 1108 annotated genomes, 45% of mature
G calls ended somewhere other than the precursor's C-terminus, median 28 codons
short and up to 321.

This matters far more for polyproteins. A precursor with three internal products
has six cleaved termini, and every one left to drift compounds along the chain.

    python3 check_cleaved_ends.py --json module.json --ann-dir coverage_eval/ann \
            --precursor G --products G_SP,G_MAT

Reports, per adjacent pair and per product-vs-precursor boundary, the
distribution of offsets in nucleotides. A clean cut is a spike at one value; a
spread is creep. An offset of +3 against a precursor end is the stop codon and
is expected.
"""
import argparse, collections, glob, json, os, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True)
    ap.add_argument("--ann-dir", required=True)
    ap.add_argument("--precursor", required=True)
    ap.add_argument("--products", required=True, help="comma-separated, 5' to 3'")
    ap.add_argument("--tolerance", type=int, default=3,
                    help="nt offsets within this of the mode count as crisp")
    args = ap.parse_args()

    prods = [p.strip() for p in args.products.split(",") if p.strip()]

    def norm(x):
        x = "".join(c for c in x.upper() if c.isalnum())
        return x.replace("MATURE", "MAT").replace("SIGNAL", "SP")

    def same(a, b):
        """The JSON key and the table's symbol are often spelled differently:
        the key G_MAT appears in output as G_mature, and G_SP as SP.

        Match whole underscore segments only. A bare suffix test would make the
        key P match the product SP."""
        if norm(a) == norm(b):
            return True
        return (norm(a.split("_")[-1]) == norm(b)
                or norm(b.split("_")[-1]) == norm(a))

    mod = json.load(open(args.json))
    decl = {}
    for m, b in mod.items():
        if not isinstance(b, dict):
            continue
        for k, e in (b.get("features") or {}).items():
            for want in prods + [args.precursor]:
                if same(k, want):
                    decl.setdefault(want, []).append(
                        (m, k, e.get("upstream_ext"), e.get("downstream_ext")))
                    break
    print("  declared extensions")
    for want in [args.precursor] + prods:
        if not decl.get(want):
            print("    %-22s no JSON feature resolves to %r" % ("", want))
        for m, k, u, d in decl.get(want, []):
            print("    %-22s %-10s upstream=%s%s downstream=%s%s"
                  % (m, k, u, " (cleaved)" if u == 0 else "",
                     d, " (cleaved)" if d == 0 else ""))

    #  gather calls; symbol column (9) is what appears in the table
    rows = []
    for t in sorted(glob.glob(os.path.join(args.ann_dir, "*.feature.tbl"))):
        if os.path.getsize(t) == 0:
            continue
        f = collections.defaultdict(list)
        for line in open(t, errors="replace"):
            c = line.rstrip("\n").split("\t")
            if len(c) < 18:
                continue
            #  a no_features_called row carries the routing decision (module,
            #  reference contig, bitscore) with empty coordinates, so that a
            #  genome which routed but called nothing still leaves a record.
            #  It is not a feature; skip it.
            try:
                lo, hi = int(c[7]), int(c[8])
            except ValueError:
                continue
            f[c[6]].append((min(lo, hi), max(lo, hi)))
        rows.append(f)
    print("\n  %d annotated genome(s)" % len(rows))

    def report(name, counts, note):
        tot = sum(counts.values())
        if not tot:
            print("\n  %s: no genomes carry both" % name)
            return
        mode, n = counts.most_common(1)[0]
        crisp = sum(v for k, v in counts.items() if abs(k - mode) <= args.tolerance)
        print("\n  %s   n=%d   %s" % (name, tot, note))
        print("    mode %+d nt in %d (%.0f%%); within +/-%d of it: %d (%.0f%%)"
              % (mode, n, 100 * n / tot, args.tolerance, crisp, 100 * crisp / tot))
        drift = tot - crisp
        if drift:
            far = sorted(k for k in counts if abs(k - mode) > args.tolerance)
            print("    CREEP: %d genome(s) (%.0f%%) outside that, offsets %d..%d nt"
                  % (drift, 100 * drift / tot, far[0], far[-1]))
        for k, v in sorted(counts.items())[:6]:
            print("      %+6d nt %5d %s" % (k, v, "#" * int(38 * v / tot)))

    #  each product's start against the previous boundary
    prev = args.precursor
    for p in prods:
        c = collections.Counter()
        for f in rows:
            if p not in f or prev not in f:
                continue
            if prev == args.precursor:
                c[f[p][0][0] - f[prev][0][0]] += 1
            else:
                c[f[p][0][0] - f[prev][0][1] - 1] += 1
        report("%s start vs %s" % (p, "precursor start" if prev == args.precursor else prev + " end"),
               c, "0 = exactly adjacent")
        prev = p

    #  the last product's end against the precursor's
    last = prods[-1]
    c = collections.Counter()
    for f in rows:
        if last in f and args.precursor in f:
            c[f[args.precursor][0][1] - f[last][0][1]] += 1
    report("%s end vs precursor end" % last, c,
           "0 or +3 expected (+3 = the precursor's stop codon)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
