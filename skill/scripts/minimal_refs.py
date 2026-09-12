#!/usr/bin/env python3
"""Propose the fewest reference contigs that close a routing gap.

Runtime routing is a BLASTn against Viral-Rep-Contigs before any profile is
consulted, so a genome too divergent from every reference is never annotated no
matter how good the profiles are. The fix is more references -- but each one is
a permanent maintenance cost, and the clusters are steeply long-tailed, so a
handful of well-chosen contigs recover most of the gap and the next hundred
recover almost nothing.

Greedy set cover over the unrouted genomes: repeatedly take the exemplar that
routes the most still-unrouted genomes, and report the curve so the cutoff is
visible rather than assumed. In Rhabdoviridae the first 5 took coverage from
75% to 91% and the 6th added 0.4%.

    python3 minimal_refs.py --unrouted unrouted.ids --clusters clu_cluster.tsv \
            --json module.json --metadata genome_metadata.tsv --max 10

TWO RULES DECIDE THE FINAL LIST, AND ONLY ONE IS ARITHMETIC.

The script applies the first: a candidate is admissible only if its taxon is one
the module has profiles for. A reference contig for a taxon with no PSSMs routes
genomes to a module that then calls nothing, which is worse than not routing
them -- it converts a visible gap into a silent one.

You decide the second. The knee is a judgement call about how many references
are reasonable to carry, and the counts below inform it; they do not make it.
Stop where the curve flattens.
"""
import argparse, collections, json, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--unrouted", required=True, help="one genome id per line")
    ap.add_argument("--clusters", required=True, help="mmseqs TSV: exemplar<TAB>member")
    ap.add_argument("--json", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--id-col", default="genome.genome_id")
    ap.add_argument("--taxon-col", default="genome.genome_name")
    ap.add_argument("--max", type=int, default=10)
    args = ap.parse_args()

    un = {l.strip() for l in open(args.unrouted) if l.strip()}
    print("  %d unrouted genome(s)" % len(un))

    mem = collections.defaultdict(set)
    for line in open(args.clusters, errors="replace"):
        c = line.rstrip("\n").split("\t")
        if len(c) >= 2:
            mem[c[0]].add(c[1])

    hdr, tax = None, {}
    for line in open(args.metadata, errors="replace"):
        c = line.rstrip("\n").split("\t")
        if hdr is None:
            hdr = c; continue
        r = dict(zip(hdr, c))
        tax[r.get(args.id_col, "").strip()] = r.get(args.taxon_col, "").strip()

    mod = json.load(open(args.json))
    covered = [m for m in mod if isinstance(mod[m], dict) and mod[m].get("features")]
    print("  modules with profiles: %s" % ", ".join(sorted(covered)))

    def admissible(ex):
        """True when the exemplar's taxon is one a module actually covers."""
        name = tax.get(ex, "").lower()
        return any(m.lower().rstrip("es").rstrip("s") in name
                   or any(w.lower() in name for w in (m,)) for m in covered)

    cand = {e: (v & un) for e, v in mem.items() if v & un}
    print("  %d cluster(s) touch the gap\n" % len(cand))

    left, chosen, skipped = set(un), [], []
    for _ in range(args.max):
        best, gain = None, 0
        for e, v in cand.items():
            g = len(v & left)
            if g > gain:
                best, gain = e, g
        if not best or not gain:
            break
        if not admissible(best):
            skipped.append((best, gain, tax.get(best, "?")))
            del cand[best]; continue
        left -= cand[best]; chosen.append((best, gain, tax.get(best, "?")))
        del cand[best]

    n = len(un)
    print("  rank  exemplar          adds   cumulative  taxon")
    cum = 0
    for i, (e, g, t) in enumerate(chosen, 1):
        cum += g
        print("  %4d  %-16s %5d   %5.1f%%     %s"
              % (i, e, g, 100.0 * cum / n, t[:44]))
    if skipped:
        print("\n  excluded -- no module covers the taxon (would route to an empty call):")
        for e, g, t in skipped[:8]:
            print("    %-16s would have added %4d   %s" % (e, g, t[:44]))

    if len(chosen) > 1:
        print("\n  marginal gain per added reference:")
        for i, (e, g, _t) in enumerate(chosen, 1):
            print("    #%-3d %4d genome(s)  %s" % (i, g, "#" * min(52, g)))
        print("\n  Stop where this flattens. That cutoff is yours to set, not the")
        print("  script's -- weigh the genomes recovered against carrying one more")
        print("  reference contig forever.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
