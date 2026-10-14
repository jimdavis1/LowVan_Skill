#!/usr/bin/env python3
"""Drop the sequences that make a collection's N-terminus ragged.

A CDS profile has to begin where the protein begins, or the annotator's hit
starts late and `scan_to_met_start` then halts at whatever internal Met it
meets first. What makes that happen is not bad homology -- it is a minority of
collection members whose N-terminus disagrees with the rest, which drops the
occupancy of the first alignment columns below the curation threshold and gets
them trimmed off. The majority then loses its own start.

Norovirus POLY: 46% of 3,411 sequences were N-terminally short, most of them by
ONE OR TWO residues, and the largest cluster was trimmed at column 349. The
polyprotein was called a median 127 codons late on 409 genomes while every
mature peptide inside it tiled correctly. Sapovirus POLY failed the same way on
a 71-member cluster, caught by `check_met_start.py` before it shipped.

Two independent tests, because they catch different populations:

  EXTENDED   another sequence in the collection extends this one's N-terminus.
             That is the definition of a ragged start and it is self-contained:
             no reference, no canonical, no genogroup assumption. Do NOT test
             against a single canonical start -- that is genogroup-specific and
             silently discards divergent lineages, the same error as writing a
             cleavage motif from one strain.

  SINGLETON  this sequence's N-terminus resembles no other sequence's. These
             are not truncations of anything, so EXTENDED cannot see them:
             mis-annotated frames, wrong-ORF calls, and genuinely unique starts.
             A handful of them in a cluster is enough to collapse the column
             occupancy for everyone else -- three such members cost 20 Sapovirus
             POLY sequences 85 residues each.

    python3 filter_nterm.py --fasta collections/M/POLY.fasta          # report
    python3 filter_nterm.py --fasta collections/M/POLY.fasta --write

`--min-support` is the one judgement. At 2 (the default) a start must be shared
with two other sequences to survive, so a lineage represented three times is
kept and a one-off is not. Raise it for a large collection, lower it only if
the taxon is genuinely sparse -- and check what comes out, because this throws
away real sequence to buy a clean terminus.
"""
import argparse, collections, os, shutil, sys


def read_fasta(path):
    out, h, s = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if h is not None:
                out.append((h, "".join(s)))
            h, s = line.rstrip("\n"), []
        else:
            s.append(line.strip())
    if h is not None:
        out.append((h, "".join(s)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--k", type=int, default=18, help="N-terminal k-mer (default 18)")
    ap.add_argument("--max-offset", type=int, default=30,
                    help="how far in to look for an extending match (default 30)")
    ap.add_argument("--support-k", type=int, default=10,
                    help="k-mer for the singleton test (default 10)")
    ap.add_argument("--max-mismatch", type=int, default=3,
                    help="mismatches allowed when matching starts (default 3)")
    ap.add_argument("--min-support", type=int, default=2,
                    help="other sequences that must share a start (default 2)")
    ap.add_argument("--no-singleton", action="store_true",
                    help="run only the EXTENDED test")
    ap.add_argument("--require-met", action="store_true",
                    help="also drop every sequence that does not begin with Met")
    ap.add_argument("--write", action="store_true",
                    help="rewrite the fasta, keeping a .preNterm backup")
    a = ap.parse_args()

    recs = read_fasta(a.fasta)
    if not recs:
        sys.exit("no sequences in %s" % a.fasta)

    #  EXTENDED: is this sequence's opening k-mer found at an INTERNAL offset
    #  of some other sequence? If so, that other sequence extends it.
    internal = collections.defaultdict(set)
    for _, x in recs:
        for k in range(1, a.max_offset + 1):
            if len(x) >= k + a.k:
                internal[x[k:k + a.k]].add(k)

    extended, offs = set(), collections.Counter()
    for i, (h, x) in enumerate(recs):
        if len(x) < a.k:
            continue
        hit = internal.get(x[:a.k])
        if hit:
            extended.add(i)
            offs[min(hit)] += 1

    #  SINGLETON: does this sequence's start resemble any other's?
    singleton = set()
    if not a.no_singleton:
        starts = [(i, x[:a.support_k]) for i, (_, x) in enumerate(recs)
                  if len(x) >= a.support_k]
        arr = [s for _, s in starts]
        for n, (i, s) in enumerate(starts):
            sup = 0
            for m, t in enumerate(arr):
                if m == n:
                    continue
                d = 0
                for p, q in zip(s, t):
                    if p != q:
                        d += 1
                        if d > a.max_mismatch:
                            break
                if d <= a.max_mismatch:
                    sup += 1
                    if sup >= a.min_support:
                        break
            if sup < a.min_support:
                singleton.add(i)

    #  A third test, off by default and needed where the full-length form is a
    #  MINORITY of the deposits. EXTENDED only catches a truncation whose longer
    #  counterpart is also in the collection; where most records are short, the
    #  short ones have nothing to be caught by and survive. Lagovirus POLY is
    #  that case: 747 records, only 93 beginning at the real start, and after
    #  EXTENDED + SINGLETON the collection was still 63% non-Met. Restricting to
    #  Met left 126 sequences with a median length of exactly 2,344 -- the
    #  literature ORF1 -- still covering both species.
    #
    #  Only for a CDS, and only after checking what it removes: a protein with a
    #  genuine non-AUG start loses everything here.
    nomet = set()
    if a.require_met:
        nomet = {i for i, (_, x) in enumerate(recs) if not x.upper().startswith("M")}

    drop = extended | singleton | nomet
    keep = [r for i, r in enumerate(recs) if i not in drop]

    print("  %s" % a.fasta)
    print("    sequences in                 %6d" % len(recs))
    print("    EXTENDED  (another extends)  %6d  (%.1f%%)"
          % (len(extended), 100.0 * len(extended) / len(recs)))
    if not a.no_singleton:
        print("    SINGLETON (start unique)     %6d  (%.1f%%)"
              % (len(singleton), 100.0 * len(singleton) / len(recs)))
    if a.require_met:
        print("    NO-MET    (does not start M) %6d  (%.1f%%)"
              % (len(nomet), 100.0 * len(nomet) / len(recs)))
    print("    kept                         %6d  (%.1f%%)"
          % (len(keep), 100.0 * len(keep) / len(recs)))
    if offs:
        shown = ", ".join("-%d:%d" % (k, v) for k, v in sorted(offs.items())[:6])
        print("    residues missing, commonest: %s" % shown)
    if len(keep) < 0.25 * len(recs):
        print("    WARNING: this drops more than three quarters of the collection."
              " Check the taxon before writing.")
    if not keep:
        sys.exit("    refusing to leave an empty collection")

    if a.write:
        shutil.copy2(a.fasta, a.fasta + ".preNterm")
        with open(a.fasta, "w") as fh:
            for h, x in keep:
                fh.write("%s\n%s\n" % (h, x))
        print("    wrote %d sequences (backup at %s.preNterm)"
              % (len(keep), os.path.basename(a.fasta)))
    else:
        print("    (report only -- rerun with --write)")


main()
