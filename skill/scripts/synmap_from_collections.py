#!/usr/bin/env python3
"""Reconstruct synonyms.tsv from a dump and the collections already built.

For modules built before collection_engine.py emitted synonyms.tsv. It does
NOT re-run the binning -- it reads which sequences ended up in which
collection and recovers the annotation strings behind them, so the shipped
collections are left exactly as they are.

That distinction matters: re-running build_collections.py to get the file
would truncate the collections and discard the homology rescue, which is the
trap this project keeps re-tripping.

    python3 synmap_from_collections.py --workdir Closterovirus \
            --module Closterovirus --dump Clos

Writes <workdir>/collections/synonyms.tsv with count/annotation/genus/bins_to,
counting FEATURES (via id_md5) rather than unique sequences, the same as the
engine.

ONE DIFFERENCE FROM THE ENGINE, AND IT CHANGES THE NUMBERS. The engine
records a string as bound when a RULE matched it, before the length,
ambiguity and cleavage filters run. This script can only see what survived
into the collection, so a sequence whose string matched CP but which was then
dropped as too short appears here as UNASSIGNED.

So the two are not interchangeable:

    engine        "this string names this protein"       vocabulary
    this script   "this string reached this collection"  post-QC membership

Both are honest and neither is wrong, but a page built from this file reports
a lower bound fraction than one built from the engine's, and must say so.
Do not compare the two directly across modules.
"""
import argparse, collections, glob, os, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--module", required=True)
    ap.add_argument("--dump", required=True)
    a = ap.parse_args()
    W, M, T = a.workdir, a.module, a.dump

    md5s = [l.strip() for l in open(os.path.join(W, "%s.uniq.md5" % T))]
    ann = [l.rstrip("\n").split("\t") for l in open(os.path.join(W, "%s.uniq.id_ann" % T))]
    if len(md5s) != len(ann):
        sys.exit("uniq.md5 and uniq.id_ann differ in length -- run check_dump.py")

    per_md5 = collections.Counter()
    for l in open(os.path.join(W, "%s.id_md5" % T), errors="replace"):
        p = l.rstrip("\n").split("\t")
        if len(p) >= 2 and p[1].strip():
            per_md5[p[1].strip()] += 1

    g2genus = {}
    gp = os.path.join(W, "%s.id_name_gs" % T)
    if os.path.exists(gp):
        for l in open(gp, errors="replace"):
            f = l.rstrip("\n").split("\t")
            if len(f) >= 4:
                g2genus[f[0]] = f[3] or "unclassified"

    def genus_of(fid):
        try:
            return g2genus.get(".".join(fid.split("|", 1)[1].split(".")[:2]), "unclassified")
        except IndexError:
            return "unclassified"

    #  feature id -> the collection it landed in
    where = {}
    cdir = os.path.join(W, "collections", M)
    for f in sorted(glob.glob(os.path.join(cdir, "*.fasta"))):
        k = os.path.basename(f)[:-6]
        if k.endswith(".outliers"):
            continue
        for l in open(f):
            if l.startswith(">"):
                where[l[1:].strip().split()[0]] = k

    syn = collections.Counter()
    for i, m in enumerate(md5s):
        fid = ann[i][0] if ann[i] else ""
        s = " ".join((ann[i][1] if len(ann[i]) > 1 else "").split())
        if not fid:
            continue
        syn[(s, genus_of(fid), where.get(fid, "UNASSIGNED"))] += per_md5.get(m, 1)

    out = os.path.join(W, "collections", "synonyms.tsv")
    with open(out, "w") as fh:
        fh.write("count\tannotation\tgenus\tbins_to\n")
        for (s, g, k), c in sorted(syn.items(), key=lambda x: -x[1]):
            fh.write("%d\t%s\t%s\t%s\n" % (c, s, g, k))
    binned = sum(c for (_s, _g, k), c in syn.items() if k != "UNASSIGNED")
    tot = sum(syn.values())
    print("  %-14s %d rows, %d of %d features binned = %.1f%% -> %s"
          % (M, len(syn), binned, tot, 100.0 * binned / tot if tot else 0.0, out))


main()
