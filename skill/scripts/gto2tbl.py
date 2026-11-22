#!/usr/bin/env python3
"""Anno_GTOs -> a directory of *.feature.tbl, for the kit's tbl-based QC.

`evaluate_coverage.py` writes GTOs, not the flat `ann/` tree that
`check_cleaved_ends.py` and the rarefaction curve expect, so this bridges
them.  Columns 7 and 8 carry REAL genomic coordinates, 1-based and inclusive,
converted from the GTO's (start, strand, length) triple -- the whole point of
the cleaved-ends check is the offset between one feature's end and the next
one's start, so a placeholder coordinate would produce a confident wrong
answer rather than an error.

Column 6 carries the MODULE'S OWN FEATURE KEY (`NSP1A`), not the deposited
gene symbol (`nsp1alpha`).  check_cleaved_ends matches the table's column 6
against the JSON key by normalising underscores and case, and `nsp1alpha` does
not normalise to `NSP1A` under any of its rules -- so emitting the symbol
makes every product report "no JSON feature resolves to it".

The annotator rewrites a gto's internal id, so the genome id comes from the
FILENAME.

Usage: gto2tbl.py <Anno_GTOs> <outdir> [module]
"""
import json
import os
import sys

gdir, outdir = sys.argv[1], sys.argv[2]
MODULE = sys.argv[3] if len(sys.argv) > 3 else None
os.makedirs(outdir, exist_ok=True)

n_g = n_f = 0
for fn in sorted(os.listdir(gdir)):
    if not fn.endswith(".gto"):
        continue
    gid = fn[:-4]
    if gid.endswith(".anno"):
        gid = gid[:-5]
    try:
        g = json.load(open(os.path.join(gdir, fn)))
    except ValueError:
        continue
    rows = []
    for f in g.get("features", []):
        mod = key = None
        for fa in f.get("family_assignments", []) or []:
            if not fa:
                continue
            #  a PSSM call is ['LOWVAN', '<Module>.<KEY>.<n>', ...]; a
            #  transcript-edited call is ['<Module>', '<KEY>', ...]
            if fa[0] == "LOWVAN" and "." in str(fa[1]):
                p = str(fa[1]).split(".")
                mod, key = p[0], p[1]
            elif len(fa) > 1 and str(fa[-1]).startswith("LowVan Transcript"):
                mod, key = str(fa[0]), str(fa[1])
        if not key or (MODULE and mod != MODULE):
            continue
        loc = (f.get("location") or [None])[0]
        if not loc or len(loc) < 4:
            continue
        contig, start, strand, ln = loc[0], int(loc[1]), loc[2], int(loc[3])
        if strand == "+":
            lo, hi = start, start + ln - 1
        else:
            lo, hi = start - ln + 1, start
            lo, hi = hi, lo          # report 5'->3' as the tbl does
        #  Column order is the annotator's, and it is load-bearing: the
        #  consumers index by position.  0 contig, 1 domain, 2 genome id,
        #  3 source, 4 type, 5 FEATURE ID, 6 symbol, 7 start, 8 stop,
        #  9 strand, 10 length, 11 module, 12 reference contig, 13 contig
        #  length, 14 accession, 15 genome name, 16 pssm, 17 annotation.
        #  Writing the key one column early put the symbol where the start
        #  belongs; check_cleaved_ends then read a coordinate as the feature
        #  name and reported "no genomes carry both" for every pair -- a
        #  clean, confident, entirely wrong answer.
        rows.append("\t".join([
            gid, "Viruses", gid, "LV", f.get("type", "CDS"),
            str(f.get("id", "")), key,
            str(lo), str(hi), strand, str(ln), mod, "", str(len(contig)), "",
            g.get("scientific_name", ""),
            "%s.%s.pssm" % (mod, key), (f.get("function") or "").strip()]))
    if not rows:
        continue
    with open(os.path.join(outdir, gid + ".feature.tbl"), "w") as fh:
        fh.write("\n".join(rows) + "\n")
    n_g += 1
    n_f += len(rows)

print("%d genome(s), %d feature row(s) -> %s" % (n_g, n_f, outdir))
