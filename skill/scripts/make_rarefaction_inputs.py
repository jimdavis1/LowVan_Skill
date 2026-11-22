#!/usr/bin/env python3
"""Build the two inputs annotation_rarefaction.py wants, for a run that
produced GTOs rather than a flat ann/ directory.

  --new : a directory of <genome>.feature.tbl.  The rarefaction parser reads
          only the filename and column 17 (the annotation string), so the rows
          written here carry the real annotation in that column and placeholders
          elsewhere.  The annotator REWRITES the gto's internal id, so the
          genome id comes from the FILENAME, which is the only link back to the
          deposit.

  --old : genome_id <TAB> product, one row per deposited protein, recovered by
          joining the dump's three line-index-aligned uniq.* files.  That
          alignment IS the join key -- the files are never sorted independently.

Usage:
  make_rarefaction_inputs.py <Anno_GTOs> <outdir> <module> <dumpdir> <dumpstem>
"""
import json
import os
import sys

gdir, outdir, MODULE, dumpdir, stem = sys.argv[1:6]
tbl = os.path.join(outdir, "tbl")
os.makedirs(tbl, exist_ok=True)

#  ---- new: one .feature.tbl per genome -------------------------------------
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
        fam = ""
        for fa in f.get("family_assignments", []) or []:
            if fa and fa[0] == "LOWVAN":
                fam = str(fa[1])
        if not fam.startswith(MODULE + "."):
            continue
        anno = (f.get("function") or "").strip()
        if not anno:
            continue
        loc = f.get("location") or [["", 0, "+", 0]]
        c, start, strand, ln = loc[0][0], loc[0][1], loc[0][2], loc[0][3]
        rows.append("\t".join([
            gid, "Viruses", gid, "LV", f.get("type", "CDS"),
            f.get("id", ""), "", str(start), str(start), strand, str(ln),
            MODULE, "", "", "", g.get("scientific_name", ""),
            fam + ".pssm", anno]))
    if not rows:
        continue
    open(os.path.join(tbl, gid + ".feature.tbl"), "w").write("\n".join(rows) + "\n")
    n_g += 1
    n_f += len(rows)
print("new: %d genome(s), %d feature row(s) -> %s" % (n_g, n_f, tbl))

#  ---- old: genome_id <TAB> deposited product -------------------------------
md5s = [l.strip() for l in open(os.path.join(dumpdir, stem + ".uniq.md5"))]
anns = [l.rstrip("\n").split("\t") for l in
        open(os.path.join(dumpdir, stem + ".uniq.id_ann"))]
if len(md5s) != len(anns):
    sys.exit("uniq.md5 (%d) and uniq.id_ann (%d) are NOT line-aligned -- the "
             "join key is broken, refusing to guess" % (len(md5s), len(anns)))
ann_of = {m: (a[1] if len(a) > 1 else "") for m, a in zip(md5s, anns)}

out = os.path.join(outdir, "old_annotations.tsv")
n = miss = 0
with open(out, "w") as fh:
    for line in open(os.path.join(dumpdir, stem + ".id_md5")):
        c = line.rstrip("\n").split("\t")
        if len(c) < 2:
            continue
        fid, md5 = c[0], c[1]
        a = ann_of.get(md5)
        if not a:
            miss += 1
            continue
        #  fig|28344.26516.peg.3  ->  28344.26516
        try:
            gid = ".".join(fid.split("|")[1].split(".")[:2])
        except IndexError:
            continue
        fh.write("%s\t%s\n" % (gid, a))
        n += 1
print("old: %d row(s), %d md5 without an annotation -> %s" % (n, miss, out))
