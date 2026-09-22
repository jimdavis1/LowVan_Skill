#!/usr/bin/env python3
"""Build the two inputs annotation_rarefaction.py needs.

  tbl/<genome_id>.feature.tbl   what this module called, one row per feature
  bvbrc_products.tsv            genome_id<TAB>product, what BV-BRC had before

The rarefaction reads only column 17 (the annotation) and the filename, but
the rows are written with all 18 columns so the files are the same shape as
the ones evaluate_coverage.py emits and can be read by the same tools.

The source side needs care: uniq.id_ann holds one annotation per UNIQUE
protein, not per feature, so it has to be joined back through id_md5 to
recover what each genome was actually annotated with. Reading uniq.id_ann
directly would give each genome the vocabulary of its representative
sequences only, which understates the source vocabulary -- the very quantity
the artifact is measuring.
"""
import collections, glob, json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
T = "Ortho"
MODULE = "Orthoflavivirus"

decl = json.load(open(os.path.join(HERE, "%s_Viral_PSSM.json" % MODULE)))[MODULE]["features"]
anno2sym = {v["anno"]: v["gene_symbol"] for v in decl.values()}

# ---- the module side, from the annotated GTOs -------------------------------
out = os.path.join(HERE, "tbl")
os.makedirs(out, exist_ok=True)
for f in glob.glob(os.path.join(out, "*.feature.tbl")):
    os.remove(f)
n = rows = 0
for f in sorted(glob.glob(os.path.join(HERE, "gto_out", "*.ann.gto"))):
    try:
        d = json.load(open(f))
    except Exception:
        continue
    gid = d.get("source_id") or d.get("id") or os.path.basename(f)[:-8]
    fam = d.get("viral_family", MODULE)
    lines = []
    for i, x in enumerate(d.get("features", []), 1):
        fn = (x.get("function") or "").strip()
        if not fn:
            continue
        loc = (x.get("location") or [["", 0, "+", 0]])[0]
        try:
            start, strand, ln = int(loc[1]), loc[2], int(loc[3])
        except Exception:
            start, strand, ln = 0, "+", 0
        lines.append("\t".join([
            d.get("ncbi_taxonomy_id", "") and str(d["ncbi_taxonomy_id"]) or "",
            "Viruses", str(gid), "LV", x.get("type", "CDS"),
            "lv|%s.%d" % (gid, i), anno2sym.get(fn, ""),
            str(start), str(start + ln), strand, str(ln),
            fam, "", "", "", d.get("scientific_name", ""), "", fn]))
    if lines:
        open(os.path.join(out, "%s.feature.tbl" % gid), "w").write("\n".join(lines) + "\n")
        n += 1
        rows += len(lines)
print("wrote %d feature.tbl file(s), %d rows -> tbl/" % (n, rows))

# ---- the source side, joined back through id_md5 ----------------------------
md5_anno = {}
md5s = [l.strip() for l in open(os.path.join(HERE, "%s.uniq.md5" % T))]
ann = [l.rstrip("\n").split("\t") for l in open(os.path.join(HERE, "%s.uniq.id_ann" % T))]
for i, m in enumerate(md5s):
    a = ann[i][1] if i < len(ann) and len(ann[i]) > 1 else ""
    if m and a:
        md5_anno[m] = " ".join(a.split())

pergen = collections.defaultdict(set)
miss = 0
for l in open(os.path.join(HERE, "%s.id_md5" % T)):
    p = l.rstrip("\n").split("\t")
    if len(p) < 2:
        continue
    fid, m = p[0], p[1].strip()
    g = re.match(r"fig\|(\d+\.\d+)\.", fid)
    if not g:
        continue
    a = md5_anno.get(m)
    if a:
        pergen[g.group(1)].add(a)
    else:
        miss += 1
with open(os.path.join(HERE, "bvbrc_products.tsv"), "w") as fh:
    for g in sorted(pergen):
        for a in sorted(pergen[g]):
            fh.write("%s\t%s\n" % (g, a))
print("wrote bvbrc_products.tsv: %d genomes, %d rows (%d features had no annotation)"
      % (len(pergen), sum(len(v) for v in pergen.values()), miss))

tbl = {os.path.basename(x)[:-12] for x in glob.glob(os.path.join(out, "*.feature.tbl"))}
print("\ngenomes in both: %d" % len(tbl & set(pergen)))
