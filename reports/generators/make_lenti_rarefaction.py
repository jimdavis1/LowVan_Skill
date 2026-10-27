#!/usr/bin/env python3
"""Rarefaction inputs for the four lentivirus modules, from work already done.

annotation_rarefaction.py wants a directory of *.feature.tbl (it reads only
column 17, the annotation) and a TSV of genome_id<TAB>deposited_product. The
26-27 Oct re-annotation already produced the module's call for every exemplar,
and the L.* dump already holds every deposited product, so neither side needs a
fresh annotation run. The 22 Oct estimate of "~11 real annotation runs each,
16-20 h" that got these pages skipped was costing the annotation twice.

The tables written here carry only the columns the consumer reads, padded to
the real width, and are kept out of the module workdirs so nothing mistakes
them for pipeline output.
"""
import json, os, glob, collections, sys

L = os.path.expanduser("~/lowvan_scratch/lenti")
E = os.path.join(L, "audit_2026-10-26e")
OUT = os.path.join(L, "pagegen", "rare")
WORK = {"HIV2": ("cov_hiv2_work", os.path.join(L, "hiv2")),
        "FIV":  ("cov_genus_work", os.path.join(L, "genus", "FIV")),
        "SRLV": ("cov_genus_work", os.path.join(L, "genus", "SRLV")),
        "EIAV": ("cov_genus_work", os.path.join(L, "genus", "EIAV"))}

#  Which taxon prefixes belong to each module, read from its own collections so
#  this cannot drift from the build.
def module_taxa(wdir, mod):
    taxa = set()
    for f in glob.glob(os.path.join(wdir, "collections", mod, "*.fasta")):
        if os.path.basename(f)[:-6].endswith(".outliers"):
            continue
        for l in open(f):
            if l.startswith(">"):
                h = l[1:].strip().split()[0]
                gid = h.split("|", 1)[1] if h.startswith("fig|") else h.split("|")[0]
                if gid.split(".")[0]:
                    taxa.add(gid.split(".")[0])
    return taxa

def main():
    #  deposited products, per genome, for the whole family dump
    md5s = [l.strip() for l in open(os.path.join(L, "L.uniq.md5"))]
    ann = [l.rstrip("\n").split("\t") for l in open(os.path.join(L, "L.uniq.id_ann"), errors="replace")]
    m2a = {m: " ".join((a[1] if len(a) > 1 else "").split()) for m, a in zip(md5s, ann)}
    dep = collections.defaultdict(list)
    for l in open(os.path.join(L, "L.id_md5"), errors="replace"):
        f = l.rstrip("\n").split("\t")
        if len(f) < 2: continue
        g = ".".join(f[0].split("|")[-1].split(".")[:2])
        s = m2a.get(f[1], "")
        if s: dep[g].append(s)

    for mod, (work, wdir) in WORK.items():
        taxa = module_taxa(wdir, mod)
        tdir = os.path.join(OUT, mod, "tbl")
        os.makedirs(tdir, exist_ok=True)
        for p in glob.glob(os.path.join(tdir, "*.feature.tbl")): os.remove(p)

        n_tbl = 0
        for p in sorted(glob.glob(os.path.join(E, work, "Anno_GTOs", "*.anno.gto"))):
            g = os.path.basename(p)[:-len(".anno.gto")]
            if g.split(".")[0] not in taxa: continue
            try: d = json.load(open(p))
            except Exception: continue
            lines = []
            for ft in d.get("features") or []:
                fam = (ft.get("family_assignments") or [[None, ""]])[0]
                key = str(fam[1] if len(fam) > 1 else "")
                if not key.startswith(mod + "."): continue
                fn = (ft.get("function") or "").strip()
                if not fn: continue
                #  pad to the real table width; the consumer reads col 17 only
                row = [""] * 18
                row[5] = ft.get("id", ""); row[11] = mod; row[16] = key; row[17] = fn
                lines.append("\t".join(row))
            if lines:
                open(os.path.join(tdir, "%s.feature.tbl" % g), "w").write("\n".join(lines) + "\n")
                n_tbl += 1

        op = os.path.join(OUT, mod, "source_products.tsv")
        rows = 0
        with open(op, "w") as fh:
            for t in sorted(glob.glob(os.path.join(tdir, "*.feature.tbl"))):
                g = os.path.basename(t)[:-len(".feature.tbl")]
                for s in dep.get(g, []):
                    fh.write("%s\t%s\n" % (g, s)); rows += 1
        print("  %-5s %3d feature tables, %5d deposited product rows, taxa %s"
              % (mod, n_tbl, rows, ",".join(sorted(taxa))))

main()
