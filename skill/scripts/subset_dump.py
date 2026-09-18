#!/usr/bin/env python3
"""Cut a per-module dump out of the class-wide Alsuviricetes dump.

Preserves the contract in references/inputs.md, in particular the line-index
alignment of the three uniq.* files -- they are re-emitted in one pass over the
class uniq.md5 so the order can not drift.
"""
import os, sys, argparse

ap = argparse.ArgumentParser()
ap.add_argument("--src",  default=".", help="dir holding the class dump")
ap.add_argument("--prefix", default="Alsu")
ap.add_argument("--out",  required=True, help="module workdir to create")
ap.add_argument("--name", required=True, help="dump name inside that workdir")
ap.add_argument("--families", default="", help="comma-separated family names")
ap.add_argument("--genera",   default="", help="comma-separated genus names")
ap.add_argument("--genome-ids", default="", help="file of genome ids, overrides the above")
a = ap.parse_args()

fams = set(x for x in a.families.split(",") if x)
gens = set(x for x in a.genera.split(",")   if x)
S, P = a.src, a.prefix

keep = set()
rows = []
for l in open(os.path.join(S, P + ".id_name_gs")):
    p = l.rstrip("\n").split("\t")
    gid, fam, gen = p[0], (p[2] if len(p) > 2 else ""), (p[3] if len(p) > 3 else "")
    if a.genome_ids: continue
    if (fams and fam in fams) or (gens and gen in gens):
        keep.add(gid); rows.append(l)
if a.genome_ids:
    want = set(x.strip() for x in open(a.genome_ids) if x.strip())
    rows = []
    for l in open(os.path.join(S, P + ".id_name_gs")):
        gid = l.split("\t")[0]
        if gid in want: keep.add(gid); rows.append(l)

os.makedirs(a.out, exist_ok=True)
O = lambda ext: os.path.join(a.out, a.name + ext)
open(O(".id_name_gs"), "w").writelines(rows)

md5s = set()
n = 0
with open(O(".id_md5"), "w") as out:
    for l in open(os.path.join(S, P + ".id_md5")):
        fid, md5 = l.rstrip("\n").split("\t")
        gid = fid.split("|", 1)[1].rsplit(".", 2)[0] if "|" in fid else ""
        if gid in keep:
            out.write(l); n += 1
            if md5: md5s.add(md5)

#  Re-pick the representative feature id from INSIDE the subset.
#
#  The class-wide uniq.id_ann names one representative feature per md5, and for
#  a genus-level subset that feature often belongs to a genome in a different
#  genus that happens to encode an identical protein. Carrying it over leaves
#  check_dump.py reporting "uniq.id_ann row resolves to None" -- 15 rows on
#  Tobamovirus, whose proteins are shared with other Virgaviridae genera. The
#  annotation string is still the class-wide one for that sequence, which is
#  correct; only the id has to be local.
local_rep = {}
for l in open(O(".id_md5")):
    fid, md5 = l.rstrip("\n").split("\t")
    if md5: local_rep[md5] = fid          # last wins, matching the upstream perl

um = open(O(".uniq.md5"), "w"); ua = open(O(".uniq.id_ann"), "w"); us = open(O(".uniq.seq"), "w")
fa = open(os.path.join(S, P + ".uniq.id_ann")); fs = open(os.path.join(S, P + ".uniq.seq"))
kept = repointed = 0
for lm in open(os.path.join(S, P + ".uniq.md5")):
    la, ls = fa.readline(), fs.readline()
    m = lm.strip()
    if m in md5s:
        fid, _, ann = la.rstrip("\n").partition("\t")
        if local_rep.get(m) and local_rep[m] != fid:
            fid = local_rep[m]; repointed += 1
        um.write(lm); ua.write("%s\t%s\n" % (fid, ann)); us.write(ls); kept += 1
for f in (um, ua, us): f.close()
if repointed:
    print("  re-pointed %d representative id(s) to a feature inside the subset" % repointed)
print("%s: %d genomes, %d features, %d unique sequences" % (a.name, len(keep), n, kept))
