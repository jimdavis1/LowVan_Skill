#!/usr/bin/env python3
"""Emit per-module metadata in the 7-column shape make_panel.py reads.

Contigs.metadata as BV-BRC ships it is three columns -- genome_id, name,
length -- and make_panel.py reads species from column 5 and length from
column 6. Fed the three-column file it sees species "?" for every genome and
length 0, so --per-species collapses the whole module into ONE bucket and the
panel is 8 genomes of whatever species happens to sort first. That is not a
stratified panel; it is a sample of the commonest species.

SPECIES KEY: the leading field of the BV-BRC genome id is the NCBI taxon id,
which is assigned at the type/species level -- 276 distinct ones across
Mastadenovirus, 15 across Siadenovirus. The strain-level genome_name is not
usable as a key (every genome is its own "species"), and word-prefix
heuristics on it merge unrelated viruses.
"""
import os, sys

W = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
lens = {}
for line in open(os.path.join(W, "Contigs.metadata")):
    f = line.rstrip("\n").split("\t")
    if len(f) > 2 and f[2].isdigit():
        lens[f[0]] = f[2]

rows = {}
for line in open(os.path.join(W, "A.id_name_gs.inscope")):
    f = line.rstrip("\n").split("\t")
    if len(f) < 4 or not f[3]:
        continue
    gid, name, fam, gen = f[0], f[1], f[2], f[3]
    rows.setdefault(gen, []).append(
        (gid, name, fam, gen, gid.split(".")[0], lens.get(gid, "0"), ""))

for gen, rs in sorted(rows.items()):
    out = os.path.join(W, "panel", "%s.metadata" % gen)
    with open(out, "w") as fh:
        for r in rs:
            fh.write("\t".join(r) + "\n")
    print("%-16s %5d genomes, %3d species keys -> %s"
          % (gen, len(rs), len({r[4] for r in rs}), out))
