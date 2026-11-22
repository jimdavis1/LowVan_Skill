#!/usr/bin/env python3
"""copy_num re-audit, per taxon.

`copy_num` makes a feature ESSENTIAL: viral_genome_quality.pl scores only
features that carry it, and flags any genome that does not call one.  That is
right for a feature every member of the module really has, and wrong for a
feature a lineage genuinely lacks -- there it flags a complete, correctly
annotated genome as defective, and the flag is a statement about the module's
reach rather than about the genome.

So the question is never "what fraction of genomes call it".  It is: **is there
a COMPLETE genome, ROUTED TO THIS MODULE, in which this feature is never
called?**

Two filters that the first version of this script did not have, and both
changed the answer:

  * a genome in the run may have routed to a DIFFERENT module.  Its features
    then carry that module's family assignments, so tallying on the bare
    feature name credits another module's calls to this one -- and, worse,
    counts its missing features against this one.  `2916858.3` routed to
    Simarterivirinae and looked like an Arteriviridae genome calling 5 of 22.

  * the SOURCE's own completeness flag.  A 13.1 kb deposit is not a complete
    arterivirus genome merely because it clears a nucleotide floor; BV-BRC
    marks it `Partial`, and a feature uncalled in a partial genome is exactly
    the case `copy_num` exists to flag.  Length alone cannot see this.

Usage:
  copynum_audit.py <Anno_GTOs> <module json> [module] [meta.tsv] [min_nt]
"""
import json
import os
import sys
from collections import defaultdict

gdir = sys.argv[1]
mjson = sys.argv[2]
MODULE = sys.argv[3] if len(sys.argv) > 3 else "Arteriviridae"
META = sys.argv[4] if len(sys.argv) > 4 else None
MIN_NT = int(sys.argv[5]) if len(sys.argv) > 5 else 12000

feats = json.load(open(mjson))[MODULE]["features"]
essential = sorted(k for k, v in feats.items() if "copy_num" in v)

#  genome id -> source completeness flag, column 12 of the meta dump
status = {}
if META and os.path.exists(META):
    for line in open(META):
        c = line.rstrip("\n").split("\t")
        if len(c) > 11:
            status[c[0]] = c[11]

tax = defaultdict(lambda: {"n": 0, "called": defaultdict(int), "ids": []})
skipped = defaultdict(int)


def calls(g):
    """{feature: True} for the features THIS module called, and the module(s)
    that actually annotated the genome."""
    mods, ks = set(), set()
    for f in g.get("features", []):
        for fa in f.get("family_assignments", []) or []:
            if fa and fa[0] == "LOWVAN" and "." in str(fa[1]):
                parts = str(fa[1]).split(".")
                mods.add(parts[0])
                if parts[0] == MODULE:
                    ks.add(parts[1])
    return ks, mods


for fn in sorted(os.listdir(gdir)):
    if not fn.endswith(".gto"):
        continue
    #  the annotator REWRITES the gto's internal id, so the filename is the
    #  only link back to the deposit the metadata describes.
    gid = fn[:-4] if fn.endswith(".gto") else fn
    if gid.endswith(".anno"):
        gid = gid[:-5]
    try:
        g = json.load(open(os.path.join(gdir, fn)))
    except ValueError:
        skipped["unreadable gto"] += 1
        continue
    ks, mods = calls(g)
    if mods and MODULE not in mods:
        skipped["routed to %s" % "/".join(sorted(mods))] += 1
        continue
    if not mods:
        skipped["no LOWVAN call at all"] += 1
        continue
    st = status.get(gid)
    if status and st != "Complete":
        skipped["source says %s" % (st or "unknown")] += 1
        continue
    if sum(len(c.get("dna", "")) for c in g.get("contigs", [])) < MIN_NT:
        skipped["under %d nt" % MIN_NT] += 1
        continue
    rec = tax[gid.split(".")[0]]
    rec["n"] += 1
    rec["ids"].append(gid)
    for k in ks:
        rec["called"][k] += 1

n = sum(r["n"] for r in tax.values())
print("%s: %d complete genome(s) routed to this module, in %d taxa" % (
    MODULE, n, len(tax)))
for why, c in sorted(skipped.items(), key=lambda x: -x[1]):
    print("   skipped %4d  %s" % (c, why))
print()

print("features called, per taxon (of %d essential / %d total):" % (
    len(essential), len(feats)))
for t, r in sorted(tax.items(), key=lambda x: len(x[1]["called"])):
    print("   %-10s %3d genome(s)  %2d feature(s)" % (t, r["n"], len(r["called"])))
print()

bad = False
for k in essential:
    zero = [(t, r["n"]) for t, r in tax.items() if r["called"].get(k, 0) == 0]
    part = [(t, r["n"], r["called"][k]) for t, r in tax.items()
            if 0 < r["called"][k] < r["n"]]
    if not zero and not part:
        continue
    bad = True
    print("== %s" % k)
    for t, c in sorted(zero, key=lambda x: -x[1]):
        print("   NEVER called  taxon %-10s %3d complete genome(s)" % (t, c))
    for t, c, g_ in sorted(part, key=lambda x: x[2] - x[1]):
        print("   partial       taxon %-10s %3d/%3d" % (t, g_, c))
    print()

if not bad:
    print("every copy_num feature is called in every complete genome "
          "routed to this module.")
