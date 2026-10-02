#!/usr/bin/env python3
"""Draw a held-out evaluation panel from the contigs the build already has.

Building the rep contigs (step 6c) downloads every contig in the taxon from
BV-BRC. The held-out panel is therefore free: take genomes from that directory,
stratified across species so one dominant species cannot carry the score, and
**exclude every genome that became a reference**.

That exclusion is the whole point and it is easy to skip. A panel containing
the genomes the references were chosen from measures routing against itself --
self-recall one level up. The same mistake in a different costume cost this
project a real measurement once already: the closterovirid transcript-edit
reference sets were extracted from the very panels the quality run scores, so
the measured RdRp miss rate was a floor rather than an estimate.

    python3 make_panel.py --contigs Contigs --metadata Contigs.metadata \\
            --close-genomes Rep-Contigs/close_genomes.json \\
            --out Panel --tsv panel.tsv --per-species 8 --max 120

Writes <out>/<genome_id>.fna plus the headerless `file, name, taxon_id` TSV
that make_gto.py wants, then:

    python3 make_gto.py     --fasta-dir Panel --metadata panel.tsv --out gto
    python3 run_gto_eval.py --gto-dir gto --repo $LOWVAN_DATA_DIR --out gto_out

Use evaluate_module.py only when comparing against a SUBMITTER's GenBank
annotation is the actual question; it goes to NCBI eutils one accession at a
time and is blocked on the network throughout.
"""
import argparse, collections, json, os, random, shutil, sys

ap = argparse.ArgumentParser()
ap.add_argument("--contigs", required=True, help="directory of <genome_id>.fna")
ap.add_argument("--metadata", required=True,
                help="Contigs.metadata from fetch_contigs.py: genome_id, name, "
                     "family, genus, species, length, accession")
ap.add_argument("--close-genomes", default=None,
                help="Rep-Contigs/close_genomes.json -- the genomes to EXCLUDE")
ap.add_argument("--out", required=True)
ap.add_argument("--tsv", required=True)
ap.add_argument("--per-species", type=int, default=8)
ap.add_argument("--max", type=int, default=120)
ap.add_argument("--min-len", type=int, default=0)
ap.add_argument("--seed", type=int, default=11)
a = ap.parse_args()

meta = {}
for line in open(a.metadata):
    f = line.rstrip("\n").split("\t")
    if not f or f[0] in ("genome_id", ""): continue
    meta[f[0]] = {"name": f[1] if len(f) > 1 else "",
                  "species": f[4] if len(f) > 4 else "?",
                  "length": int(f[5]) if len(f) > 5 and f[5].isdigit() else 0,
                  "accession": f[6] if len(f) > 6 else ""}
if not meta: sys.exit("no rows read from %s" % a.metadata)

refacc, refgid = set(), set()
if a.close_genomes:
    for v in json.load(open(a.close_genomes)).values():
        g = v.get("genome_ids")
        if g: refacc.add(g); refgid.add(g)

pool = collections.defaultdict(list)
skipped_ref = 0
for gid, m in meta.items():
    if not os.path.exists(os.path.join(a.contigs, gid + ".fna")): continue
    if m["accession"] in refacc or gid in refgid:
        skipped_ref += 1; continue
    if a.min_len and m["length"] < a.min_len: continue
    pool[m["species"]].append(gid)

if not pool: sys.exit("every genome was excluded -- check --close-genomes")

random.seed(a.seed)
picked = []
for sp in sorted(pool, key=lambda s: -len(pool[s])):
    picked += random.sample(pool[sp], min(a.per_species, len(pool[sp])))
picked = picked[:a.max]

os.makedirs(a.out, exist_ok=True)
rows = []
for gid in picked:
    shutil.copy(os.path.join(a.contigs, gid + ".fna"),
                os.path.join(a.out, gid + ".fna"))
    rows.append("%s.fna\t%s\t%s" % (gid, meta[gid]["name"], gid.split(".")[0]))
open(a.tsv, "w").write("\n".join(rows) + "\n")

print("panel: %d genome(s) across %d species, drawn from %d available"
      % (len(picked), len({meta[g]["species"] for g in picked}), sum(map(len, pool.values()))))
if a.close_genomes:
    print("       %d reference genome(s) excluded -- the panel is held out" % skipped_ref)
else:
    print("       NO --close-genomes given: the panel may contain the genomes the\n"
          "       references were chosen from, which is not a held-out measurement")
print("       wrote %s/ and %s" % (a.out, a.tsv))
