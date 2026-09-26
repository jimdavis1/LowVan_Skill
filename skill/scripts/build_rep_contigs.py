#!/usr/bin/env python3
"""Emit the module's rep contigs, greedily largest-cluster-first.

`repcontig_budget.py` answers whether a taxon *can* be routed inside the
budget; it is a go/no-go measurement and writes nothing. This writes the
references that decision implies, from contigs already on disk, so the set
that ships is the same set that was measured.

The selection rule is the one the budget script models: cluster the genomes at
--mi, sort clusters by size, take the exemplar of each until --budget is
spent. Largest-first is what makes a small reference set cover a large
fraction of the taxon; picking one per species instead is how a 141-species
family ends up needing 141 references.

    python3 build_rep_contigs.py --contigs Contigs --meta Beta_MP.meta.tsv \
            --module Trivirinae --budget 25 --out Rep-Contigs

Writes <Module>.<n>.dna (one contig each, header ">ACCESSION Genome name")
and close_genomes.json, which is what the module JSON's close_genomes block
is built from.
"""
import os, sys, json, glob, argparse, tempfile, subprocess, collections

ap = argparse.ArgumentParser()
ap.add_argument("--contigs", required=True, help="directory of per-genome FASTA")
ap.add_argument("--meta", required=True, help="tsv: genome_id, name, family, genus, ...")
ap.add_argument("--module", required=True)
ap.add_argument("--out", default="Rep-Contigs")
ap.add_argument("--budget", type=int, default=25)
ap.add_argument("--mi", type=float, default=0.7)
ap.add_argument("--min-len", type=int, default=0)
ap.add_argument("--threads", type=int, default=6)
a = ap.parse_args()

name = {}
for l in open(a.meta):
    p = l.rstrip("\n").split("\t")
    if p and p[0]:
        name[p[0]] = p[1] if len(p) > 1 else p[0]

#  genome id -> (accession, description, sequence).  The accession is on the
#  contig header BV-BRC served, not in the metadata table, which is why this
#  reads the contigs rather than re-querying.
gseq, gacc = {}, {}
for f in sorted(glob.glob(os.path.join(a.contigs, "*"))):
    g = os.path.basename(f)
    for ext in (".fna", ".fa", ".fasta"):
        if g.endswith(ext):
            g = g[: -len(ext)]
            break
    acc, parts = None, []
    for line in open(f):
        if line.startswith(">"):
            if acc is None:
                acc = line[1:].split()[0]
        else:
            parts.append(line.strip())
    s = "".join(parts)
    if len(s) < a.min_len:
        continue
    #  --meta SELECTS, it does not merely name. It used to be read only into
    #  the name lookup while the candidate set came entirely from the --contigs
    #  glob, so passing a subset metadata changed nothing -- three different
    #  per-genus subsets of Bromoviridae returned the identical 474/683 and the
    #  same 25 accessions, which is how the gap was found. A split test is
    #  worthless if the restriction is ignored.
    if name and g not in name:
        continue
    gseq[g], gacc[g] = s, acc or g
if not gseq:
    sys.exit("no contigs in %s" % a.contigs)

tmp = tempfile.mkdtemp(prefix="repc.")
fa = os.path.join(tmp, "g.fna")
with open(fa, "w") as fh:
    for g, s in gseq.items():
        fh.write(">%s\n%s\n" % (g, s))
pre = os.path.join(tmp, "c")
r = subprocess.run(["mmseqs", "easy-cluster", fa, pre, os.path.join(tmp, "t"),
                    "--min-seq-id", str(a.mi), "-c", "0.5", "--cov-mode", "1",
                    "--threads", str(a.threads)], capture_output=True, text=True)
if not os.path.exists(pre + "_cluster.tsv"):
    sys.exit("mmseqs produced no clustering:\n" + r.stderr[-2000:])

clu = collections.defaultdict(list)
for l in open(pre + "_cluster.tsv"):
    rep, mem = l.rstrip("\n").split("\t")
    clu[rep].append(mem)

order = sorted(clu.items(), key=lambda kv: (-len(kv[1]), kv[0]))
chosen = order[: a.budget]
N = len(gseq)
covered = sum(len(v) for _, v in chosen)

os.makedirs(a.out, exist_ok=True)
for old in glob.glob(os.path.join(a.out, "%s.*.dna" % a.module)):
    os.remove(old)

close = {}
print("%-28s %d genomes, %d clusters at %.0f%% identity"
      % (a.module, N, len(clu), 100 * a.mi))
for i, (rep, mem) in enumerate(chosen, 1):
    fn = "%s.%d.dna" % (a.module, i)
    with open(os.path.join(a.out, fn), "w") as fh:
        fh.write(">%s %s\n" % (gacc[rep], name.get(rep, rep)))
        s = gseq[rep]
        for j in range(0, len(s), 60):
            fh.write(s[j:j + 60] + "\n")
    close[fn] = {"genome_ids": gacc[rep], "genome_name": name.get(rep, rep)}
    print("   %-28s %-12s %5d nt  covers %4d genome(s)  %s"
          % (fn, gacc[rep], len(gseq[rep]), len(mem), name.get(rep, rep)[:34]))

json.dump(close, open(os.path.join(a.out, "close_genomes.json"), "w"),
          indent=1, sort_keys=True)
print("\n   %d reference(s) cover %d/%d genomes = %.1f%%"
      % (len(chosen), covered, N, 100.0 * covered / N))
if len(clu) > a.budget:
    print("   %d cluster(s) left uncovered, holding %d genome(s) -- name this in "
          "the coverage audit" % (len(clu) - a.budget, N - covered))
