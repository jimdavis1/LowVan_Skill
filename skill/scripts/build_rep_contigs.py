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
import os, re, sys, json, glob, argparse, tempfile, subprocess, collections

ap = argparse.ArgumentParser()
ap.add_argument("--contigs", required=True, help="directory of per-genome FASTA")
ap.add_argument("--meta", required=True, help="tsv: genome_id, name, family, genus, ...")
ap.add_argument("--module", required=True)
ap.add_argument("--out", default="Rep-Contigs")
ap.add_argument("--budget", type=int, default=25)
ap.add_argument("--mi", type=float, default=0.7)
ap.add_argument("--min-len", type=int, default=0)
ap.add_argument("--threads", type=int, default=6)
ap.add_argument("--segment-level", action="store_true",
                help="SEGMENTED, REASSORTING taxa: cluster individual segments and choose "
                     "genomes by how many uncovered segment clusters each one adds, instead "
                     "of clustering whole concatenated genomes. See the note below.")
a = ap.parse_args()

name = {}
for l in open(a.meta):
    p = l.rstrip("\n").split("\t")
    if p and p[0]:
        name[p[0]] = p[1] if len(p) > 1 else p[0]

#  genome id -> (accession, description, sequence).  The accession is on the
#  contig header BV-BRC served, not in the metadata table, which is why this
#  reads the contigs rather than re-querying.
gseq, gacc, grecs = {}, {}, {}
corrupt = []
for f in sorted(glob.glob(os.path.join(a.contigs, "*"))):
    g = os.path.basename(f)
    for ext in (".fna", ".fa", ".fasta"):
        if g.endswith(ext):
            g = g[: -len(ext)]
            break
    #  Keep the contigs SEPARATE. This used to join every record in the file
    #  into one string and write one contig, which for a segmented virus
    #  fuses the whole genome under the first segment's header: a Rotavirus
    #  reference came out as 18 kb of "28875.58294_seg1". Routing still
    #  worked -- BLASTn finds the right HSP in a concatenated subject -- so
    #  nothing failed, but the shipped file said something untrue about the
    #  genome, and it disagreed with the influenza rep contigs, which carry
    #  one record per segment. Clustering below still uses the concatenation,
    #  because that is the genome-level comparison it wants.
    acc, recs, hdr, parts = None, [], None, []
    for line in open(f):
        if line.startswith(">"):
            if hdr is not None:
                recs.append((hdr, "".join(parts)))
            hdr, parts = line[1:].rstrip("\n"), []
            if acc is None:
                acc = hdr.split()[0]
        else:
            parts.append(line.strip())
    if hdr is not None:
        recs.append((hdr, "".join(parts)))
    s = "".join(x[1] for x in recs)
    if len(s) < a.min_len:
        continue
    #  A reference whose sequence is missing one of the four bases ENTIRELY is
    #  not a sequence. BV-BRC carries such records -- 48 of 710 Coltivirus
    #  segment records hold no cytosine at all, every one of them described
    #  "Modified Microbial Nucleic Acid." -- and they are poison here
    #  specifically, because they cluster happily with EACH OTHER. The greedy
    #  cover then spends references on them and reports the genomes they
    #  "cover" as routed: Coltivirus measured 25 references and 93.4% with
    #  them in, and 23 references and 100.0% with them out. Ten of those 25
    #  references were corrupt records covering other corrupt records.
    #  The build repo shipped eight such references across six modules before
    #  this check existed.
    comp = set(s.upper()) & set("ACGT")
    if len(comp) < 4:
        corrupt.append((g, len(s), "".join(sorted(set("ACGT") - comp))))
        continue
    #  --meta SELECTS, it does not merely name. It used to be read only into
    #  the name lookup while the candidate set came entirely from the --contigs
    #  glob, so passing a subset metadata changed nothing -- three different
    #  per-genus subsets of Bromoviridae returned the identical 474/683 and the
    #  same 25 accessions, which is how the gap was found. A split test is
    #  worthless if the restriction is ignored.
    if name and g not in name:
        continue
    #  `acc` is the first contig's header token, which for an unsegmented
    #  genome is a real accession (NC_002023). For a genome this project
    #  merged from segments it is our own "<gid>_seg1" label, and shipping
    #  that as the reference's identity says the whole genome is segment 1.
    gseq[g], gacc[g], grecs[g] = s, re.sub(r"_seg\d+$", "", acc or g), recs
if corrupt:
    print("  %d record(s) EXCLUDED: sequence is missing a base entirely "
          "(not a sequence; they cluster with each other and inflate coverage)"
          % len(corrupt))
    for g, ln, miss in corrupt[:8]:
        print("     %-22s %7d nt  no %s" % (g, ln, miss))
    if len(corrupt) > 8:
        print("     ... and %d more" % (len(corrupt) - 8))
if not gseq:
    sys.exit("no contigs in %s" % a.contigs)

#  --------------------------------------------------------------------------
#  WHAT GETS CLUSTERED, and why it is not always the genome.
#
#  The default concatenates each genome and clusters that. For a REASSORTING
#  segmented virus that is the wrong model, and it fails in a way that looks
#  like a hopeless taxon rather than a bad measurement: a concatenated
#  rotavirus is a mosaic, so every fresh combination of otherwise-common
#  segments becomes its own cluster. Rotavirus gave 1,574 clusters for 6,648
#  genomes, 74% of them singletons, and a projected reach of 49%. Orbivirus
#  gives 54.7% at 30 references by this rule.
#
#  What actually happens at runtime is different. Every incoming record of a
#  segmented virus is ONE segment; a rep contig is a whole genome carrying all
#  of them; so an incoming segment routes if it resembles ANY segment of ANY
#  chosen reference. Measured that way the same Orbivirus set reaches 100% at
#  30, and Rotavirus 100% at 33.
#
#  --segment-level switches to that model: cluster segments, then choose
#  genomes greedily by how many uncovered SEGMENT clusters each one adds. Use
#  it for any taxon whose records arrive one segment at a time.
#  --------------------------------------------------------------------------
tmp = tempfile.mkdtemp(prefix="repc.")
fa = os.path.join(tmp, "g.fna")
with open(fa, "w") as fh:
    if a.segment_level:
        for g, recs in grecs.items():
            for i, (_h, cs) in enumerate(recs, 1):
                if len(cs) >= 300:
                    fh.write(">%s@@%d\n%s\n" % (g, i, cs))
    else:
        for g, s in gseq.items():
            fh.write(">%s\n%s\n" % (g, s))
pre = os.path.join(tmp, "c")
r = subprocess.run(["mmseqs", "easy-cluster", fa, pre, os.path.join(tmp, "t"),
                    "--min-seq-id", str(a.mi), "-c", "0.5", "--cov-mode", "1",
                    "--threads", str(a.threads)], capture_output=True, text=True)
if not os.path.exists(pre + "_cluster.tsv"):
    sys.exit("mmseqs produced no clustering:\n" + r.stderr[-2000:])

if a.segment_level:
    seg2clu = {}
    for l in open(pre + "_cluster.tsv"):
        rep, mem = l.rstrip("\n").split("\t")[:2]
        seg2clu[mem] = rep
    gclu = collections.defaultdict(set)
    for seg, rep in seg2clu.items():
        gclu[seg.split("@@")[0]].add(rep)
    covered, chosen_g, pool = set(), [], dict(gclu)
    while len(chosen_g) < a.budget:
        best = max(pool, key=lambda k: len(pool[k] - covered), default=None)
        if best is None or not (pool[best] - covered):
            break
        covered |= pool[best]; chosen_g.append(best); pool.pop(best)
    #  how many genomes each reference routes, for the printout: a genome is
    #  routed by the first chosen reference sharing a segment cluster with it
    chosen = []
    claimed = set()
    for g in chosen_g:
        mem = [h for h, v in gclu.items() if h not in claimed and (v & gclu[g])]
        claimed |= set(mem)
        chosen.append((g, mem))
    routed = len(claimed)
    print("%-28s %d genomes, %d SEGMENT clusters at %.0f%% identity (--segment-level)"
          % (a.module, len(gclu), len(set(seg2clu.values())), 100 * a.mi))
else:
    clu = collections.defaultdict(list)
    for l in open(pre + "_cluster.tsv"):
        rep, mem = l.rstrip("\n").split("\t")
        clu[rep].append(mem)
    order = sorted(clu.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    chosen = order[: a.budget]
    routed = sum(len(v) for _, v in chosen)
    print("%-28s %d genomes, %d clusters at %.0f%% identity"
          % (a.module, len(gseq), len(clu), 100 * a.mi))

N = len(gseq)
covered = routed

os.makedirs(a.out, exist_ok=True)
for old in glob.glob(os.path.join(a.out, "%s.*.dna" % a.module)):
    os.remove(old)

close = {}
for i, (rep, mem) in enumerate(chosen, 1):
    fn = "%s.%d.dna" % (a.module, i)
    with open(os.path.join(a.out, fn), "w") as fh:
        #  one record per source contig; a single-contig genome is unchanged
        for ci, (chdr, cseq) in enumerate(grecs[rep], 1):
            cacc = chdr.split()[0]
            fh.write(">%s %s\n" % (cacc, name.get(rep, rep)))
            for j in range(0, len(cseq), 60):
                fh.write(cseq[j:j + 60] + "\n")
    close[fn] = {"genome_ids": gacc[rep], "genome_name": name.get(rep, rep)}
    print("   %-28s %-12s %5d nt  covers %4d genome(s)  %s"
          % (fn, gacc[rep], len(gseq[rep]), len(mem), name.get(rep, rep)[:34]))

json.dump(close, open(os.path.join(a.out, "close_genomes.json"), "w"),
          indent=1, sort_keys=True)
print("\n   %d reference(s) cover %d/%d genomes = %.1f%%"
      % (len(chosen), covered, N, 100.0 * covered / N))
if a.segment_level:
    if N - covered:
        print("   %d genome(s) share no segment cluster with any reference -- name "
              "this in the coverage audit" % (N - covered))
elif len(clu) > a.budget:
    print("   %d cluster(s) left uncovered, holding %d genome(s) -- name this in "
          "the coverage audit" % (len(clu) - a.budget, N - covered))
