#!/usr/bin/env python3
"""Project mature-peptide boundaries onto polyprotein-only genomes.

Most picornavirus genomes in BV-BRC carry one "polyprotein" CDS and no
mat_peptides, so a mature-peptide collection is built from the minority of
genomes somebody annotated -- 2 to 5 sequences per feature in Avihepatovirus,
against 196 polyproteins. This recovers the rest.

References are genomes whose polyprotein contains each of its own collection
mat_peptides as an exact substring: their cut coordinates are a fact of the
source annotation, not an inference. Every other polyprotein is globally
aligned to its nearest reference (blastp bitscore), and each reference cut is
carried across the alignment. Nothing is predicted; a cut is only ever moved
along an alignment.

Accepted only when
  * the target is >= --min-id identical to its reference over the alignment
  * the extracted segment is inside the feature's JSON min_len..max_len
  * the segment is >= --min-seg-id identical to the reference segment
so a cut drifting across a divergent region is refused, not shipped.

--split KEY:NEWKEY:REGEX cuts KEY in two at the first REGEX match (the NEWKEY
product starts at the match). Used where the source annotation lumps two
products that a documented site separates, e.g. DHAV 2A2|2A3.

Writes <workdir>/Extracted/<G>.<KEY>.faa, which the collection merge picks up.
"""
import argparse, collections, json, os, re, subprocess, sys, tempfile
from Bio import Align
from Bio.Align import substitution_matrices

ap = argparse.ArgumentParser()
ap.add_argument("--workdir", required=True)
ap.add_argument("--genus", required=True)
ap.add_argument("--dump", default="G")
ap.add_argument("--min-id", type=float, default=0.50)
ap.add_argument("--min-seg-id", type=float, default=0.35)
ap.add_argument("--split", action="append", default=[])
ap.add_argument("--first-product", default="",
                help="KEY:NEXTKEY -- KEY is the polyprotein N-terminus up to the start of "
                     "NEXTKEY. For a leader, which few genomes annotate but every polyprotein "
                     "carries. Derived, not projected: the coordinates are the reference's own.")
ap.add_argument("--threads", type=int, default=4)
a = ap.parse_args()
W, G = a.workdir, a.genus
C = os.path.join(W, "collections", G)

def fasta(p):
    out, h = collections.OrderedDict(), None
    for l in open(p):
        l = l.rstrip("\n")
        if l.startswith(">"): h = l[1:].split()[0]; out[h] = []
        elif h: out[h].append(l.strip())
    return {k: "".join(v) for k, v in out.items()}

mod = json.load(open(os.path.join(W, "%s_Viral_PSSM.json" % G)))[G]["features"]
skip = {"POLY", "POOL", "PRECURSORS"}
coll = {}
for f in os.listdir(C):
    k = f[:-6]
    #  only keys the module declares -- a stale collection (e.g. L before the Ldr rename) is ignored
    if f.endswith(".fasta") and k not in skip and k in mod: coll[k] = set(fasta(os.path.join(C, f)).values())
poly = fasta(os.path.join(C, "POLY.fasta"))
splits = []
for s in a.split:
    k, nk, rx = s.split(":", 2); splits.append((k, nk, re.compile(rx)))

# genome -> sequences, from the dump
md5seq = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(os.path.join(W, a.dump + ".uniq.seq"))}
bygen = collections.defaultdict(set)
for l in open(os.path.join(W, a.dump + ".id_md5")):
    fid, m = l.rstrip("\n").split("\t")
    gid = fid.split("|")[1].rsplit(".", 2)[0]
    if m in md5seq: bygen[gid].add(md5seq[m])

polyset = set(poly.values())
refs = {}                     # poly seq -> {key: (start, end)}
for gid, seqs in bygen.items():
    for p in seqs & polyset:
        cuts = {}
        for k, cs in coll.items():
            for s in seqs & cs:
                i = p.find(s)
                if i >= 0 and p.find(s, i + 1) < 0: cuts[k] = (i, i + len(s))
        if a.first_product:
            k0, k1 = a.first_product.split(":")
            if k1 in cuts and cuts[k1][0] > 0: cuts[k0] = (0, cuts[k1][0])
        if len(cuts) >= 5 and (p not in refs or len(cuts) > len(refs[p])): refs[p] = cuts
print("  %d polyproteins, %d reference polyproteins with >=5 placed mat_peptides" % (len(polyset), len(refs)))
if not refs: sys.exit("no references")
print("  features placed on references: %s" % dict(collections.Counter(k for c in refs.values() for k in c)))

# nearest reference per target, by blastp bitscore
td = tempfile.mkdtemp()
ref_ids = {"r%d" % i: p for i, p in enumerate(refs)}
tgt_ids = {"t%d" % i: p for i, p in enumerate(sorted(polyset))}
open(td + "/r.faa", "w").write("".join(">%s\n%s\n" % kv for kv in ref_ids.items()))
open(td + "/t.faa", "w").write("".join(">%s\n%s\n" % kv for kv in tgt_ids.items()))
subprocess.run(["makeblastdb", "-in", td + "/r.faa", "-dbtype", "prot", "-out", td + "/r"], check=True, capture_output=True)
bl = subprocess.run(["blastp", "-query", td + "/t.faa", "-db", td + "/r", "-outfmt", "6 qseqid sseqid bitscore",
                     "-max_target_seqs", "2000", "-evalue", "1e-20", "-num_threads", str(a.threads)],
                    check=True, capture_output=True, text=True).stdout
hits = collections.defaultdict(dict)
for l in bl.splitlines():
    q, s, b = l.split("\t"); b = float(b)
    hits[q][s] = max(b, hits[q].get(s, 0))

al = Align.PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                           open_gap_score=-10, extend_gap_score=-0.5)

ALPHA = set(al.substitution_matrix.alphabet)
def clean(x): return "".join(ch if ch in ALPHA else "X" for ch in x.upper())

def mapping(r, t):
    aln = al.align(clean(r), clean(t))[0]
    m = [None] * (len(r) + 1)
    same = alnlen = 0
    for (rs, re_), (ts, te) in zip(*aln.aligned):
        for d in range(re_ - rs):
            m[rs + d] = ts + d
            same += r[rs + d] == t[ts + d]
        alnlen += re_ - rs
    # unaligned ref residues take the next aligned target position
    nxt = len(t)
    for i in range(len(r), -1, -1):
        if m[i] is None: m[i] = nxt
        else: nxt = m[i]
    return m, same / float(max(len(r), len(t)))

def pid(x, y):
    a2 = al.align(clean(x), clean(y))[0]
    same = sum(x[i] == y[j] for (rs, re_), (ts, te) in zip(*a2.aligned) for i, j in zip(range(rs, re_), range(ts, te)))
    return same / float(max(len(x), len(y)))

out = collections.defaultdict(dict)   # key -> {seq: header}
why = collections.Counter()
for tid, t in tgt_ids.items():
    if tid not in hits: why["no reference hit"] += 1; continue
    ranked = sorted(hits[tid], key=lambda x: -hits[tid][x])
    cache = {}
    for k in sorted({k for c in refs.values() for k in c}):
        #  nearest reference that actually places this feature
        rid = next((x for x in ranked if k in refs[ref_ids[x]]), None)
        if rid is None: why[k + " no reference places it"] += 1; continue
        r = ref_ids[rid]
        if rid not in cache:
            cache[rid] = (list(range(len(r) + 1)), 1.0) if t == r else mapping(r, t)
        m, ident = cache[rid]
        if ident < a.min_id: why[k + " polyprotein below %.0f%% id" % (100 * a.min_id)] += 1; continue
        s, e = refs[r][k]
        seg, rseg = t[m[s]:m[e]], r[s:e]
        lo, hi = mod[k]["min_len"] if k in mod else 0, mod[k]["max_len"] if k in mod else 10 ** 6
        if not seg or not (lo <= len(seg) <= hi) or "X" in seg or "*" in seg: why[k + " length/ambiguity"] += 1; continue
        if seg != rseg and pid(rseg, seg) < a.min_seg_id: why[k + " segment id"] += 1; continue
        out[k][seg] = "%s|%s|ref=%s|id=%.2f" % (G, k, rid, ident)

# splits: applied to the projected segments and to the existing collection alike
for k, nk, rx in splits:
    pool = dict(out.pop(k, {}))
    for s in coll.get(k, ()): pool.setdefault(s, "%s|%s|collection" % (G, k))
    for s, h in pool.items():
        mt = rx.search(s)
        if not mt or mt.start() == 0: why[k + " split motif absent"] += 1; continue
        out[k][s[:mt.start()]] = h + "|split"
        out[nk][s[mt.start():]] = h.replace("|%s|" % k, "|%s|" % nk) + "|split"

os.makedirs(os.path.join(W, "Extracted"), exist_ok=True)
print("\n  %-6s %8s %8s %8s" % ("key", "extract", "new", "in coll"))
for k in sorted(out):
    p = os.path.join(W, "Extracted", "%s.%s.faa" % (G, k))
    with open(p, "w") as f:
        for i, (s, h) in enumerate(sorted(out[k].items())):
            f.write(">ext%d|%s\n%s\n" % (i + 1, h, s))
    new = len(set(out[k]) - coll.get(k, set()))
    print("  %-6s %8d %8d %8d" % (k, len(out[k]), new, len(coll.get(k, ()))))
print("\n  refused: %s" % dict(why))
