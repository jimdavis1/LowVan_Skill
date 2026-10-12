#!/usr/bin/env python3
"""Test EVERY sequence in every collection against each feature's canonical.

    python3 audit_contamination.py --workdir . --module Orbivirus
    python3 audit_contamination.py ... --write      # then REBUILD

`qc_cross_feature` compares one cluster master at a time, so a collection
cross-contaminated by a handful of mislabelled deposits is fixed one cluster
per rebuild -- Orbivirus went 10 flags, then 4, then more, because each
re-clustering exposes the next layer. This asks the question once of every
sequence: 49 contaminating VP2 records that would have taken roughly five more
rebuild rounds to surface came out in a single pass.

**Compare to the canonical, never to the collection you are auditing.** The
obvious test -- blast the suspect cluster against every collection and take
the best hit -- does not work, because the collections contaminate each other:
a true Orbivirus VP2 scored 1881 against the VP2 collection and 1872 against
VP3, a 9-bit margin that decides nothing. Against the canonicals the same
comparison returns 400-2000 bits.

The canonical for a feature is the median-length member of its LARGEST
cluster, which is the most densely-supported form of that protein. A sequence
is reassigned only when it beats its current feature's canonical by --margin
bits, so a genuinely divergent species -- whose VP2 may resemble no canonical
closely -- is left where the source put it rather than dragged to whichever
canonical happens to be nearest.
"""
import argparse, collections, os, subprocess, tempfile

ap = argparse.ArgumentParser()
ap.add_argument("--workdir", default=".")
ap.add_argument("--module", required=True)
ap.add_argument("--margin", type=float, default=50.0,
                help="bits the winner must beat the incumbent by (default 50)")
ap.add_argument("--write", action="store_true")
a = ap.parse_args()

A = os.path.join(a.workdir, "Alignments", a.module)
C = os.path.join(a.workdir, "collections", a.module)


def readfa(p):
    out, h, parts = [], None, []
    for l in open(p, errors="replace"):
        if l.startswith(">"):
            if h is not None: out.append((h.rstrip("\n"), "".join(parts)))
            h, parts = l, []
        else: parts.append(l.strip())
    if h is not None: out.append((h.rstrip("\n"), "".join(parts)))
    return out


keys = sorted(d for d in os.listdir(A) if os.path.isdir(os.path.join(A, d)))
tmp = tempfile.mkdtemp(prefix="contam.")
with open(tmp + "/ref.faa", "w") as r:
    for k in keys:
        d = os.path.join(A, k, "corrected_alis")
        if not os.path.isdir(d): d = os.path.join(A, k, "alis")
        best, bn = None, 0
        for fn in os.listdir(d):
            if fn.endswith(".fa"):
                recs = readfa(os.path.join(d, fn))
                if len(recs) > bn: bn, best = len(recs), recs
        if not best: continue
        med = sorted(best, key=lambda x: len(x[1].replace("-", "")))[len(best) // 2]
        r.write(">%s\n%s\n" % (k, med[1].replace("-", "")))

n_q = 0
with open(tmp + "/q.faa", "w") as q:
    for k in keys:
        p = os.path.join(C, k + ".fasta")
        if not os.path.exists(p): continue
        for i, (h, s) in enumerate(readfa(p)):
            q.write(">%s@@%d\n%s\n" % (k, i, s)); n_q += 1
subprocess.run(["makeblastdb", "-in", tmp + "/ref.faa", "-dbtype", "prot",
                "-out", tmp + "/ref"], check=True, capture_output=True)
r = subprocess.run(["blastp", "-query", tmp + "/q.faa", "-db", tmp + "/ref",
                    "-outfmt", "6 qseqid sseqid bitscore", "-num_threads", "6",
                    "-max_target_seqs", "20"], capture_output=True, text=True)
best = collections.defaultdict(dict)
for l in r.stdout.splitlines():
    qid, sid, bits = l.split("\t"); b = float(bits)
    if b > best[qid].get(sid, 0): best[qid][sid] = b

move = collections.Counter(); nohit = 0
plan = collections.defaultdict(dict)
for qid, d in best.items():
    k, i = qid.split("@@"); i = int(i)
    top = max(d.items(), key=lambda x: x[1])
    if top[0] == k: continue
    if top[1] - d.get(k, 0.0) < a.margin: continue
    move[(k, top[0])] += 1; plan[k][i] = top[0]
print("  %d sequences tested against %d canonicals" % (n_q, len(keys)))
if not move:
    print("  no cross-feature contamination above %.0f bits" % a.margin); raise SystemExit(0)
print("  %-6s -> %-6s %s" % ("from", "to", "sequences"))
for (f, t), n in move.most_common(): print("  %-6s -> %-6s %6d" % (f, t, n))
print("  %d sequences in total" % sum(move.values()))
if not a.write:
    print("  (report only -- rerun with --write)"); raise SystemExit(0)
add = collections.defaultdict(list)
for k, m in plan.items():
    p = os.path.join(C, k + ".fasta"); recs = readfa(p); keep = []
    for i, (h, s) in enumerate(recs):
        if i in m: add[m[i]].append((h, s))
        else: keep.append((h, s))
    with open(p, "w") as f:
        for h, s in keep: f.write("%s\n%s\n" % (h, s))
for t, recs in add.items():
    with open(os.path.join(C, t + ".fasta"), "a") as f:
        for h, s in recs: f.write("%s\n%s\n" % (h, s))
print("  written -- rebuild every affected feature before measuring anything")
