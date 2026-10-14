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

Every CLUSTER of every feature contributes a reference (its median-length
member), and a sequence is scored against the best cluster of each feature.
Using one reference per feature -- the largest cluster's median -- makes the
reference set the best-sequenced lineage rather than the feature, and a
divergent species then matches whichever lineage reference is nearest. A
sequence
is reassigned only when it beats its current feature's canonical by --margin
bits, so a genuinely divergent species -- whose VP2 may resemble no canonical
closely -- is left where the source put it rather than dragged to whichever
canonical happens to be nearest.
"""
import argparse, collections, json, os, subprocess, tempfile

ap = argparse.ArgumentParser()
ap.add_argument("--workdir", default=".")
ap.add_argument("--module", required=True)
ap.add_argument("--margin", type=float, default=50.0,
                help="bits the winner must beat the incumbent by (default 50)")
ap.add_argument("--write", action="store_true")
a = ap.parse_args()

A = os.path.join(a.workdir, "Alignments", a.module)
C = os.path.join(a.workdir, "collections", a.module)


def feature_types(workdir, module):
    """key -> feature_type from the module JSON, {} if it is not written yet."""
    p = os.path.join(workdir, "%s_Viral_PSSM.json" % module)
    if not os.path.exists(p):
        return {}
    try:
        d = json.load(open(p))
    except ValueError:
        return {}
    f = (d.get(module) or {}).get("features") or {}
    return {k: (v.get("feature_type") or "CDS") for k, v in f.items()}


#  A mature peptide IS a piece of its parent CDS and matches it at 100% by
#  construction, so "move VP8* into VP4" is not a finding -- it is the
#  definition of VP8*. Rotavirus reported 5 of them. qc_cross_feature learned
#  this already; the same rule belongs here, where the consequence is worse
#  because --write would delete the derived feature's collection into its
#  parent.
FTYPE = feature_types(a.workdir, a.module)


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

#  ONE reference per cluster, not one per feature.
#
#  Taking only the largest cluster's median makes the reference set the
#  best-sequenced LINEAGE, not the feature. On Seadornavirus every one of the
#  twelve canonicals came from Banna virus, because Banna is 598 of the 1,202
#  records and so wins every largest-cluster; the Kadipiro, Liao ning and
#  Mangshi proteins then matched whichever Banna canonical happened to be
#  closest. The audit reported 99 moves -- 14% of the module -- including a
#  perfectly symmetric VP5->VP6 24 / VP6->VP5 24, which is the signature of an
#  artefact rather than a mislabel: real contamination runs one way. Applying
#  them would have shuffled three species' proteins at random.
#
#  Every cluster contributes a reference, tagged with its feature, and a
#  sequence is scored against the BEST cluster of each feature. A divergent
#  species is then judged against its own lineage's cluster wherever one
#  exists, and only loses to another feature if it genuinely resembles that
#  feature more than anything its own feature contains.
nref = 0
with open(tmp + "/ref.faa", "w") as r:
    for k in keys:
        d = os.path.join(A, k, "corrected_alis")
        if not os.path.isdir(d): d = os.path.join(A, k, "alis")
        if not os.path.isdir(d): continue
        for fn in sorted(os.listdir(d)):
            if not fn.endswith(".fa"): continue
            recs = readfa(os.path.join(d, fn))
            if not recs: continue
            med = sorted(recs, key=lambda x: len(x[1].replace("-", "")))[len(recs) // 2]
            r.write(">%s@@%s\n%s\n" % (k, fn[:-3], med[1].replace("-", "")))
            nref += 1

#  Observed length range per feature, from its own collection. A move has to
#  be consistent with it: an Orbivirus VP4 of 586 residues was proposed for
#  VP2, whose members run 870-1300, on a 51-bit hit that was noise -- it had
#  scored ZERO against every VP4 cluster, so anything at all beat its
#  incumbent. Identity alone cannot catch that; length can.
span, qlen = {}, {}
n_q = 0
with open(tmp + "/q.faa", "w") as q:
    for k in keys:
        p = os.path.join(C, k + ".fasta")
        if not os.path.exists(p): continue
        recs = readfa(p)
        L = sorted(len(x[1]) for x in recs)
        if L:
            lo = L[int(0.02 * (len(L) - 1))]
            hi = L[int(0.98 * (len(L) - 1))]
            span[k] = (lo * 0.80, hi * 1.20)
        for i, (h, s) in enumerate(recs):
            q.write(">%s@@%d\n%s\n" % (k, i, s)); n_q += 1
        qlen[k] = [len(x[1]) for x in recs]
subprocess.run(["makeblastdb", "-in", tmp + "/ref.faa", "-dbtype", "prot",
                "-out", tmp + "/ref"], check=True, capture_output=True)
r = subprocess.run(["blastp", "-query", tmp + "/q.faa", "-db", tmp + "/ref",
                    "-outfmt", "6 qseqid sseqid bitscore", "-num_threads", "6",
                    "-max_target_seqs", "20"], capture_output=True, text=True)
#  collapse cluster-level hits to the best score per FEATURE
best = collections.defaultdict(dict)
for l in r.stdout.splitlines():
    qid, sid, bits = l.split("\t"); b = float(bits)
    feat = sid.split("@@")[0]
    if b > best[qid].get(feat, 0): best[qid][feat] = b

move = collections.Counter(); blocked = collections.Counter()
parent = collections.Counter(); nohit = 0
plan = collections.defaultdict(dict)
for qid, d in best.items():
    k, i = qid.split("@@"); i = int(i)
    top = max(d.items(), key=lambda x: x[1])
    if top[0] == k: continue
    if top[1] - d.get(k, 0.0) < a.margin: continue
    #  the destination has to be able to hold a sequence this long
    if FTYPE.get(k) == "mat_peptide" and FTYPE.get(top[0], "CDS") != "mat_peptide":
        parent[(k, top[0])] += 1; continue
    lo, hi = span.get(top[0], (0, 10 ** 9))
    if not (lo <= qlen[k][i] <= hi):
        blocked[(k, top[0])] += 1; continue
    move[(k, top[0])] += 1; plan[k][i] = top[0]
print("  %d sequences tested against %d cluster references across %d features"
      % (n_q, nref, len(keys)))
if not move:
    print("  no cross-feature contamination above %.0f bits" % a.margin); raise SystemExit(0)
print("  %-6s -> %-6s %s" % ("from", "to", "sequences"))
for (f, t), n in move.most_common(): print("  %-6s -> %-6s %6d" % (f, t, n))
print("  %d sequences in total" % sum(move.values()))
if parent:
    print("  %d mature-peptide sequence(s) matched their parent CDS, as they must: %s"
          % (sum(parent.values()),
             ", ".join("%s->%s %d" % (a_, b_, n) for (a_, b_), n in parent.most_common(6))))
if blocked:
    print("  %d further hit the margin but were BLOCKED on length: %s"
          % (sum(blocked.values()),
             ", ".join("%s->%s %d" % (a_, b_, n) for (a_, b_), n in blocked.most_common(6))))
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
