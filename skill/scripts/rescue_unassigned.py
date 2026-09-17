#!/usr/bin/env python3
"""Homology rescue for sequences no triage rule matched.

references/annotation-triage.md: blast every unmatched protein against the named
collections OF ITS OWN MODULE and adopt at >= 80% identity over >= 60% of the
query. Below that, leave it out -- not adopted, and not swept into an
uncharacterized bag either, because that asserts a homology the number does not
support.

The threshold here MUST equal the one filter_unchar.py uses, or a protein 85%
identical to a named feature is adopted into that feature *and* left in the
grab-bag, and both profiles fire on one locus.

Adoption routes to the FEATURE KEY, never straight to the fasta, so the
collection's length window still applies afterwards and a rescued fragment
lands in <FEAT>.outliers.fasta like any other fragment.
"""
import os, re, sys, subprocess, collections, argparse, tempfile

ap = argparse.ArgumentParser()
ap.add_argument("--workdir", default=".")
ap.add_argument("--module",  required=True)
ap.add_argument("--dump",    required=True)
ap.add_argument("--min-id",  type=float, default=80.0)
ap.add_argument("--min-qcov",type=float, default=0.60)
ap.add_argument("--expected", required=True,
                help="KEY:lo-hi,KEY:lo-hi ... the same windows build_collections used")
ap.add_argument("--write", action="store_true")
a = ap.parse_args()

W = a.workdir
EXPECTED = {}
for part in a.expected.split(","):
    k, r = part.split(":"); lo, hi = r.split("-")
    EXPECTED[k] = (int(lo), int(hi))

cdir = os.path.join(W, "collections", a.module)

# the sequences that matched no rule
g2g = {}
for l in open(os.path.join(W, a.dump + ".id_name_gs")):
    p = l.rstrip("\n").split("\t"); g2g[p[0]] = p[3] or "unclassified"
seq = dict(l.rstrip("\n").split("\t") for l in open(os.path.join(W, a.dump + ".uniq.seq")))

# rebuild the unassigned set: anything in the dump that is in no collection file
assigned = set()
for fn in os.listdir(cdir):
    if not fn.endswith(".fasta"): continue
    for line in open(os.path.join(cdir, fn)):
        # NB the feature id itself contains a pipe (fig|1213422.5.CDS.2), so the
        # md5 is the LAST field, not field 2.
        if line.startswith(">"): assigned.add(line.rstrip("\n").rsplit("|", 1)[-1])

md5 = [l.strip() for l in open(os.path.join(W, a.dump + ".uniq.md5"))]
ann = [l.rstrip("\n").split("\t") for l in open(os.path.join(W, a.dump + ".uniq.id_ann"))]
todo = []
for m, (fid, an) in zip(md5, ann):
    if m in assigned or not seq.get(m): continue
    gid = fid.split("|", 1)[1].rsplit(".", 2)[0]
    todo.append((fid, g2g.get(gid, "unclassified"), m, an, seq[m]))
sys.stderr.write("%d sequences matched no rule\n" % len(todo))
if not todo: sys.exit(0)

tmp = tempfile.mkdtemp()
qf = os.path.join(tmp, "q.faa")
with open(qf, "w") as f:
    for fid, gen, m, an, s in todo:
        f.write(">%s|%s|%s\n%s\n" % (fid, gen, m, s))

# subject = the in-range named collections only
sf = os.path.join(tmp, "s.faa")
with open(sf, "w") as f:
    for key in EXPECTED:
        p = os.path.join(cdir, key + ".fasta")
        if not os.path.exists(p): continue
        for line in open(p):
            f.write(line.replace(">", ">%s@" % key, 1) if line.startswith(">") else line)

subprocess.run(["makeblastdb", "-in", sf, "-dbtype", "prot", "-out", os.path.join(tmp, "db")],
               check=True, stdout=subprocess.DEVNULL)
out = subprocess.run(["blastp", "-query", qf, "-db", os.path.join(tmp, "db"),
                      "-outfmt", "6 qseqid sseqid pident length qlen bitscore",
                      "-evalue", "1e-5", "-max_target_seqs", "5", "-num_threads", "8"],
                     check=True, capture_output=True, text=True).stdout

best = {}
for line in out.splitlines():
    q, s, pid, ln, qlen, bits = line.split("\t")
    pid, ln, qlen, bits = float(pid), int(ln), int(qlen), float(bits)
    if pid < a.min_id or ln / qlen < a.min_qcov: continue
    if q not in best or bits > best[q][2]:
        best[q] = (s.split("@")[0], pid, bits, ln / qlen)

info = {("%s|%s|%s" % (f, g, m)): (f, g, m, an, s) for f, g, m, an, s in todo}
assert len(info) == len(todo), "query id collision"
adopt = collections.defaultdict(list); rej = []
counts = collections.Counter()
for q, (key, pid, bits, qcov) in best.items():
    fid, gen, m, an, s = info[q]
    lo, hi = EXPECTED[key]
    tgt = key if lo <= len(s) <= hi else key + ".outliers"
    adopt[tgt].append(">%s|%s|%s\n%s\n" % (fid, gen, m, s))
    counts[(an, gen, key)] += 1
for q in info:
    if q not in best: rej.append(info[q])

print("adopted %d of %d at >=%.0f%% identity over >=%.0f%% of the query"
      % (sum(len(v) for v in adopt.values()), len(todo), a.min_id, 100 * a.min_qcov))
print("\n%-6s %-34s %-16s %s" % ("n", "annotation", "genus", "adopted as"))
for (an, gen, key), c in counts.most_common(30):
    print("%-6d %-34s %-16s %s" % (c, an[:34], gen, key))
print("\nleft out (below threshold, deliberately not swept into a grab-bag): %d" % len(rej))
rc = collections.Counter((x[3], x[1]) for x in rej)
for (an, gen), c in rc.most_common(12):
    print("   %-4d %-34s %s" % (c, an[:34], gen))

if a.write:
    for tgt, recs in adopt.items():
        p = os.path.join(cdir, tgt + ".fasta")
        with open(p, "a") as f: f.writelines(recs)
        print("appended %d to %s" % (len(recs), os.path.basename(p)))
    with open(os.path.join(W, "collections", "RESCUED.tsv"), "w") as f:
        f.write("count\tannotation\tgenus\tadopted_as\n")
        for (an, gen, key), c in counts.most_common():
            f.write("%d\t%s\t%s\t%s\n" % (c, an, gen, key))
else:
    print("\n(dry run -- rerun with --write)")
