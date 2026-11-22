#!/usr/bin/env python3
"""Homology rescue for sequences no triage rule matched.

references/annotation-triage.md: blast every unmatched protein against the named
collections OF ITS OWN MODULE and adopt at >= 80% identity over >= 85% of the
query.

THE COVERAGE FLOOR IS 0.85, RAISED FROM 0.60 ON THE CURATOR'S INSTRUCTION,
1 October 2026: "bring the coverage up to 85% if we ever do the rescue again,
i don't want to slurp back in all those short proteins."

Read what this does and does not do before relying on it. `ln / qlen` is
coverage of the QUERY, so it rejects a partial ALIGNMENT -- a protein matching
a collection member over only part of itself. It does NOT by itself reject a
short PROTEIN: a truncated sequence aligns over nearly all of itself and so
scores a high query coverage by construction. Raising the floor to 0.85 still
tightens adoption considerably, because most borderline hits are partial
alignments; it is not a truncation filter.

What actually bounds truncation is the feature's length window, applied after
adoption. Measured over the six Alsuviricetes modules rescued at 0.60: of 2,196
sequences adopted into collections, 20 (0.9%) fell below their collection's
previous MINIMUM length, worst Bromoviridae REP1 at 727 aa against a prior
floor of 908. Those 20 are what invalidated 13 features' min_len/max_len. If
short proteins must be excluded outright rather than merely bounded, the
instrument is a floor at the collection's existing minimum, not this flag. Below that, leave it out -- not adopted, and not swept into an
uncharacterized bag either, because that asserts a homology the number does not
support.

The threshold here MUST equal the one filter_unchar.py uses, or a protein 85%
identical to a named feature is adopted into that feature *and* left in the
grab-bag, and both profiles fire on one locus.

Adoption routes to the FEATURE KEY, never straight to the fasta, so the
collection's length window still applies afterwards and a rescued fragment
lands in <FEAT>.outliers.fasta like any other fragment.
"""
import os
import sys, re, sys, subprocess, collections, argparse, tempfile

ap = argparse.ArgumentParser()
ap.add_argument("--workdir", default=".")
ap.add_argument("--module",  required=True)
ap.add_argument("--dump",    required=True)
ap.add_argument("--min-id",  type=float, default=80.0)
ap.add_argument("--min-qcov",type=float, default=0.85)
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
    p = l.rstrip("\n").split("\t")
    if len(p) >= 4: g2g[p[0]] = p[3] or "unclassified"
#  BV-BRC writes uniq.seq and uniq.id_ann with a TRAILING TAB, so every line
#  splits into three fields, not two. dict() over that raises
#  "dictionary update sequence element #0 has length 3" and this script died
#  on every Alsuviricetes dump -- which is why no taxon in that class had ever
#  been rescued. Take the first two fields and ignore the rest.
seq = {}
for l in open(os.path.join(W, a.dump + ".uniq.seq")):
    f = l.rstrip("\n").split("\t")
    if len(f) >= 2 and f[1]:
        seq[f[0]] = f[1]

md5 = [l.strip() for l in open(os.path.join(W, a.dump + ".uniq.md5"))]
ann = [l.rstrip("\n").split("\t") for l in open(os.path.join(W, a.dump + ".uniq.id_ann"))]
if not (len(ann) == len(md5)):
    sys.exit("uniq.md5 and uniq.id_ann differ in length -- run check_dump.py")

#  rebuild the unassigned set: anything in the dump that is in no collection
#  file. build_collections.py writes the FEATURE ID as the fasta header
#  (>fig|28347.83.peg.1), not the md5, so the md5 cannot be recovered from the
#  header at all -- the old rsplit("|") took "28347.83.peg.1" and compared it
#  to md5s, which never matched, leaving `assigned` permanently empty. Had the
#  parse bugs above been fixed on their own, this would have re-rescued every
#  sequence that was ALREADY binned and duplicated it into the collections.
#  Match on the feature id instead, which is what the header actually holds.
#  THIS SCRIPT COULD NOT SEE ITS OWN OUTPUT, so it was NOT IDEMPOTENT.
#  Adoption below writes the header as ">fid|genus|md5", but this loop read
#  back split()[0] and compared it to a BARE feature id. A rescued record
#  therefore never appeared in `assigned`, and a second run re-adopted every
#  sequence the first run had already placed.
#
#  Measured on Arteriviridae, 19 Nov 2026: a second rescue over collections
#  that already held the first one's output duplicated 245 records -- 91 into
#  GP5, 37 into ORF1A, 25 into GP2, 20 into GP3, 17 into GP4. Nothing failed.
#  The collections simply grew, and the duplicated members would have been
#  double-weighted in every alignment and profile built from them.
#
#  UNASSIGNED_TRACKING.tsv over-reporting is a SEPARATE defect with the same
#  symptom, and it is fixed at the bottom of this file, not here. See the
#  comment on _rewrite_tracking(): this script never rewrote that file, so the
#  install gate counted sequences the rescue had already binned. 143 were
#  listed as uncalled when 109 of them were in collections -- which is what
#  sent the curator looking for a missing-sequence defect that did not exist.
#
#  Take the first two "|" fields of a fig id, which is the feature id itself,
#  and fall back to the whole token for any other id shape.
def _fid(tok):
    p = tok.split("|")
    return "|".join(p[:2]) if len(p) >= 2 and p[0] == "fig" else tok

assigned = set()
for fn in os.listdir(cdir):
    if not fn.endswith(".fasta"): continue
    for line in open(os.path.join(cdir, fn)):
        if line.startswith(">"): assigned.add(_fid(line[1:].split()[0].strip()))

#  PROJECT RULE: a record the source annotates as partial or incomplete is
#  DISCARDED, and that has to hold HERE as well. This script re-derives the
#  unassigned set from the dump, so without this test it hands every discarded
#  fragment straight back: on Rotavirus it adopted 36,800 of them at >=80%
#  identity, because they ARE genuine VP4/VP7/VP6 by homology -- they are
#  simply not whole. Shared with the binner so the two cannot disagree; they
#  did, and Sapovirus `Partial Capsid Protein` was binned out and rescued back
#  in the same run.
from annotation_rules import PARTIAL  # noqa: E402

#  PROJECT RULE, the same argument one step further: a string the taxon script
#  sent to NOT_MODELLED is a CURATOR'S DECISION that the module does not call
#  that product, and this script rebuilds its pool from "everything in the dump
#  that is in no collection", which includes every one of them. Without this
#  test the rescue hands them straight back, and it does so at high identity
#  precisely because they ARE homologous to a named feature -- which is why
#  they were excluded rather than ignored.
#
#  Measured on Arteriviridae, 18 Nov 2026: pp1ab, `polyprotein 1aTF` and
#  `polyprotein 1aN` are the -1 and -2 ribosomal frameshift products of the
#  replicase and are declared `special: transcript_edit` with no PSSM. The
#  rescue adopted 31 of them into ORF1A, because a 1,221 aa 1aTF product shares
#  its entire N-terminus with pp1a and so clears both the identity and the
#  query-coverage floors. A frameshift fusion inside the ORF1a collection is
#  exactly the contamination step 2a exists to prevent, and nothing downstream
#  would have shown it.
NOTMOD = set()
_nm = os.path.join(W, "collections", "NOT_MODELLED.tsv")
if os.path.exists(_nm):
    for _l in open(_nm):
        _p = _l.rstrip("\n").split("\t")
        if len(_p) >= 2 and _p[0] != "features":
            NOTMOD.add(" ".join(_p[1].split()).lower())

todo, n_partial, n_notmod = [], 0, 0
for i, m in enumerate(md5):
    fid = ann[i][0] if ann[i] else ""
    an = ann[i][1] if len(ann[i]) > 1 else ""
    if not fid or fid in assigned or not seq.get(m): continue
    if PARTIAL.search(" ".join((an or "").split())):
        n_partial += 1; continue
    if " ".join((an or "").split()).lower() in NOTMOD:
        n_notmod += 1; continue
    gid = fid.split("|", 1)[1].rsplit(".", 2)[0]
    todo.append((fid, g2g.get(gid, "unclassified"), m, an, seq[m]))
sys.stderr.write("%d sequences matched no rule; %d discarded as partial-annotated; "
                 "%d held out as NOT_MODELLED by the taxon script\n"
                 % (len(todo), n_partial, n_notmod))
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

def _rewrite_tracking():
    """Recompute UNASSIGNED_TRACKING.tsv's `not_already_called` column.

    ONLY collection_engine.py ever wrote this file, and it writes it BEFORE
    this script runs. Adoption below moves sequences into the collections but
    left the file saying they were uncalled, so the install gate -- which gates
    on exactly this column -- reported a pool that no longer existed, for every
    module that was ever rescued. Nothing failed and no number went down, which
    is why it survived: an inflated uncalled count looks like unfinished
    curation, and the obvious response is to go looking for the missing
    sequences.

    Measured on Arteriviridae, 19 Nov 2026: the gate reported 143 uncalled
    across 4 strings when 109 of them were already binned -- 41 `glycoprotein`
    and 32 `envelope glycoprotein` in GP5, and so on. Only the 34
    `hypothetical protein` records were genuinely uncalled. Acting on the
    inflated number duplicated 245 records into the collections before the
    arithmetic was checked.

    `features` and `uniq_seqs` describe the annotation group and do not change
    when a sequence is adopted, so they are preserved byte-for-byte; only the
    third column is recomputed, against the collections as they now stand.
    """
    p = os.path.join(W, "collections", "UNASSIGNED_TRACKING.tsv")
    if not os.path.exists(p):
        return
    called = set()
    for fn in os.listdir(cdir):
        if not fn.endswith(".fasta"): continue
        cur = []
        for line in open(os.path.join(cdir, fn)):
            if line.startswith(">"):
                if cur: called.add("".join(cur)); cur = []
            else: cur.append(line.strip())
        if cur: called.add("".join(cur))
    #  annotation -> the unique sequences deposited under it
    by_ann = collections.defaultdict(set)
    for i, m in enumerate(md5):
        an = " ".join((ann[i][1] if len(ann[i]) > 1 else "").split())
        if m in seq: by_ann[an].add(seq[m])
    rows, changed = [], 0
    with open(p) as fh:
        hdr = fh.readline().rstrip("\n")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4: rows.append(line.rstrip("\n")); continue
            an = f[-1]
            new = len(by_ann.get(" ".join(an.split()), set()) - called)
            if str(new) != f[2]: changed += 1
            rows.append("\t".join([f[0], f[1], str(new), an]))
    with open(p, "w") as fh:
        fh.write(hdr + "\n")
        for r in rows: fh.write(r + "\n")
    print("UNASSIGNED_TRACKING.tsv: recomputed uncalled counts, %d row(s) changed"
          % changed)


if a.write:
    for tgt, recs in adopt.items():
        p = os.path.join(cdir, tgt + ".fasta")
        with open(p, "a") as f: f.writelines(recs)
        print("appended %d to %s" % (len(recs), os.path.basename(p)))
    with open(os.path.join(W, "collections", "RESCUED.tsv"), "w") as f:
        f.write("count\tannotation\tgenus\tadopted_as\n")
        for (an, gen, key), c in counts.most_common():
            f.write("%d\t%s\t%s\t%s\n" % (c, an, gen, key))
    _rewrite_tracking()
else:
    print("\n(dry run -- rerun with --write)")
