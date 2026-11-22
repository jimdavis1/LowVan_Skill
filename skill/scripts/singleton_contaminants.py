#!/usr/bin/env python3
"""Find collection members that belong to ANOTHER feature, however few there are.

`audit_bins.py` has a foreign-subset test, and it is the right test, but it
only reports a group once it reaches `SUBSET_MIN = 5` -- the point at which
foreign members are numerous enough to survive clustering and become a profile
of their own.

**A single foreign member can become a profile too.**  `build_leftover_pssms.py
--min-seqs 1` duplicates the row and builds a profile from one sequence, which
is correct and deliberate: a lineage-specific protein with one deposit must not
be dropped silently.  The cost is that a single MISLABELLED deposit gets the
same treatment, and the resulting profile carries no conservation information
at all -- every column is that one sequence -- so it matches its true family
far better than it matches the feature it was filed under.

Measured on Arteriviridae, 21 Nov 2026.  `fig|28344.22115.peg.9` was deposited
as `Envelope protein`, is 51 aa, and sits in the E collection.  E (ORF2b) is
~70 aa; ORF5a is ~51.  It is 98% identical over its whole length to ORF5a
deposits and matches nothing in E but itself.  It became `Arteriviridae.E.19`,
and that profile called E on the ORF5a reading frame in three genomes -- E and
ORF5a on byte-identical coordinates, a kilobase from where E belongs.  One
sequence in 588.

So: for every member of every collection, find its best hit in every OTHER
collection, and report any member that matches a foreign feature at high
identity across most of its own length.  No count floor.

Usage: singleton_contaminants.py <collections dir> [id] [cov] [A:B,C:D]
"""
import collections
import os
import shutil
import subprocess
import sys
import tempfile

CDIR = sys.argv[1]
MIN_ID = float(sys.argv[2]) if len(sys.argv) > 2 else 90.0
MIN_COV = float(sys.argv[3]) if len(sys.argv) > 3 else 0.85

#  PAIRS THE SUBJECT-COVERAGE TEST CANNOT SEPARATE, as A:B, comma separated.
#
#  Subject coverage distinguishes a contaminant (the SAME protein as something
#  else, so the hit covers both) from a product (PART of something else, so the
#  hit covers the query and a fraction of the subject).  That works whenever
#  the product is clearly shorter than its precursor, which is the usual case:
#  arterivirus nsp7alpha covers 6% of ORF1A.
#
#  It fails when a product is nearly as long as its precursor.  Measured on
#  Birnaviridae, 21 Nov 2026: pVP2 is VP2 plus the C-terminal structural
#  peptides -- 512 aa against 441 in IBDV -- so VP2 covers 86% of pVP2 and
#  clears an 85% floor.  All 398 VP2 members were reported, and the right
#  answer is that VP2 IS the matured pVP2, which is what the module models.
#  Lowering the floor would hide real contaminants; raising it would not fix
#  this pair.  So the pair is declared, and the declaration is short because
#  only near-length maturation steps need it.
PROTECT = set()
for _p in (sys.argv[4] if len(sys.argv) > 4 else "").split(","):
    if ":" in _p:
        _a, _b = _p.split(":", 1)
        PROTECT.add((_a.strip(), _b.strip()))
        PROTECT.add((_b.strip(), _a.strip()))

BV = "/Applications/BV-BRC.app/runtime/bin"
if os.path.isdir(BV):
    os.environ["PATH"] = BV + ":" + os.environ.get("PATH", "")


def read_fasta(p):
    out, h, s = {}, None, []
    for line in open(p):
        if line.startswith(">"):
            if h:
                out[h] = "".join(s)
            h, s = line[1:].split()[0].strip(), []
        else:
            s.append(line.strip())
    if h:
        out[h] = "".join(s)
    return out


feats = {}
for fn in sorted(os.listdir(CDIR)):
    if not fn.endswith(".fasta") or ".outliers" in fn:
        continue
    feats[fn[:-6]] = read_fasta(os.path.join(CDIR, fn))

owner, paths = {}, {}
tmp = tempfile.mkdtemp(prefix="singleton_")
for f, seqs in feats.items():
    for h in seqs:
        owner["%s__%s" % (f, h)] = f

print("%d collection(s), %d sequence(s); reporting a member whose best FOREIGN "
      "hit is >=%.0f%% identity over >=%.0f%% of the query\n"
      % (len(feats), len(owner), MIN_ID, 100 * MIN_COV))

#  ONE SEARCH PER FEATURE, AGAINST EVERY OTHER FEATURE.
#
#  The first version of this built a single database of everything and blasted
#  everything against it.  That is quadratic in the WRONG term: GP5 has 2,375
#  members which all hit each other, so almost every alignment computed and
#  reported was a within-collection hit that the next line threw away.  It ran
#  31 minutes of CPU without finishing.  Capping -max_target_seqs instead would
#  have been worse than slow: for a 2,375-member collection the top hits are
#  all same-feature, so the one foreign hit being looked for falls off the end
#  and the sweep reports CLEAN.
#
#  Excluding the queried feature from its own subject database asks the
#  question directly, and -max_target_seqs 5 is then ample.
best = {}
for f, seqs in sorted(feats.items()):
    if not seqs:
        continue
    q = os.path.join(tmp, "q.faa")
    db = os.path.join(tmp, "s.faa")
    with open(q, "w") as fh:
        for h, s_ in seqs.items():
            if s_:
                fh.write(">%s__%s\n%s\n" % (f, h, s_))
    n_sub = 0
    with open(db, "w") as fh:
        for g, other in feats.items():
            if g == f:
                continue
            for h, s_ in other.items():
                if s_:
                    fh.write(">%s__%s\n%s\n" % (g, h, s_))
                    n_sub += 1
    if not n_sub:
        continue
    subprocess.run(["makeblastdb", "-in", db, "-dbtype", "prot",
                    "-out", os.path.join(tmp, "db")],
                   capture_output=True, text=True)
    r = subprocess.run(
        ["blastp", "-query", q, "-db", os.path.join(tmp, "db"),
         "-outfmt", "6 qseqid sseqid pident length qlen slen bitscore",
         "-evalue", "1e-3", "-max_target_seqs", "5",
         "-task", "blastp-fast", "-num_threads", "8"],
        capture_output=True, text=True)
    hits = 0
    for line in r.stdout.splitlines():
        qid, sid, pid, aln, qlen, slen, bits = line.split("\t")
        pid, aln, qlen, slen = float(pid), int(aln), int(qlen), int(slen)
        bits = float(bits)
        #  THE HIT MUST COVER BOTH SEQUENCES, not just the query.
        #
        #  Requiring query coverage alone reports every mature peptide in the
        #  module, because a mature peptide IS a 100%-identical substring of
        #  its own precursor -- nsp7alpha covers 149 of its own 149 residues
        #  against ORF1A and is not a contaminant, it is the definition.  On
        #  the first run that buried the one real finding under eleven
        #  tautologies.
        #
        #  Subject coverage separates them without a hand-maintained list of
        #  precursor/product pairs: nsp7alpha covers ~6% of the 2,400 aa
        #  ORF1A, while the mislabelled ORF5a in the E bin covers 100% of the
        #  51 aa ORF5a deposits it matches.  A contaminant is the SAME protein
        #  as something else; a product is PART of something else.
        if pid < MIN_ID or aln < MIN_COV * qlen or aln < MIN_COV * slen:
            continue
        if (owner[qid], owner[sid]) in PROTECT:
            continue
        hits += 1
        if qid not in best or bits > best[qid][3]:
            best[qid] = (owner[sid], pid, aln, bits, qlen, sid)
    print("  %-12s %5d member(s) vs %5d foreign -> %d member(s) over threshold"
          % (f, len(seqs), n_sub, len({k for k in best if owner[k] == f})))

shutil.rmtree(tmp, ignore_errors=True)
print()

by = collections.defaultdict(list)
for q, v in best.items():
    by[(owner[q], v[0])].append((q, v))

if not by:
    print("  no member of any collection matches a foreign feature at that "
          "threshold.")
for (src, tgt), rows in sorted(by.items(), key=lambda x: -len(x[1])):
    print("== %s -> %s   %d member(s) of %d" % (src, tgt, len(rows), len(feats[src])))
    for q, v in sorted(rows, key=lambda x: -x[1][1])[:8]:
        print("     %-34s %3d aa  %.1f%% id over %d  best: %s"
              % (q.split("__", 1)[1], v[4], v[1], v[2], v[5].split("__", 1)[1]))
    if len(rows) > 8:
        print("     ... %d more" % (len(rows) - 8))
    print()
