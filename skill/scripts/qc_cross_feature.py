#!/usr/bin/env python3
"""Find clusters that were binned under the wrong feature.

BV-BRC annotation strings are the only thing that decides which collection a
protein lands in, so a mislabelled source feature puts a protein under the wrong
key, and the resulting PSSM then out-scores the correct one at the *other*
protein's locus. The symptom at runtime is two features called on identical
coordinates.

This blasts every cluster's master sequence against a database of all the
module's other collections and flags any cluster whose best hit is a different
feature. Run it before building the JSON.

    python3 qc_cross_feature.py --workdir . --module Alpharhabdovirinae
    python3 qc_cross_feature.py --workdir .            # every module

Requires blastp and makeblastdb on PATH.
"""

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict


def read_fasta(path):
    hdr, seq = None, []
    for line in open(path, errors="replace"):
        if line.startswith(">"):
            if hdr is not None:
                yield hdr, "".join(seq)
            hdr, seq = line[1:].rstrip("\n"), []
        else:
            seq.append(line.strip())
    if hdr is not None:
        yield hdr, "".join(seq)


def masters(workdir, module):
    """One representative per cluster: the alignment master, ungapped."""
    out = []
    base = os.path.join(workdir, "Alignments", module)
    if not os.path.isdir(base):
        return out
    for feat in sorted(os.listdir(base)):
        for sub in ("corrected_alis", "reclustered_alis", "truncated_alis"):
            d = os.path.join(base, feat, sub)
            if not os.path.isdir(d):
                continue
            for fa in sorted(glob.glob(os.path.join(d, "*.fa"))):
                num = os.path.basename(fa)[:-3]
                for hdr, seq in read_fasta(fa):
                    ident = re.sub(r"/\d+-\d+", "", hdr.split()[0])
                    out.append((feat, num, ident, seq.replace("-", "").replace(".", "").upper()))
                    break
    return out


def run_module(workdir, module, pident, cov, tmp):
    reps = masters(workdir, module)
    if not reps:
        print("  %s: no alignments found" % module)
        return []

    # database of every collection, tagged by feature
    dbf = os.path.join(tmp, "db.faa")
    collected = set()
    n = 0
    with open(dbf, "w") as out:
        for fa in sorted(glob.glob(os.path.join(workdir, "collections", module, "*.fasta"))):
            feat = os.path.basename(fa)[:-6]
            # "<FEAT>.outliers.fasta" is a length-split subset of <FEAT>, not a
            # separate feature -- fold it back or every feature flags itself.
            if feat.endswith(".outliers"):
                feat = feat[:-len(".outliers")]
            collected.add(feat)
            for hdr, seq in read_fasta(fa):
                out.write(">%s@%s\n%s\n" % (feat, hdr.split()[0], seq))
                n += 1
    if not n:
        print("  %s: no collections/ found -- skipping" % module)
        return []
    subprocess.run(["makeblastdb", "-in", dbf, "-dbtype", "prot", "-out",
                    os.path.join(tmp, "db")], capture_output=True)

    qf = os.path.join(tmp, "q.faa")
    with open(qf, "w") as out:
        for feat, num, ident, seq in reps:
            out.write(">%s|%s|%s\n%s\n" % (feat, num, ident, seq))

    res = subprocess.run(
        ["blastp", "-query", qf, "-db", os.path.join(tmp, "db"), "-outfmt",
         "6 qseqid sseqid pident length bitscore", "-evalue", "1e-5",
         "-max_target_seqs", "500", "-num_threads", "4"],
        capture_output=True, text=True)

    # query -> feature -> (best pident, aligned length, subject)
    best = defaultdict(dict)
    for line in res.stdout.splitlines():
        q, s, pid, alen, bits = line.split("\t")
        f, sid = s.split("@", 1)
        pid, alen = float(pid), int(alen)
        cur = best[q].get(f)
        if cur is None or pid > cur[0]:
            best[q][f] = (pid, alen, sid)

    flags = []
    for feat, num, ident, seq in reps:
        # Derived features (G_MAT, G_SP...) have no collection of their own and
        # are *supposed* to match their parent -- nothing to check.
        if feat not in collected:
            continue
        q = "%s|%s|%s" % (feat, num, ident)
        hits = best.get(q, {})
        cross = [(p, a, s, f) for f, (p, a, s) in hits.items()
                 if f != feat and f in collected]
        if not cross:
            continue
        pid, alen, sid, of = max(cross)
        # Near-identity across two feature labels means one of the two source
        # annotations is wrong. Which one is a judgement call, so report the
        # collision rather than asserting a direction. Collisions against the
        # uncharacterised grab-bag are a weaker signal -- the same protein
        # simply got binned both specifically and generically -- so they are
        # separated out rather than mixed in with true mislabels.
        if pid >= pident and alen >= cov * len(seq):
            kind = "UNCHAR" if "UNCHAR" in (feat, of) else "MISLABEL"
            flags.append((kind, module, feat, num, ident, len(seq), of, pid, alen, sid))
    return flags


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--module", help="one module (default: all under Alignments/)")
    ap.add_argument("--pident", type=float, default=90.0,
                    help="flag at or above this %% identity to another feature (default 90)")
    ap.add_argument("--cov", type=float, default=0.60,
                    help="minimum aligned fraction of the master (default 0.60)")
    args = ap.parse_args()

    workdir = os.path.abspath(args.workdir)
    adir = os.path.join(workdir, "Alignments")
    mods = [args.module] if args.module else sorted(
        d for d in os.listdir(adir) if os.path.isdir(os.path.join(adir, d)))

    tmp = tempfile.mkdtemp(prefix="xfeat.")
    allflags = []
    try:
        for m in mods:
            allflags += run_module(workdir, m, args.pident, args.cov, tmp)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    if not allflags:
        print("\nNo cross-feature contamination found.")
        return 0
    hard = [f for f in allflags if f[0] == "MISLABEL"]
    soft = [f for f in allflags if f[0] == "UNCHAR"]

    def show(rows, title):
        if not rows:
            return
        print("\n%s (%d)\n" % (title, len(rows)))
        print("  %-20s %-10s %-5s %-24s %5s %-10s %7s  %s"
              % ("module", "binned as", "clus", "master", "len", "collides", "%id", "with"))
        for _k, m, feat, num, ident, ln, of, pid, alen, sid in sorted(rows, key=lambda x: -x[7]):
            print("  %-20s %-10s %-5s %-24s %5d %-10s %6.1f%%  %s"
                  % (m, feat, num, ident[:24], ln, of, pid, sid[:26]))

    show(hard, "MISLABEL -- two named features share near-identical sequence")
    if hard:
        print("\n  One of the two source annotations is wrong. Whichever cluster is")
        print("  mislabelled will call its feature on the other protein's locus.")
        print("  Decide from length and gene order, then drop or move that cluster.")
    show(soft, "UNCHAR overlap -- same protein binned both specifically and generically")
    if soft:
        print("\n  Lower severity: fold these into the specific feature, or accept the")
        print("  overlap if the UNCHAR profile is a genuine grab-bag.")
    return 1 if hard else 0


if __name__ == "__main__":
    sys.exit(main())
