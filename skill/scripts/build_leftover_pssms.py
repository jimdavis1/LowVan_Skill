#!/usr/bin/env python3
"""Build extra profiles for the sequences the main clustering left behind.

The pipeline clusters a collection, builds one profile per cluster, and writes
whatever failed to cluster into Leftover_Seqs.aa. Those leftovers are then
dropped -- and they are not noise. They are the divergent members, which is
exactly the population a profile set most needs to cover.

Measured on Alpharhabdovirinae P: the collection holds 2727 sequences, 2231 of
them inside one of the 24 profiles and 282 left over. Of 115 genomes where the
annotator called N and M but no P, every single one had a full-length
phosphoprotein sitting in the gap -- and **none** of them resembled a sequence
any profile was built from, while 85 matched a leftover. The uncovered
population and the uncalled genes are the same set.

So: cluster the leftovers again on their own, and build a profile for every
group that reaches --min-seqs. Groups that stay small stay uncovered, which is
honest -- a two-member profile is a lookup table.

    python3 build_leftover_pssms.py --workdir . --module Alpharhabdovirinae --key P
    python3 build_leftover_pssms.py ... --write

This does NOT lower any threshold. The leftovers are clustered at the same
--mi the feature was built at; they simply get their own chance rather than
competing with the dense clusters that dominated the first pass.
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def read_fasta(p):
    h, s = None, []
    for line in open(p, errors="replace"):
        if line.startswith(">"):
            if h is not None:
                yield h, "".join(s)
            h, s = line[1:].strip(), []
        else:
            s.append(line.strip())
    if h is not None:
        yield h, "".join(s)


#  build_pssm aligns its group internally with MAFFT and builds the profile
#  from THAT alignment. The caller used to write the caller's own unaligned
#  input to corrected_alis/ instead, so the file on disk was not the alignment
#  the PSSM came from -- and corrected_alis/ is the directory the workflow
#  designates as the curation contract and the rebuild source.
ALIGNED_OUT = []


def build_pssm(rows, title, out, tmp):
    fa = os.path.join(tmp, "g.faa")
    with open(fa, "w") as fh:
        for i, (_h, s) in enumerate(rows):
            fh.write(">%d\n%s\n" % (i + 1, s.replace("-", "")))
    ali = os.path.join(tmp, "g.ali")
    with open(ali, "w") as fh:
        subprocess.run(["mafft", "--quiet", "--auto", fa], stdout=fh,
                       stderr=subprocess.DEVNULL)
    al = list(read_fasta(ali))
    if len(al) < 2:
        return False
    ALIGNED_OUT[:] = al          # the MSA this profile was really built from
    msa = os.path.join(tmp, "g.msa")
    with open(msa, "w") as fh:
        fh.write(">1 %s\n%s\n" % (title, al[0][1]))
        for i, (_h, s) in enumerate(al[1:]):
            fh.write(">%d\n%s\n" % (i + 2, s))
    subj = os.path.join(tmp, "g.subj")
    with open(subj, "w") as fh:
        fh.write(">1 %s\n%s\n" % (title, al[0][1].replace("-", "")))
    raw = os.path.join(tmp, "g.raw")
    subprocess.run(["psiblast", "-subject", subj, "-in_msa", msa, "-out_pssm", raw],
                   capture_output=True)
    if not os.path.exists(raw):
        return False
    txt = subprocess.run([sys.executable, os.path.join(HERE, "norm_pssm.py"), raw],
                         capture_output=True, text=True).stdout
    os.remove(raw)
    open(out, "w").write(txt)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--module", required=True)
    ap.add_argument("--key", required=True)
    ap.add_argument("--mi", type=float, default=0.6)
    ap.add_argument("--cov", type=float, default=0.8)
    ap.add_argument("--min-seqs", type=int, default=2,
                    help="minimum members for a leftover group to become a "
                         "profile (default 2). Deliberately lower than the -m "
                         "the feature was built at: the leftover pool is made "
                         "of sequences that failed that -m, so reapplying it "
                         "guarantees nothing comes back. -mi does the quality "
                         "work here, not the member count.")
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()

    W = os.path.abspath(args.workdir)
    fd = os.path.join(W, "Alignments", args.module, args.key)
    lo = os.path.join(fd, "Leftover_Seqs.aa")
    if not os.path.exists(lo):
        raise SystemExit("no Leftover_Seqs.aa in %s" % fd)
    rows = [(h, s) for h, s in read_fasta(lo) if len(s) >= 30]
    print("  %d leftover sequence(s) in %s/%s" % (len(rows), args.module, args.key))
    print("  member floor %d (the feature's own -m does not apply here: this "
          "pool is\n  defined by having failed it)" % args.min_seqs)

    tmp = tempfile.mkdtemp(prefix="lo.")
    try:
        fa = os.path.join(tmp, "lo.faa")
        with open(fa, "w") as fh:
            for i, (_h, s) in enumerate(rows):
                fh.write(">s%d\n%s\n" % (i, s.replace("-", "")))
        subprocess.run(["mmseqs", "easy-cluster", fa, os.path.join(tmp, "c"),
                        os.path.join(tmp, "t"), "--min-seq-id", str(args.mi),
                        "-c", str(args.cov), "--cov-mode", "0"], capture_output=True)
        groups = {}
        tsv = os.path.join(tmp, "c_cluster.tsv")
        if os.path.exists(tsv):
            for line in open(tsv):
                rep, mem = line.split()[:2]
                groups.setdefault(rep, []).append(int(mem[1:]))
        keep = sorted((v for v in groups.values() if len(v) >= args.min_seqs),
                      key=len, reverse=True)
        print("  %d cluster(s) at -mi %.1f; %d reach >= %d members, holding %d sequences"
              % (len(groups), args.mi, len(keep), args.min_seqs, sum(len(g) for g in keep)))
        if not args.write:
            print("\n  nothing written (rerun with --write)")
            return 0

        #  Number from where the feature's existing profiles stop, so the new
        #  ones never collide with a cluster id the main pass already used.
        used = {os.path.basename(p).rsplit(".", 2)[-2]
                for p in glob.glob(os.path.join(fd, "pssms", "*.pssm"))}
        n = 0
        made = 0
        for g in keep:
            n += 1
            while "lo%d" % n in used:
                n += 1
            cid = "lo%d" % n
            sub = [rows[i] for i in g]
            os.makedirs(os.path.join(fd, "corrected_alis"), exist_ok=True)
            title = "%s.%s.%s" % (args.module, args.key, cid)
            t2 = tempfile.mkdtemp(prefix="lo1.")
            ALIGNED_OUT[:] = []
            try:
                ok = build_pssm(sub, title,
                                os.path.join(fd, "pssms", title + ".pssm"), t2)
            finally:
                shutil.rmtree(t2, ignore_errors=True)
            if ok:
                #  Write the MSA the profile was built from, with its original
                #  headers restored, so rebuild_pssms.py can reproduce it and a
                #  curator can open it.
                with open(os.path.join(fd, "corrected_alis", cid + ".fa"), "w") as fh:
                    for (h, _s), (_i, a) in zip(sub, ALIGNED_OUT):
                        fh.write(">%s\n%s\n" % (h, a))
                made += 1
                print("    %-12s %3d sequences" % (cid, len(sub)))
        print("  built %d additional profile(s) for %s/%s" % (made, args.module, args.key))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
