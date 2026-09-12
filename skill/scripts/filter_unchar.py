#!/usr/bin/env python3
"""Keep the UNCHAR grab-bag disjoint from the named features.

UNCHAR collects proteins no binning rule named. When one of them is also
near-identical to a member of a named collection in the same module, both
profiles hit the same locus and the annotator emits two features on identical
coordinates. That is always a defect, and it is not hypothetical: declaring
UNCHAR for Novirhabdovirus and Dichorhavirus immediately produced

    DUPLICATE P and Novir_unchar both at 1481-2149
    DUPLICATE N and dicho_unchar both at 228-1580

Exact-sequence de-duplication does not catch it, because the overlapping
proteins come from different genomes and differ by a few residues. This uses
blastp.

The grab-bag is the side that yields: a protein already covered by a named
feature does not belong in it.

The 80% default is the same line build_collections.py uses to decide that an
unnamed protein *is* a named feature. The two must agree. If this were looser,
a protein 85% identical to a named feature would be adopted into that feature
AND left in the grab-bag, and both profiles would fire on the same locus --
the duplicate-call defect this script exists to prevent.

    python3 filter_unchar.py --workdir .            # report
    python3 filter_unchar.py --workdir . --write    # rewrite the UNCHAR fastas

Run after build_collections.py and before building PSSMs. Requires blastp.
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile


def read_fasta(path):
    hdr, seq = None, []
    for line in open(path, errors="replace"):
        if line.startswith(">"):
            if hdr is not None:
                yield hdr, "".join(seq)
            hdr, seq = line.rstrip("\n"), []
        else:
            seq.append(line.strip())
    if hdr is not None:
        yield hdr, "".join(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--pident", type=float, default=80.0,
                    help="drop an UNCHAR member at or above this %% identity "
                         "to a named feature (default 80)")
    ap.add_argument("--cov", type=float, default=0.60,
                    help="minimum aligned fraction of the UNCHAR member (default 0.60)")
    ap.add_argument("--write", action="store_true")
    args = ap.parse_args()

    workdir = os.path.abspath(args.workdir)
    total = 0
    #  Matches UNCHAR.fasta and the UNC1..UNCn homology groups that replace it.
    #  Every one of them is a grab-bag member for this purpose: a protein
    #  already covered by a NAMED feature does not belong in any of them.
    for uf in sorted(glob.glob(os.path.join(workdir, "collections", "*", "UNC*.fasta"))):
        module = os.path.basename(os.path.dirname(uf))
        named = [p for p in glob.glob(os.path.join(workdir, "collections", module, "*.fasta"))
                 if not os.path.basename(p).startswith("UNC")]
        members = list(read_fasta(uf))
        if not members or not named:
            continue
        tmp = tempfile.mkdtemp(prefix="unchar.")
        try:
            db = os.path.join(tmp, "named.faa")
            with open(db, "w") as fh:
                for p in named:
                    feat = os.path.basename(p)[:-6].replace(".outliers", "")
                    for h, s in read_fasta(p):
                        fh.write(">%s@%s\n%s\n" % (feat, h[1:].split()[0], s))
            subprocess.run(["makeblastdb", "-in", db, "-dbtype", "prot",
                            "-out", os.path.join(tmp, "db")], capture_output=True)
            q = os.path.join(tmp, "q.faa")
            with open(q, "w") as fh:
                for i, (h, s) in enumerate(members):
                    fh.write(">u%d\n%s\n" % (i, s))
            r = subprocess.run(
                ["blastp", "-query", q, "-db", os.path.join(tmp, "db"), "-outfmt",
                 "6 qseqid sseqid pident length", "-evalue", "1e-5",
                 "-max_target_seqs", "5", "-num_threads", "4"],
                capture_output=True, text=True)
            drop = {}
            for line in r.stdout.splitlines():
                qid, sid, pid, alen = line.split("\t")
                i = int(qid[1:])
                if float(pid) >= args.pident and int(alen) >= args.cov * len(members[i][1]):
                    drop.setdefault(i, sid.split("@")[0])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

        if not drop:
            print("  %-24s %4d members, none overlap a named feature" % (module, len(members)))
            continue
        total += len(drop)
        by_feat = {}
        for i, f in drop.items():
            by_feat[f] = by_feat.get(f, 0) + 1
        print("  %-24s %4d members, %d overlap: %s"
              % (module, len(members), len(drop),
                 ", ".join("%s x%d" % (k, v) for k, v in sorted(by_feat.items()))))
        if args.write:
            keep = [(h, s) for i, (h, s) in enumerate(members) if i not in drop]
            with open(uf, "w") as fh:
                for h, s in keep:
                    fh.write("%s\n" % h)
                    for j in range(0, len(s), 60):
                        fh.write(s[j:j + 60] + "\n")
            print("  %-24s rewritten with %d members" % ("", len(keep)))

    print("\n  %d UNCHAR member(s) %s"
          % (total, "removed" if args.write else "would be removed (rerun with --write)"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
