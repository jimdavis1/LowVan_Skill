#!/usr/bin/env python3
"""Find out why a feature is called in fewer genomes than its neighbours.

Step 11 gives a per-protein call rate. Where one core protein trails the others,
the gene is either absent from those genomes or invisible to the profiles, and
those demand opposite responses. Guessing wastes the effort: the prior
expectation for Alpharhabdovirinae P was genome quality, and it was wrong.

Method: take genomes that called the target's genomic NEIGHBOURS but not the
target, cut the interval between them, ORF-call it in six frames, and ask two
questions of the longest ORF -- is it this protein at all, and is it closer to a
sequence that HAS a profile or one that does not.

    python3 analyze_feature_gap.py --workdir . --module Alpharhabdovirinae \
            --key P --left N --right M --ann-dir coverage_eval/ann

Four outcomes, four conclusions:

    interval short or ambiguous      genome quality; nothing to fix
    no ORF                           the gene really is absent
    ORF closest to a CLUSTERED member    threshold or coverage cutoff is wrong
    ORF closest to a LEFTOVER            the profile set misses that variant

The last is the common case, and the fix is build_leftover_pssms.py rather than
touching any threshold.

ORF-calling also recovers proteins the source does not have: 21 of 115
Alpharhabdovirinae phosphoproteins appear in no BV-BRC protein record. Pass
--save-extracted to write those to Extracted/<Module>.<KEY>.faa, which
build_collections.py merges back into the collection.
"""
import argparse, collections, glob, os, re, shutil, subprocess, sys, tempfile

CODON = {}
_b = "TCAG"
_a = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_i = 0
for x in _b:
    for y in _b:
        for z in _b:
            CODON[x + y + z] = _a[_i]; _i += 1


def rc(s):
    return s.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def orfs(seq, lo, hi):
    out = []
    for s in (seq, rc(seq)):
        for fr in range(3):
            p = "".join(CODON.get(s[j:j + 3].upper(), "X") for j in range(fr, len(s) - 2, 3))
            for m in re.finditer(r"M[^*]*", p):
                if lo <= len(m.group()) <= hi:
                    out.append(m.group())
    return sorted(out, key=len, reverse=True)


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


def blast_best(q, subj, tmp, tag):
    db = os.path.join(tmp, tag)
    subprocess.run(["makeblastdb", "-in", subj, "-dbtype", "prot", "-out", db],
                   capture_output=True)
    r = subprocess.run(["blastp", "-query", q, "-db", db, "-outfmt",
                        "6 qseqid bitscore", "-evalue", "1e-5", "-num_threads", "4",
                        "-max_target_seqs", "1"], capture_output=True, text=True)
    best = {}
    for line in r.stdout.splitlines():
        k, b = line.split("\t")
        best[k] = max(best.get(k, 0.0), float(b))
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--module", required=True)
    ap.add_argument("--key", required=True)
    ap.add_argument("--left", required=True, help="feature key immediately 5' of the target")
    ap.add_argument("--right", required=True, help="feature key immediately 3' of the target")
    ap.add_argument("--ann-dir", required=True, help="directory of *.feature.tbl")
    ap.add_argument("--genome-dir", required=True, help="directory of the genome fastas")
    ap.add_argument("--min-aa", type=int, default=80)
    ap.add_argument("--max-aa", type=int, default=1200)
    ap.add_argument("--save-extracted", action="store_true")
    args = ap.parse_args()

    W = os.path.abspath(args.workdir)
    fd = os.path.join(W, "Alignments", args.module, args.key)

    #  genomes that called both flanks and not the target
    sel = []
    for t in sorted(glob.glob(os.path.join(args.ann_dir, "*.feature.tbl"))):
        if os.path.getsize(t) == 0:
            continue
        f = collections.defaultdict(list); mod = None
        for line in open(t, errors="replace"):
            c = line.rstrip("\n").split("\t")
            if len(c) < 18:
                continue
            mod = c[11]; f[c[6]].append((int(c[7]), int(c[8])))
        if mod != args.module or args.key in f:
            continue
        if args.left not in f or args.right not in f:
            continue
        lo = max(max(a, b) for a, b in f[args.left])
        hi = min(min(a, b) for a, b in f[args.right])
        if hi - lo < 150:
            continue
        sel.append((os.path.basename(t)[:-12], lo, hi))
    print("  %d genome(s) called %s and %s but not %s"
          % (len(sel), args.left, args.right, args.key))
    if not sel:
        return 0

    tmp = tempfile.mkdtemp(prefix="gap.")
    try:
        q = os.path.join(tmp, "orfs.faa")
        stats = []
        with open(q, "w") as fh:
            for g, lo, hi in sel:
                p = os.path.join(args.genome_dir, g + ".fna")
                if not os.path.exists(p):
                    continue
                seq = "".join(l.strip() for l in open(p) if not l.startswith(">"))
                seg = seq[min(lo, hi):max(lo, hi) - 1]
                amb = sum(1 for c in seg.upper() if c not in "ACGT")
                o = orfs(seg, args.min_aa, args.max_aa)
                stats.append((g, len(seg), amb, len(o[0]) if o else 0))
                if o:
                    fh.write(">%s\n%s\n" % (g, o[0]))
        L = sorted(x[1] for x in stats)
        O = sorted(x[3] for x in stats if x[3])
        print("    interval nt   : min %d median %d max %d" % (L[0], L[len(L)//2], L[-1]))
        print("    with ambiguity: %d" % sum(1 for x in stats if x[2] > 0))
        print("    contain an ORF: %d of %d" % (len(O), len(stats)))
        if O:
            print("    ORF aa        : min %d median %d max %d" % (O[0], O[len(O)//2], O[-1]))
        if not O:
            print("\n  VERDICT: the gene is absent. Nothing to build.")
            return 0

        coll = os.path.join(W, "collections", args.module, args.key + ".fasta")
        hit = blast_best(q, coll, tmp, "coll") if os.path.exists(coll) else {}
        print("    match the %s collection: %d" % (args.key, len(hit)))

        lof = os.path.join(fd, "Leftover_Seqs.aa")
        clf = os.path.join(tmp, "clustered.faa")
        with open(clf, "w") as fh:
            for a in glob.glob(os.path.join(fd, "corrected_alis", "*.fa")) + \
                     glob.glob(os.path.join(fd, "reclustered_alis", "*.fa")):
                for h, s in read_fasta(a):
                    fh.write(">%s\n%s\n" % (h, s.replace("-", "")))
        lo_b = blast_best(q, lof, tmp, "lo") if os.path.exists(lof) else {}
        cl_b = blast_best(q, clf, tmp, "cl") if os.path.getsize(clf) else {}
        qs = {h for h, _s in read_fasta(q)}
        a = sum(1 for k in qs if k in lo_b and (k not in cl_b or lo_b[k] > cl_b[k]))
        b = sum(1 for k in qs if k in cl_b and (k not in lo_b or cl_b[k] >= lo_b[k]))
        print("\n    closest to a LEFTOVER  : %d" % a)
        print("    closest to a CLUSTERED : %d" % b)
        print("    neither                : %d" % (len(qs) - a - b))
        print("\n  VERDICT: %s" % (
            "the profile set does not cover these variants -- run "
            "build_leftover_pssms.py" if a > b else
            "a profile exists but does not fire -- check bit_cutoff and "
            "coverage_cutoff" if b else
            "too divergent to place; leave them out"))

        if args.save_extracted:
            allp = set()
            for line in open(os.path.join(W, "Rhabdo.uniq.seq"), errors="replace") \
                    if os.path.exists(os.path.join(W, "Rhabdo.uniq.seq")) else []:
                f2 = line.rstrip("\n").split("\t")
                if len(f2) > 1:
                    allp.add(f2[1].strip())
            new = [(h, s) for h, s in read_fasta(q) if s not in allp]
            if new:
                d = os.path.join(W, "Extracted"); os.makedirs(d, exist_ok=True)
                out = os.path.join(d, "%s.%s.faa" % (args.module, args.key))
                with open(out, "w") as fh:
                    fh.write("; %d %s sequences ORF-called from the %s-%s interval of\n"
                             "; genomes where the module called neither, and absent from\n"
                             "; the source protein export entirely.\n"
                             % (len(new), args.key, args.left, args.right))
                    for h, s in new:
                        fh.write(">extracted|%s %s\n%s\n" % (h, args.key, s))
                print("\n  wrote %d recovered sequence(s) to %s" % (len(new), out))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
