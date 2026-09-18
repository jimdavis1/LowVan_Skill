#!/usr/bin/env python3
"""Find clusters whose start codon is wrong for the feature as a whole.

`fasta-cluster-pssm-2.pl` already evaluates N-termini, but only *within* a
cluster: it asks whether the members agree with each other. A cluster whose
members were all annotated from the same wrong upstream Met agrees with itself
perfectly and passes. Nothing compares a cluster's start against the feature's
other clusters, and that is the case that does real damage.

Betaflexiviridae_MP CP cluster 18: 11 sequences, all 249-250 aa, against a
modal CP of 193. The pipeline logged it `18.fa 250 25 10 2 OK`. Its profile
anchored 57 codons upstream of every other CP profile, landed inside the
movement protein ORF in a different frame, hit a stop, and was rejected for
coverage -- while scoring 451 bits and beating the clean profiles at 263 and
247 that would have called the protein correctly. The coat protein went
uncalled on 26 genomes, 17% of the genus, and every other check passed.

    python3 qc_nterm_cross_cluster.py --workdir . --module <M> --key CP
    python3 qc_nterm_cross_cluster.py ... --write      # trim to the right Met

The test: take each cluster's master, align it to the master of the reference
cluster (the largest one at the feature's modal length), and look at where the
reference starts inside it. If the reference's residue 1 lands well inside the
query, the query carries an N-terminal extension the rest of the feature does
not have. `--write` trims the alignment to that column, but ONLY when there is
a Met there -- without one, the right start is a judgement for the curator and
the script refuses.
"""
import os, re, sys, glob, argparse, subprocess, tempfile, statistics

def read_fasta(path):
    recs, h, s = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if h is not None: recs.append((h, "".join(s)))
            h, s = line.rstrip("\n"), []
        else: s.append(line.strip())
    if h is not None: recs.append((h, "".join(s)))
    return recs

def ungapped(s): return re.sub(r"[-.]", "", s)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--module", required=True)
    ap.add_argument("--key", required=True)
    ap.add_argument("--min-extension", type=int, default=20,
                    help="flag a cluster whose master starts this many residues "
                         "or more upstream of the reference cluster (default 20)")
    ap.add_argument("--write", action="store_true",
                    help="trim flagged alignments to the internal Met and leave "
                         "a note; refuses when there is no Met at the offset")
    a = ap.parse_args()

    fd = os.path.join(os.path.abspath(a.workdir), "Alignments", a.module, a.key)
    alis = []
    for sub in ("corrected_alis", "reclustered_alis"):
        alis += sorted(glob.glob(os.path.join(fd, sub, "*.fa")))
    if not alis:
        sys.exit("no alignments in %s" % fd)

    info = {}
    for f in alis:
        recs = read_fasta(f)
        if not recs: continue
        lens = [len(ungapped(s)) for _, s in recs]
        info[f] = {"n": len(recs), "med": statistics.median(lens),
                   "master": ungapped(recs[0][1]), "recs": recs}

    #  Reference = the cluster at the feature's modal length with the most
    #  members. Modal length is taken over sequences, not over clusters, so one
    #  large aberrant cluster cannot define the mode.
    all_len = [l for v in info.values() for l in [v["med"]] * v["n"]]
    mode = statistics.median(all_len)
    cands = [f for f, v in info.items() if abs(v["med"] - mode) <= 0.1 * mode]
    if not cands:
        sys.exit("no cluster near the modal length %.0f" % mode)
    ref = max(cands, key=lambda f: info[f]["n"])
    print("  %s/%s: %d clusters, modal length %.0f, reference %s (n=%d, %d aa)"
          % (a.module, a.key, len(info), mode, os.path.basename(ref),
             info[ref]["n"], len(info[ref]["master"])))

    tmp = tempfile.mkdtemp(prefix="nterm.")
    rf = os.path.join(tmp, "ref.faa")
    open(rf, "w").write(">ref\n%s\n" % info[ref]["master"])

    flagged = []
    for f, v in sorted(info.items()):
        if f == ref: continue
        if v["med"] < mode + a.min_extension: continue
        qf = os.path.join(tmp, "q.faa")
        open(qf, "w").write(">q\n%s\n" % v["master"])
        r = subprocess.run(["blastp", "-query", rf, "-subject", qf, "-outfmt",
                            "6 qstart qend sstart send pident length"],
                           capture_output=True, text=True)
        if not r.stdout.strip(): continue
        qs, qe, ss, se, pid, ln = r.stdout.split("\n")[0].split("\t")
        off = int(ss) - int(qs)          # residues the query has before the reference starts
        if off < a.min_extension: continue
        #  The offset depends on which cluster happens to be the reference and
        #  on exactly where BLAST starts the alignment, so it can be a residue
        #  or two out. Look for a Met in a small window around it rather than
        #  demanding one at the exact position: with a 193 aa reference this
        #  cluster's start came out at 58, which is a Met, and with a 198 aa
        #  one at 59, which is not -- the same cluster, the same right answer.
        res, best = "?", None
        for d in (0, -1, 1, -2, 2, -3, 3):
            i = off + d
            if 0 <= i < len(v["master"]) and v["master"][i] == "M":
                best = i; break
        if best is not None:
            off, res = best, "M"
        elif off < len(v["master"]):
            res = v["master"][off]
        flagged.append((f, v, off, res, float(pid)))

    if not flagged:
        print("  no cluster starts upstream of the feature's other clusters.")
        return 0

    print("\n  N-TERMINAL EXTENSION (%d) -- these clusters start upstream of the"
          "\n  rest of the feature, and their profiles will anchor there:" % len(flagged))
    for f, v, off, res, pid in flagged:
        met = "Met" if res == "M" else "%s -- NOT a Met" % res
        print("    %-22s n=%-3d %4.0f aa (mode %.0f)  starts %d residues early"
              % (os.path.basename(f), v["n"], v["med"], mode, off))
        print("    %-22s reference aligns from its residue %d, which is %s, at %.0f%% identity"
              % ("", off + 1, met, pid))

    if not a.write:
        print("\n  nothing written (rerun with --write)")
        return 1

    for f, v, off, res, pid in flagged:
        if res != "M":
            print("  REFUSED %s: residue %d is %s, not a Met. The right start is a"
                  " curator's call." % (os.path.basename(f), off + 1, res))
            continue
        #  Map the ungapped offset onto an alignment column using the master row.
        master_row = v["recs"][0][1]
        n, col = 0, None
        for i, ch in enumerate(master_row):
            if ch not in "-.":
                n += 1
                if n == off + 1: col = i; break
        if col is None:
            print("  REFUSED %s: could not map offset to a column" % os.path.basename(f))
            continue
        with open(f, "w") as out:
            for h, s in v["recs"]: out.write("%s\n%s\n" % (h, s[col:]))
        print("  trimmed %s by %d columns; rebuild its PSSM before installing"
              % (os.path.basename(f), col))
    return 1

sys.exit(main())
