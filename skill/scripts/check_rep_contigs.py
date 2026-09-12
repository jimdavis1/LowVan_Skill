#!/usr/bin/env python3
"""Check that the representative contigs actually route the taxon's genomes.

Routing is step one of annotation: the incoming contig is BLASTn'd against
Viral-Rep-Contigs/ and the best hit decides which PSSM directory gets used. A
genome with no hit above the minimum contig bitscore is rejected outright and
gets no annotation at all, no matter how good the PSSMs are.

This is easy to under-provision, because BLASTn is a *nucleotide* search. Within
a genus it is generous; across genera of a divergent group it can return nothing
at all. Rice yellow stunt virus had zero hits against five Betarhabdovirinae rep
contigs -- one per genus was not enough, and no PSSM ever ran.

Give it the accessions you expect the module to cover:

    python3 check_rep_contigs.py --repdir Rep-Contigs --acc-file accs.txt
    python3 check_rep_contigs.py --repdir Rep-Contigs --acc-file accs.txt --mcb 150

Anything reported UNROUTED needs another rep contig near it.
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def fetch_fasta(accs):
    url = "%s/efetch.fcgi?%s" % (EUTILS, urllib.parse.urlencode(
        {"db": "nuccore", "id": ",".join(accs), "rettype": "fasta", "retmode": "text"}))
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=120) as fh:
                return fh.read().decode("utf8", "replace")
        except Exception:
            if attempt == 2:
                raise
            time.sleep(3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repdir", required=True, help="directory of *.dna rep contigs")
    ap.add_argument("--acc-file", help="accessions to test, one per line")
    ap.add_argument("--acc", help="comma-separated accessions")
    ap.add_argument("--mcb", type=float, default=150.0,
                    help="minimum contig bitscore the annotator requires (default 150)")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    accs = []
    if args.acc:
        accs += [a.strip() for a in args.acc.split(",") if a.strip()]
    if args.acc_file:
        accs += [l.split("#")[0].strip() for l in open(args.acc_file)
                 if l.split("#")[0].strip()]
    if not accs:
        sys.exit("give --acc or --acc-file")

    reps = sorted(f for f in os.listdir(args.repdir) if f.endswith(".dna"))
    if not reps:
        sys.exit("no *.dna files in %s" % args.repdir)
    print("%d rep contig files, %d test accessions\n" % (len(reps), len(accs)))

    tmp = tempfile.mkdtemp(prefix="reprout.")
    db = os.path.join(tmp, "rep.fa")
    # tag every sequence with the rep-contig file it came from, so a hit names
    # the file (and therefore the module) that would be chosen
    with open(db, "w") as out:
        for r in reps:
            for line in open(os.path.join(args.repdir, r), errors="replace"):
                if line.startswith(">"):
                    out.write(">%s|%s" % (r, line[1:]))
                else:
                    out.write(line)
    subprocess.run(["makeblastdb", "-in", db, "-dbtype", "nucl",
                    "-out", os.path.join(tmp, "rep")], capture_output=True)

    q = os.path.join(tmp, "q.fa")
    with open(q, "w") as out:
        out.write(fetch_fasta(accs))

    res = subprocess.run(
        ["blastn", "-query", q, "-db", os.path.join(tmp, "rep"),
         "-outfmt", "6 qseqid sseqid pident bitscore", "-evalue", "1",
         "-max_target_seqs", "20", "-num_threads", str(args.threads)],
        capture_output=True, text=True)

    best = {}
    for line in res.stdout.splitlines():
        qid, sid, pid, bits = line.split("\t")
        acc = qid.split(".")[0]
        b = float(bits)
        if acc not in best or b > best[acc][0]:
            best[acc] = (b, sid.split("|")[0], float(pid))

    routed, weak, unrouted = [], [], []
    for acc in accs:
        a = acc.split(".")[0]
        if a not in best:
            unrouted.append((a, None, 0.0, 0.0))
        else:
            b, rep, pid = best[a]
            (routed if b >= args.mcb else weak).append((a, rep, b, pid))

    print("  %-12s %-30s %9s %7s" % ("accession", "routes to", "bits", "%id"))
    for a, rep, b, pid in routed:
        print("  %-12s %-30s %9.0f %6.1f%%" % (a, rep, b, pid))
    for a, rep, b, pid in weak:
        print("  %-12s %-30s %9.0f %6.1f%%   BELOW --mcb %g" % (a, rep, b, pid, args.mcb))
    for a, _r, _b, _p in unrouted:
        print("  %-12s %-30s %9s %7s" % (a, "-- NO BLASTn HIT --", "0", "-"))

    print("\n  routed %d/%d" % (len(routed), len(accs)))
    if weak or unrouted:
        print("  %d genome(s) would be REJECTED before any PSSM runs." % (len(weak) + len(unrouted)))
        print("  Add a rep contig from the same genus as each one, then re-check.")
        return 1
    print("  every test genome routes above the bitscore floor")
    return 0


if __name__ == "__main__":
    sys.exit(main())
