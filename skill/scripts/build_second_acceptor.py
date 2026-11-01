#!/usr/bin/env python3
"""Build Splice-Variants references for a second product of a SHARED DONOR.

build_splice_refs.py places both exons with tblastn. That fails when exon B is
too short to score: Amdoparvovirus NS3 is 71 aa with a ~61 aa first exon, so
exon B is about 10 residues and tblastn never reports it. The NS3 query then
yields 7 references from 300 genomes while the NS2 query -- same donor, longer
second exon -- yields 227.

The parvovirus NS transcripts share ONE donor and differ only in the acceptor
(measured on Amdoparvovirus, 1 Nov 2026: NS2 and NS3 both cut at SD last nt
179 in 1961072.4, and accept at 1839 and 1532 respectively). So take the donor
from the reference set that built cleanly, then scan every downstream AG for
the acceptor that reproduces the second product.

    python3 build_second_acceptor.py \
        --from Amdoparvovirus/NS2.fasta --query NS3.faa \
        --contigs cov/Amdoparvovirus/contigs --out Amdoparvovirus/NS3.fasta \
        --name NS3

A reference is kept only if the spliced translation starts at M, carries no
internal stop, lands within --len-tol of the query length, and is at least
--min-id identical to the query over a global alignment of equal-length
sequences. Everything else is dropped and counted; the drops are the evidence
that the acceptor is real rather than the first AG that happens to fit.
"""
import argparse, os, re, sys

B = 'TCAG'
AAs = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
TAB = {a + b + c: AAs[i] for i, (a, b, c) in enumerate(
    (x, y, z) for x in B for y in B for z in B)}


def tr(s):
    return "".join(TAB.get(s[i:i + 3], 'X') for i in range(0, len(s) - 2, 3))


def read_fasta(path):
    out, cur = [], None
    for line in open(path):
        if line.startswith(">"):
            cur = [line[1:].rstrip(), []]
            out.append(cur)
        elif cur:
            cur[1].append(line.strip())
    return [(h, "".join(s)) for h, s in out]


def ident(a, b):
    """Fraction identical over the shorter sequence, no gaps.

    The two products come from the same first exon, so they are already in
    frame with each other; an alignment would buy nothing and would let a
    shifted second exon score as a match."""
    n = min(len(a), len(b))
    return sum(1 for i in range(n) if a[i] == b[i]) / float(n) if n else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="src", required=True,
                    help="reference fasta that already built, same donor")
    ap.add_argument("--query", required=True, help="protein FASTA, one sequence")
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--name", default="FEAT")
    ap.add_argument("--len-tol", type=float, default=0.10)
    ap.add_argument("--min-id", type=float, default=0.55)
    ap.add_argument("--max-intron", type=int, default=3000)
    ap.add_argument("--intron", default="",
                    help="comma-separated intron lengths to accept, +/- "
                         "--intron-tol. Without this the scan takes the best "
                         "exon-B identity anywhere downstream, and on "
                         "Amdoparvovirus NS3 that chose a 428 nt intron in 135 "
                         "of 173 genomes while every reference that matched "
                         "its own genome's deposit had 1352. Measure the "
                         "intron on the references that built from two HSPs, "
                         "or that reproduce a deposit, and pass it here.")
    ap.add_argument("--intron-tol", type=int, default=3)
    ap.add_argument("--donor-window", type=int, default=0,
                    help="also try donors this many nt either side of the "
                         "source donor. 0 keeps the shared donor, which is "
                         "the point of this script; use a window when the two "
                         "products differ at the donor as well.")
    ap.add_argument("--deposits",
                    help="the feature's collection fasta, headers "
                         "<genome_id>|<FEATURE>_<n>. A candidate that equals "
                         "this genome's own deposited protein is accepted "
                         "outright.")
    a = ap.parse_args()

    q = read_fasta(a.query)[0][1].rstrip("*")
    lo, hi = len(q) * (1 - a.len_tol), len(q) * (1 + a.len_tol)

    want = [int(x) for x in a.intron.split(",") if x.strip()]
    deposits = {}
    if a.deposits:
        for h, s in read_fasta(a.deposits):
            deposits.setdefault(h.split("|")[0], set()).add(s)
    kept, drops = [], {"no AG": 0, "start": 0, "stop": 0, "length": 0,
                       "identity": 0, "intron": 0}
    seen = set()
    for hdr, seq in read_fasta(a.src):
        gid = hdr.split()[0]
        if gid in seen:
            continue
        seen.add(gid)
        sd0 = int(re.search(r"SD:\d+-\d+;(\d+)", hdr).group(1))
        mine = deposits.get(gid, ())

        best = None
        for sd in range(sd0 - a.donor_window, sd0 + a.donor_window + 1):
            if sd < 3 or sd >= len(seq):
                continue
            head = seq[:sd]
            #  every AG in the window, not the first that fits: the first AG
            #  after the donor reproduced NS1's own frame in early testing.
            for i in range(sd, min(len(seq) - 2, sd + a.max_intron)):
                if seq[i:i + 2] != "AG":
                    continue
                if want and not any(abs((i + 3) - sd - 1 - w) <= a.intron_tol
                                    for w in want):
                    drops["intron"] += 1
                    continue
                prot = tr(head + seq[i + 2:])
                stop = prot.find("*")
                if stop < 0:
                    continue
                prot = prot[:stop]
                if not prot.startswith("M"):
                    drops["start"] += 1
                    continue
                #  A product that IS this genome's own deposited protein needs
                #  no further argument, and is preferred over any scored
                #  candidate. Protoparvovirus NS2 splices at a GC--AG donor,
                #  which build_splice_refs will not choose while a GT sits in
                #  the same window, so this is the only rule that finds it.
                if prot in mine:
                    best = (2.0, sd, i + 3, prot)
                    break
                if not (lo <= len(prot) <= hi):
                    drops["length"] += 1
                    continue
                #  Score the SECOND exon only. Exon A is shared with the
                #  product this donor already built, so whole-protein identity
                #  is ~0.77 for almost any acceptor that yields the right
                #  length, and picking on it chose a different intron in
                #  nearly every genome (428, 608, 346 nt against the 1333 nt
                #  the two-HSP build measured). Exon B is short -- about 10
                #  residues for Amdoparvovirus NS3 -- but it is the only part
                #  that discriminates.
                k = sd // 3
                pid = ident(prot[k:], q[k:])
                if pid < a.min_id:
                    drops["identity"] += 1
                    continue
                if best is None or pid > best[0]:
                    #  SA first nt is 1-based and is the base after the AG
                    best = (pid, sd, i + 3, prot)
            if best and best[0] == 2.0:
                break
        if best is None:
            drops["no AG"] += 1
            continue
        pid, sd, sa, prot = best
        kept.append((gid, sd, sa, seq, prot, pid))

    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        for gid, sd, sa, seq, prot, pid in kept:
            #  the slice came from the OTHER product and runs to ITS stop, so
            #  cut it back to this one's: exon B is 3*len - sd nt, then the
            #  stop codon. Leaving the tail on would hand the caller a
            #  reference whose 3' end belongs to a different CDS.
            end = (sa - 1) + (3 * len(prot) - sd) + 3
            fh.write(">%s SD:%d-%d;%d SA:%d-%d;%d\n%s\n"
                     % (gid, sd - 10, sd + 10, sd, sa - 11, sa + 3, sa,
                        seq[:end]))

    lens = sorted(len(p) for _, _, _, _, p, _ in kept)
    ids = sorted(pid for *_, pid in kept if pid < 2.0)
    conf = sum(1 for *_, pid in kept if pid == 2.0)
    print("  %s references: %d written (%d reproduce their own genome's "
          "deposit) -> %s" % (a.name, len(kept), conf, a.out))
    print("  dropped: " + ", ".join("%d %s" % (v, k) for k, v in drops.items()))
    if ids:
        print("  query %d aa; spliced product min %d med %d max %d; "
              "identity to query min %.2f med %.2f"
              % (len(q), lens[0], lens[len(lens) // 2], lens[-1],
                 ids[0], ids[len(ids) // 2]))


if __name__ == "__main__":
    sys.exit(main())
