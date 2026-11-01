#!/usr/bin/env python3
"""Build Splice-Variants references from the DEPOSITS themselves.

build_splice_refs.py works from one query protein and has to decide, for every
other genome, where the junction is. That decision is what goes wrong: it
nudges onto the nearest canonical GT..AG, and for Protoparvovirus NS2 the
nearest GT is not the donor (the donor is GC..AG), while for Protoparvovirus
NS3 the nudge window itself moved the donor 24 nt and produced a 104 aa
product against a 107 aa deposit.

Here there is no decision to make. For each deposited protein, tblastn it
against ITS OWN genome, take the two plus-strand HSPs that run in query order,
and place the junction so that the spliced translation IS the deposit. A
reference is written only when it reproduces the deposit exactly, so every
reference is confirmed by construction. The output is bounded by how many
genomes deposited the protein, which is the honest bound: see
[[splice-ref-junction-is-not-proof]].

    python3 build_refs_from_deposits.py \\
        --deposits Protoparvovirus/collections/Protoparvovirus/NS3.fasta \\
        --contigs cov/Protoparvovirus/contigs \\
        --out Protoparvovirus/Splice-Variants/Protoparvovirus/NS3.fasta \\
        --name NS3
"""
import argparse, collections, os, subprocess, sys, tempfile

B = 'TCAG'
AAs = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
TAB = {a + b + c: AAs[i] for i, (a, b, c) in enumerate(
    (x, y, z) for x in B for y in B for z in B)}


def tr(s):
    return "".join(TAB.get(s[i:i + 3], 'X') for i in range(0, len(s) - 2, 3))


COMP = str.maketrans("ACGTRYKMBVDHNacgtrykmbvdhn", "TGCAYRMKVBHDNtgcayrmkvbhdn")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    out, cur = [], None
    for line in open(path):
        if line.startswith(">"):
            cur = [line[1:].rstrip(), []]
            out.append(cur)
        elif cur:
            cur[1].append(line.strip())
    #  contig fasta from the dump is soft-masked lower case; the
    #  codon table is upper case, so every translation came back as
    #  a run of X and no junction could ever match a deposit.
    return [(h, "".join(s).upper()) for h, s in out]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deposits", required=True,
                    help="collection fasta, headers <genome_id>|<FEAT>_<n>")
    ap.add_argument("--contigs", required=True, help="dir of <genome_id>.fna")
    ap.add_argument("--out", required=True)
    ap.add_argument("--name", default="FEAT")
    ap.add_argument("--slack", type=int, default=18,
                    help="how far, in nt, each HSP end may be moved to put the "
                         "junction where the deposit says it is. tblastn trims "
                         "to the last well-matching codon, so the reported "
                         "boundary is a few codons short of the real one.")
    a = ap.parse_args()

    kept, drops = [], collections.Counter()
    seen = set()
    for hdr, prot in read_fasta(a.deposits):
        gid = hdr.split("|")[0]
        prot = prot.rstrip("*")
        if (gid, prot) in seen:
            continue
        seen.add((gid, prot))
        fna = os.path.join(a.contigs, gid + ".fna")
        if not os.path.exists(fna):
            drops["no contig file"] += 1
            continue
        contigs = dict((h.split()[0], s) for h, s in read_fasta(fna))

        with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as fh:
            fh.write(">q\n%s\n" % prot)
            qf = fh.name
        try:
            r = subprocess.run(
                ["tblastn", "-query", qf, "-subject", fna, "-seg", "no",
                 "-evalue", "10", "-outfmt",
                 "6 qstart qend sstart send sseqid bitscore"],
                capture_output=True, text=True)
        finally:
            os.unlink(qf)
        rows = [l.split("\t") for l in r.stdout.strip().split("\n") if l]
        hsps = [(int(x[0]), int(x[1]), int(x[2]), int(x[3]), x[4],
                 float(x[5])) for x in rows]
        #  One strand only: a spliced CDS does not change strand between its
        #  exons, so a hit on the other strand belongs to another gene. An
        #  ambisense genome puts half its genes on the minus strand -- the
        #  Scindoambidensovirus capsid among them -- so when the minus-strand
        #  hits win, flip the contig and work in its coordinates throughout.
        plus = [h for h in hsps if h[3] > h[2]]
        minus = [h for h in hsps if h[3] < h[2]]
        if sum(h[5] for h in minus) > sum(h[5] for h in plus):
            contigs = {k: rc(v) for k, v in contigs.items()}
            hsps = [(h[0], h[1], len(contigs[h[4]]) - h[2] + 1,
                     len(contigs[h[4]]) - h[3] + 1, h[4], h[5]) for h in minus]
        else:
            hsps = plus
        if len(hsps) < 2:
            drops["fewer than two HSPs on one strand"] += 1
            continue

        hit = None
        for i, A in enumerate(hsps):
            for Bb in hsps:
                if Bb is A or Bb[4] != A[4]:
                    continue
                #  exon A first in the query AND first on the genome
                if not (A[0] < Bb[0] and A[2] < Bb[2] and A[3] < Bb[3]):
                    continue
                seq = contigs.get(A[4])
                if not seq:
                    continue
                #  exon A must start at the deposit's own residue 1
                start = A[2] - 3 * (A[0] - 1)
                if start < 1:
                    continue
                #  step 1, not 3: the junction need not be codon-aligned. The
                #  deposited Protoparvovirus NS3 of 10796.157 splits residue
                #  86 across the intron -- tblastn reports q86 as the last
                #  residue of exon A and the first of exon B -- and a step of
                #  3 could never place it.
                for da in range(-a.slack, a.slack + 1):
                    for db in range(-a.slack, a.slack + 1):
                        ea = A[3] + da          # last nt of exon A
                        sb = Bb[2] + db         # first nt of exon B
                        if not (start <= ea < sb <= len(seq)):
                            continue
                        exonA = seq[start - 1:ea]
                        exonB = seq[sb - 1:]
                        p = tr(exonA + exonB)
                        stop = p.find("*")
                        if stop < 0:
                            continue
                        if p[:stop] != prot:
                            continue
                        #  the slice the caller is given runs from the start
                        #  codon to this product's own stop
                        end = sb - 1 + (3 * len(prot) - len(exonA)) + 3
                        hit = (A[4], start, end, len(exonA), sb - start + 1)
                        break
                    if hit:
                        break
                if hit:
                    break
            if hit:
                break
        if hit is None:
            drops["no junction reproduces the deposit"] += 1
            continue
        cid, start, end, sd, sa = hit
        kept.append((gid, contigs[cid][start - 1:end], sd, sa, len(prot)))

    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        for gid, seq, sd, sa, _ in kept:
            fh.write(">%s SD:%d-%d;%d SA:%d-%d;%d\n%s\n"
                     % (gid, sd - 10, sd + 10, sd, sa - 11, sa + 3, sa, seq))

    intr = collections.Counter(sa - sd - 1 for _, _, sd, sa, _ in kept)
    lens = sorted(L for *_, L in kept)
    print("  %s references: %d written, every one reproducing its own genome's "
          "deposit -> %s" % (a.name, len(kept), a.out))
    if drops:
        print("  dropped: " + ", ".join("%d %s" % (v, k)
                                        for k, v in drops.most_common()))
    if kept:
        print("  product min %d med %d max %d; intron %s"
              % (lens[0], lens[len(lens) // 2], lens[-1],
                 ", ".join("%d nt x%d" % kv for kv in intr.most_common(4))))


if __name__ == "__main__":
    sys.exit(main())
