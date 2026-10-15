#!/usr/bin/env python3
"""Build the ORF1ab transcript-editing reference set for an astrovirus module.

ORF1ab is reached by a -1 ribosomal frameshift at a slippery heptamer, so it is
not collinear with the genome and carries no PSSM. `get_transcript_edited_features.pl`
blastn's a NUCLEOTIDE reference carrying the extra base against the genome, fills
the gap in the subject from the query and translates, so every genome gets its
own protein rather than a copy of the reference.

The construction here is **self-validating by design**. Rather than deriving the
insertion point from the heptamer's phase and hoping, it searches a window of
candidate insertion points and keeps only the one whose translation reproduces
the genome's own deposited nsp1ab protein byte for byte. A genome that cannot
reproduce its own protein donates no reference. That is the skill's "self round
trip" requirement applied at build time rather than afterwards.

The astrovirus heptamer is AAAAAAC (A_AAA_AAC); it is reported, not assumed --
it is visible in the junction window of most Mamastrovirus genomes -- but the
search does not depend on finding it, and records whether it was there.

Trailing `*`: nsp1ab really terminates at a stop codon, so the stop is kept in
the reference and the protein carries a trailing `*`, per the convention in
references/special-features.md (V, W, ORF1ab, ssGP, TF keep it; NSP12 and Gn,
which end at a cleavage, do not).
"""
import argparse, collections, os, subprocess, sys

CODON = {}
_b = "TCAG"
_a = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_i = 0
for _x in _b:
    for _y in _b:
        for _z in _b:
            CODON[_x + _y + _z] = _a[_i]; _i += 1


def translate(nt):
    return "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def read_fasta(path):
    h, s, out = None, [], []
    for line in open(path):
        if line.startswith(">"):
            if h is not None:
                out.append((h, "".join(s)))
            h = line[1:].strip(); s = []
        else:
            s.append(line.strip())
    if h is not None:
        out.append((h, "".join(s)))
    return out


def contig(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--proteins", required=True,
                    help="FASTA of deposited nsp1ab proteins, id = genome id")
    ap.add_argument("--out", required=True)
    ap.add_argument("--window", type=int, default=60,
                    help="nt either side of the HSP junction to search")
    ap.add_argument("--heptamer", default="AAAAAAC")
    a = ap.parse_args()

    prots = read_fasta(a.proteins)
    kept, stats = [], collections.Counter()
    for gid, prot in prots:
        gid = gid.split()[0]
        cp = os.path.join(a.contigs, gid + ".fna")
        if not os.path.exists(cp):
            stats["no contig"] += 1; continue
        g = contig(cp)
        open("/tmp/_te_q.faa", "w").write(">q\n%s\n" % prot)
        subprocess.run(["makeblastdb", "-in", cp, "-dbtype", "nucl", "-out", "/tmp/_te_db"],
                       capture_output=True)
        r = subprocess.run(["tblastn", "-query", "/tmp/_te_q.faa", "-db", "/tmp/_te_db",
                            "-comp_based_stats", "0", "-evalue", "1e-5",
                            "-outfmt", "6 qstart qend sstart send sframe bitscore"],
                           capture_output=True, text=True)
        hs = [x.split("\t") for x in r.stdout.splitlines()]
        hs = [h for h in hs if len(h) == 6 and int(h[4]) > 0]      # plus strand only
        if len(hs) < 2:
            stats["<2 plus-strand HSPs"] += 1; continue
        hs.sort(key=lambda x: int(x[0]))
        q1s, s1s = int(hs[0][0]), int(hs[0][2])
        junction = int(hs[0][3])
        #  0-based genome index of the first base of the protein's first codon
        start = (s1s - 1) - 3 * (q1s - 1)
        if start < 0:
            stats["start before contig"] += 1; continue
        target = prot.rstrip("*")
        hit = None
        for p in range(max(start + 3, junction - a.window), min(len(g), junction + a.window)):
            #  -1 slip: the ribosome re-reads one base, which is the same
            #  reading as duplicating that base in the template.
            cand = g[start:p] + g[p - 1] + g[p:]
            tr = translate(cand)
            stop = tr.find("*")
            if stop < 0:
                continue
            if tr[:stop] == target:
                hit = (p, cand[:(stop + 1) * 3])
                break
        if hit is None:
            stats["no insertion reproduces the protein"] += 1; continue
        p, ref = hit
        if set(ref) - set("ACGT"):
            stats["ambiguous base in construct"] += 1; continue
        hept = a.heptamer in g[max(0, p - 10):p + 10]
        stats["heptamer %s in window" % a.heptamer if hept else "no heptamer"] += 1
        kept.append((gid, ref, p, hept))
        stats["kept"] += 1

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        for gid, ref, p, hept in kept:
            fh.write(">%s\n%s\n" % (gid, ref))
    print("  references written: %d -> %s" % (len(kept), a.out))
    for k, v in stats.most_common():
        print("     %-44s %d" % (k, v))
    if kept:
        L = sorted(len(r) // 3 - 1 for _g, r, _p, _h in kept)
        print("  protein length from the references: median %d, range %d-%d"
              % (L[len(L) // 2], L[0], L[-1]))


main()
