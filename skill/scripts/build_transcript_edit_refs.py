#!/usr/bin/env python3
"""Build a transcript-editing reference set for a frameshifted polyprotein.

A frameshift product is not collinear with the genome, so it carries no PSSM.
`get_transcript_edited_features.pl` blastn's a NUCLEOTIDE reference carrying the
extra base against the genome, fills the gap in the subject from the query and
translates, so each genome gets its own protein rather than a copy of the
reference. References live at
`<LOWVAN_DATA_DIR>/Transcript-Editing/<Module>/<FEAT>.fasta`.

**The construction is self-validating.** Rather than deriving the insertion
point from the slippery heptamer's phase and hoping, it scans candidate
insertion points and keeps only one whose translation reproduces a protein the
genome itself deposits. A genome that cannot reproduce its own protein donates
no reference. That is the skill's "self round trip" applied at build time.

Two validators, because the pool for the first is usually small:

  --fused     the deposited FUSION protein (nsp1ab, ORF1ab). Strongest: the
              construct must equal it exactly.
  --n-part    the deposited pre-shift protein (nsp1a, ORF1a), and
  --c-part    the deposited post-shift protein (ORF1b, the polymerase).
              Used together: the construct must START with the pre-shift
              protein's first --prefix residues, which are shared because both
              products initiate at the same Met, and END with the post-shift
              protein exactly, which is what pins the insertion point.

The second validator matters. Building Astroviridae, Mamastrovirus deposits 217
fused proteins and gains 49 from pairs, but **Avastrovirus deposits 24 and gains
131** -- a six-fold denser reference set, against a correction gate that demands
95% identity.

Locating the protein in the genome is done by exact match in a three-frame
translation, not by tblastn: it is faster, needs no BLAST database per genome,
and an exact match is a stronger anchor than an alignment.

The trailing `*`: a product that really terminates at a stop keeps the stop
codon in the reference (V, W, ORF1ab, ssGP, TF). One that ends at a cleavage
does not (NSP12, Gn). This builder is for the first kind and keeps it.
"""
import argparse, collections, os, sys

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
                out.append((h.split()[0], "".join(s)))
            h = line[1:].strip(); s = []
        else:
            s.append(line.strip())
    if h is not None:
        out.append((h.split()[0], "".join(s)))
    return out


def contig(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def locate(genome, protein, probe=30):
    """0-based nt index of the protein's first codon, by exact 3-frame match."""
    frag = protein[:probe]
    for fr in range(3):
        aa = translate(genome[fr:])
        i = aa.find(frag)
        if i >= 0:
            return fr + i * 3
    return -1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--fused", help="FASTA of deposited fusion proteins, id = genome id")
    ap.add_argument("--n-part", help="FASTA of deposited pre-shift proteins")
    ap.add_argument("--c-part", help="FASTA of deposited post-shift proteins")
    ap.add_argument("--out", required=True)
    ap.add_argument("--window", type=int, default=400,
                    help="nt downstream of the start to search for the slip (default 400 "
                         "before the pre-shift protein's own stop)")
    ap.add_argument("--prefix", type=int, default=200,
                    help="residues of the pre-shift protein the construct must start with")
    ap.add_argument("--heptamer", default="", help="optional, reported not required")
    a = ap.parse_args()

    if not a.fused and not (a.n_part and a.c_part):
        sys.exit("need --fused, or both --n-part and --c-part")

    fused = dict(read_fasta(a.fused)) if a.fused else {}
    npart = dict(read_fasta(a.n_part)) if a.n_part else {}
    cpart = dict(read_fasta(a.c_part)) if a.c_part else {}
    gids = list(fused) + [g for g in sorted(set(npart) & set(cpart)) if g not in fused]

    kept, stats = [], collections.Counter()
    for gid in gids:
        cp = os.path.join(a.contigs, gid + ".fna")
        if not os.path.exists(cp):
            stats["no contig"] += 1; continue
        g = contig(cp)
        anchor = fused.get(gid) or npart.get(gid)
        start = locate(g, anchor)
        if start < 0:
            stats["protein not found in a 3-frame translation"] += 1; continue

        target = fused.get(gid, "").rstrip("*")
        want_c = cpart.get(gid, "").rstrip("*")
        want_n = (npart.get(gid, "") or "")[:a.prefix]

        hit = None
        for p in range(start + 3, min(len(g), start + 3 + a.window * 10)):
            cand = g[start:p] + g[p - 1] + g[p:]
            tr = translate(cand)
            stop = tr.find("*")
            if stop < 0:
                continue
            prot = tr[:stop]
            if target:
                ok = prot == target
            else:
                ok = bool(want_c) and prot.endswith(want_c) and \
                     (not want_n or prot.startswith(want_n))
            if ok:
                hit = (p, cand[:(stop + 1) * 3])
                break
        if hit is None:
            stats["no insertion reproduces the deposited protein"] += 1; continue
        p, ref = hit
        if set(ref) - set("ACGT"):
            stats["ambiguous base in construct"] += 1; continue
        if a.heptamer:
            stats["heptamer present" if a.heptamer in g[max(0, p - 12):p + 12]
                  else "heptamer absent"] += 1
        stats["kept (%s)" % ("fused" if target else "n+c pair")] += 1
        kept.append((gid, ref))

    d = os.path.dirname(a.out)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(a.out, "w") as fh:
        for gid, ref in kept:
            fh.write(">%s\n%s\n" % (gid, ref))
    print("  references written: %d -> %s" % (len(kept), a.out))
    for k, v in stats.most_common():
        print("     %-46s %d" % (k, v))
    if kept:
        L = sorted(len(r) // 3 - 1 for _g, r in kept)
        print("  protein length from the references: median %d, range %d-%d"
              % (L[len(L) // 2], L[0], L[-1]))


main()
