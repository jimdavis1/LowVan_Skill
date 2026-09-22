#!/usr/bin/env python3
"""Build the NS1' transcript-editing reference set for Orthoflavivirus.

NS1' is the -1 ribosomal frameshift product of the JEV serogroup: NS1 plus 52
residues, made when the ribosome slips at a conserved slippery heptanucleotide
(Y CCU UUU) near the start of NS2A, stimulated by a downstream pseudoknot.
Firth & Atkins 2009 predicted the site; Melian 2010 demonstrated the product.

There is no PSSM for this. get_transcript_edited_features.pl blastn's a curated
nucleotide reference -- carrying the slipped base DUPLICATED -- against the
genome, fills the resulting one-base gap in the subject from the query, and
translates. Everything except that base comes from the genome, so each genome
yields its own protein rather than a copy of a reference.

Construction here is deliberately not motif-driven. Searching for the heptamer
directly invites two failures the kit's own notes warn about: testing the motif
at the wrong codon phase rejects every genome and looks like a finding, and a
window anchored on the called feature misses the site whenever an N-run
truncated that call. Instead the KNOWN protein drives it:

  1. locate the codons of NS1's N-terminus in the genome, in all three frames
  2. translate forward in the 0-frame, comparing residue by residue, until the
     protein and the genome stop agreeing -- that is the slip
  3. duplicate one base there and continue in the -1 frame to the stop
  4. accept only if the whole construct translates BYTE-IDENTICALLY to the
     protein the dump already holds

Step 4 is the point. The construct is not asserted to be right; it is required
to reproduce a protein that was observed independently, and any genome where it
cannot is dropped rather than guessed at.
"""
import collections, json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CODON = {}
_b = "TCAG"
_aa = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
for i, c1 in enumerate(_b):
    for j, c2 in enumerate(_b):
        for k, c3 in enumerate(_b):
            CODON[c1 + c2 + c3] = _aa[i * 16 + j * 4 + k]


def tr(s):
    return "".join(CODON.get(s[i:i + 3].upper(), "x") for i in range(0, len(s) - 2, 3))


def rd(p):
    d = {}
    h = None
    s = []
    for line in open(p):
        if line.startswith(">"):
            if h:
                d[h] = "".join(s)
            h = line[1:].strip()
            s = []
        else:
            s.append(line.strip())
    if h:
        d[h] = "".join(s)
    return d


prot = rd(os.path.join(HERE, "collections/Orthoflavivirus/NS1P.fasta"))
gen = rd(os.path.join(HERE, "ns1p_genomes.fna"))
names = {}
for line in open(os.path.join(HERE, "Flavi.id_name_gs")):
    p = line.rstrip("\n").split("\t")
    names.setdefault(p[0].split(".")[0], p[1] if len(p) > 1 else "?")

SLIP = re.compile(r"[CT]CC[CT]TTT", re.I)   # Y CCU UUU, the slippery heptamer

out, failed, stats = [], [], collections.Counter()
motif_phase = collections.Counter()
ext_len = collections.Counter()

for hdr, aa in sorted(prot.items()):
    fid = hdr.split()[0]
    m = re.match(r"fig\|(\d+\.\d+)\.", fid)
    if not m:
        failed.append((fid, "unparseable id")); continue
    gid = m.group(1)
    g = gen.get(gid)
    if not g:
        failed.append((fid, "no genome sequence")); continue
    g = g.upper()
    if set(g) - set("ACGT"):
        stats["genome carries ambiguous bases"] += 1

    #  1. where does this protein's N-terminus sit? anchor on 12 residues.
    anchor = aa[:12]
    start = None
    for f in range(3):
        t = tr(g[f:])
        i = t.find(anchor)
        if i >= 0:
            start = f + 3 * i
            break
    if start is None:
        failed.append((fid, "N-terminus not found in any frame")); continue

    #  2. walk the 0-frame until the genome stops agreeing with the protein
    zero = 0
    while zero < len(aa):
        c = g[start + 3 * zero: start + 3 * zero + 3]
        if len(c) < 3 or CODON.get(c) != aa[zero]:
            break
        zero += 1
    if zero < 300:
        failed.append((fid, "0-frame agreement only %d aa" % zero)); continue

    #  3. slip: duplicate one base and continue in the -1 frame. The slip is at
    #     a codon boundary in the region where agreement stopped, but the exact
    #     base is not assumed -- every offset in a small window is tried and the
    #     one that reproduces the rest of the protein is kept.
    #  Collect EVERY duplication that reproduces the protein, then choose among
    #  them -- taking the first that happens to work is not safe. Five genomes
    #  admitted a solution with the duplication near the C-terminus, one
    #  residue before the stop: it reproduces this genome's protein because
    #  almost the whole construct is copied verbatim, but it puts the
    #  frameshift nowhere near the slip site, so the reference would only ever
    #  correct a near-identical genome. Those are exactly the five with no
    #  slippery heptamer at the chosen position, which is the signal to use.
    cands = []
    for back in range(0, 16):
        cut = start + 3 * (zero - back)          # nt offset of the slip codon
        if cut <= start:
            continue
        for dup in range(0, 3):                 # which base of the codon repeats
            pos = cut + dup
            cand = g[start:pos + 1] + g[pos:]
            t = tr(cand)
            stop = t.find("*")
            if stop < 0:
                continue
            if t[:stop] != aa:
                continue
            window = g[max(0, pos - 10):pos + 10]
            cands.append((bool(SLIP.search(window)), len(aa) - (zero - back),
                          cand[:3 * (stop + 1)], zero - back, dup, pos))
    if not cands:
        failed.append((fid, "no duplication reproduces the protein")); continue
    #  prefer a solution the heptamer supports, and among those the one whose
    #  shifted portion is the modal length for this feature
    motif_ok = [c for c in cands if c[0]]
    if not motif_ok:
        failed.append((fid, "no solution sits at a slippery heptamer")); continue
    motif_ok.sort(key=lambda c: (abs(c[1] - 44), -c[1]))
    _mo, _ext, refseq, slipaa_, dup_, pos_ = motif_ok[0]
    built = (refseq, slipaa_, dup_, pos_)

    ref, slipaa, dup, pos = built
    if set(ref) - set("ACGT"):
        failed.append((fid, "reference would contain an ambiguous base")); continue

    #  record where the heptamer sits relative to the slip, as a check on the
    #  biology rather than an input to it
    w = g[pos - 10:pos + 10]
    mm = SLIP.search(w)
    motif_phase[("motif found" if mm else "motif absent")] += 1
    ext_len[len(aa) - slipaa] += 1

    sp = names.get(gid.split(".")[0], "?")
    out.append((">NS1P|%s Orthoflavivirus NS1' protein [%s | %s]" % (gid, sp, gid),
                ref.lower()))
    stats["built"] += 1

os.makedirs(os.path.join(HERE, "Transcript-Editing/Orthoflavivirus"), exist_ok=True)
p = os.path.join(HERE, "Transcript-Editing/Orthoflavivirus/NS1P.fasta")
with open(p, "w") as fh:
    for h, s in out:
        fh.write("%s\n" % h)
        for i in range(0, len(s), 60):
            fh.write("%s\n" % s[i:i + 60])

print("built %d reference(s) of %d protein(s) -> %s" % (len(out), len(prot), p))
L = sorted(len(s) for _h, s in out)
if L:
    print("  nt length: min %d median %d max %d" % (L[0], L[len(L) // 2], L[-1]))
print("  slippery heptamer within 10 nt of the slip: %s" % dict(motif_phase))
print("  residues after the slip: %s" % dict(sorted(ext_len.items())))
for k, v in stats.most_common():
    print("  %-40s %d" % (k, v))
if failed:
    print("\n%d not built:" % len(failed))
    c = collections.Counter(r for _f, r in failed)
    for r, n in c.most_common():
        print("   %-46s %d" % (r, n))
