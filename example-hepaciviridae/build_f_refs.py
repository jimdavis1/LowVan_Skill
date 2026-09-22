#!/usr/bin/env python3
"""Build the F (ARFP / core+1) transcript-editing reference set for Hepaciviridae.

F is the +1 frameshift product of the core coding region: the ribosome
initiates at core's AUG, reads ten codons in frame 0, then shifts +1 between
codons 9 and 11 over a run of ten A residues and reads to a stop.

    Xu 2001      PMID 11447125   frameshift, EMBO J
    Walewski 2001 PMID 11350035  the antigen, RNA
    Baril 2005   PMID 15755749   argues non-AUG initiation instead

The mechanism is contested; the construction here does not depend on
resolving it. Measured over 906 records, 844 (93%) are reproduced by deleting
exactly one base at codon 10-11, so the edit is one base whichever way the
ribosome gets there, and the annotator reproduces the observed protein either
way.

This is the OPPOSITE direction to NS1'. A -1 slip makes the reference one base
LONGER than the genome; a +1 slip makes it one base SHORTER. The reference is
therefore built by deletion, and get_transcript_edited_features.pl handles it
by dropping the subject base where the query gaps.

Construction is protein-driven, as for NS1': locate core's start, walk the
0-frame until the genome stops agreeing with the known protein, then try every
single-base deletion in a window and keep only one that reproduces the protein
byte-identically.
"""
import collections, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, c1 in enumerate(_b):
    for j, c2 in enumerate(_b):
        for k, c3 in enumerate(_b):
            CODON[c1 + c2 + c3] = _aa[i * 16 + j * 4 + k]


def tr(s):
    return "".join(CODON.get(s[i:i + 3].upper(), "x") for i in range(0, len(s) - 2, 3))


def rd(p):
    d = {}; h = None; s = []
    for line in open(p):
        if line.startswith(">"):
            if h: d[h] = "".join(s)
            h = line[1:].strip(); s = []
        else:
            s.append(line.strip())
    if h: d[h] = "".join(s)
    return d


prot = rd(os.path.join(HERE, "collections/Hepaciviridae/F.fasta"))
gen = {}
for l in open(os.path.join(HERE, "f_genomes.tsv")):
    p = l.rstrip("\n").split("\t")
    if len(p) >= 2:
        gen[p[0]] = p[1].upper()
names = {}
for l in open(os.path.join(HERE, "Hepaci.id_name_gs")):
    p = l.rstrip("\n").split("\t")
    if len(p) > 1:
        names.setdefault(p[0].split(".")[0], p[1])

#  The slippery run. Xu 2001 is often quoted as "a stretch of 10A", which
#  reads as a homopolymer and is not what the genomes show: the site is an
#  A-RICH region, CAAAGAAAAACCAAA in 435 of 844 records, whose longest
#  homopolymer is five. Testing for A{7,} matched 1 record of 844 and would
#  have looked like the site being absent. A{5,} matches 90% and A{3,} 99.7%.
SLIP = re.compile(r"A{5,}")

out, failed, stats = [], [], collections.Counter()
where = collections.Counter()

for hdr, aa in sorted(prot.items()):
    fid = hdr.split()[0]
    m = re.match(r"fig\|(\d+\.\d+)\.", fid)
    if not m:
        failed.append("unparseable id"); continue
    gid = m.group(1)
    g = gen.get(gid)
    if not g:
        failed.append("no genome"); continue

    start = None
    for f in range(3):
        t = tr(g[f:])
        i = t.find(aa[:10])
        if i >= 0:
            start = f + 3 * i; break
    if start is None:
        failed.append("N-terminus not found"); continue

    zero = 0
    while zero < len(aa) and CODON.get(g[start + 3 * zero: start + 3 * zero + 3]) == aa[zero]:
        zero += 1
    if zero < 5:
        failed.append("0-frame agreement <5 aa"); continue

    #  every single-base deletion in a window, keep those that reproduce the
    #  protein exactly, prefer one sitting in the poly-A slippery run
    cands = []
    for back in range(0, 8):
        cut = start + 3 * (zero - back)
        if cut <= start:
            continue
        for off in range(0, 3):
            pos = cut + off
            cand = g[start:pos] + g[pos + 1:]
            t = tr(cand)
            stop = t.find("*")
            if stop < 0 or t[:stop] != aa:
                continue
            win = g[max(0, pos - 8):pos + 8]
            cands.append((bool(SLIP.search(win)), zero - back, cand[:3 * (stop + 1)], pos))
    if not cands:
        failed.append("no single-base deletion reproduces the protein"); continue
    inslip = [c for c in cands if c[0]]
    pick = (inslip or cands)[0]
    if not inslip:
        stats["chosen without a poly-A run at the site"] += 1
    ref = pick[2]
    if set(ref) - set("ACGT"):
        failed.append("reference would contain an ambiguous base"); continue

    where[pick[1]] += 1
    sp = names.get(gid.split(".")[0], "?")
    out.append((">F|%s Hepaciviridae F protein [%s | %s]" % (gid, sp, gid), ref.lower()))
    stats["built"] += 1

d = os.path.join(HERE, "Transcript-Editing/Hepaciviridae")
os.makedirs(d, exist_ok=True)
p = os.path.join(d, "F.fasta")
with open(p, "w") as fh:
    for h, s in out:
        fh.write(h + "\n")
        for i in range(0, len(s), 60):
            fh.write(s[i:i + 60] + "\n")

print("built %d reference(s) of %d protein(s) -> %s" % (len(out), len(prot), p))
L = sorted(len(s) for _h, s in out)
if L:
    print("  nt length: min %d median %d max %d" % (L[0], L[len(L) // 2], L[-1]))
print("  codon at which the frame shifts: %s" % dict(sorted(where.items())[:6]))
for k, v in stats.most_common():
    print("  %-44s %d" % (k, v))
if failed:
    print("\n%d not built:" % len(failed))
    for r, n in collections.Counter(failed).most_common():
        print("   %-46s %d" % (r, n))
