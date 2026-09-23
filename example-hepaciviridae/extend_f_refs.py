#!/usr/bin/env python3
"""Extend the F reference set to genomes that carry no F annotation.

The 844 references built from annotated records are a genotype-1a set: of the
337 that carry a genotype label, 319 are 1a and genotypes 2-6 have one or two
each. HCV genotypes differ by ~30% at the nucleotide level, so a 1a reference
cannot correct a genotype-3 genome through a 95% identity gate -- which is why
64% of HCV genomes came back as "Uncorrected F protein encoding region"
rather than a protein.

3,550 near-complete HCV genomes have no F reference, and they span 1a, 1b,
2a/2b/2c, 3a, 4a/4d and 6a/6k/6m.

Those genomes have no annotated F protein, so the byte-identity validation
used for the first 844 is unavailable. The construction RULE is what carries
over, and it was established on those 844: delete one base at the A-rich
slippery site in codon 10-11 of core, read the +1 frame to the stop. Four
checks stand in for the missing protein, and a genome failing any is skipped:

  1. core's own N-terminal decapeptide must be present and in frame
  2. the deleted base must sit in a run of >= 4 A's
  3. the product must translate with no internal stop
  4. its length must fall inside the range the 844 observed, 121-178 aa

This is weaker evidence than reproducing a known protein and is recorded as
such. It is not a claim that these genomes express F; it is a reference that
lets the annotator report the F reading frame where one is found.
"""
import collections, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CODON = {}; _b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            CODON[a + b + c] = _aa[i * 16 + j * 4 + k]

def tr(s): return "".join(CODON.get(s[i:i+3].upper(), "x") for i in range(0, len(s)-2, 3))

def rd(p):
    d = {}; h = None; s = []
    for line in open(p):
        if line.startswith(">"):
            if h: d[h] = "".join(s)
            h = line[1:].strip(); s = []
        else: s.append(line.strip())
    if h: d[h] = "".join(s)
    return d

existing = rd(os.path.join(HERE, "Transcript-Editing/Hepaciviridae/F.fasta"))
have = {h.split("|")[1].split()[0] for h in existing}
#  Constraints taken from the 844 protein-validated references, not guessed.
#  The first version of this script used range(8,14) and A{4,} and took the
#  first match; codon 8 always won, and 2,723 of 2,788 references cut two
#  codons early. They translated to MSTNPKPQ*E*K where core is MSTNPKPQ*R*K --
#  residue 9 wrong in every one -- and sat on a 3-A run where the validated
#  site has five. Nothing downstream would have caught it: the products had no
#  internal stop and fell inside the length range.
LO, HI = 121, 178                      # length range the 844 observed
CODONS = (10, 11)                      # where all 844 validated edits sit
SLIP = re.compile(r"A{5,}")            # 90% of validated sites; A{4,} admits the wrong one

#  core's decapeptide, from the reference set itself rather than assumed
deca = collections.Counter()
for h, s in existing.items():
    deca[tr(s)[:10]] += 1
#  core's real decapeptides, taken from the core collection rather than from
#  the reference set, so a bad reference cannot widen the accepted set
_core = rd(os.path.join(HERE, "collections/Hepaciviridae/C.fasta"))
STARTS = {w for w, n in collections.Counter(v[:10] for v in _core.values()).most_common() if n >= 5}
print("core decapeptides seen in the existing references: %d forms" % len(STARTS))

names = {}
for l in open(os.path.join(HERE, "Hepaci.id_name_gs")):
    p = l.rstrip("\n").split("\t")
    if len(p) > 1: names.setdefault(p[0], p[1])

gen = {}
for l in open(os.path.join(HERE, "f_pool_genomes.tsv")):
    p = l.rstrip("\n").split("\t")
    if len(p) >= 2: gen[p[0]] = p[1].upper()

out = []; why = collections.Counter(); lens = []
for gid, g in sorted(gen.items()):
    if gid in have: why["already has a reference"] += 1; continue
    st = None
    for fr in range(3):
        t = tr(g[fr:])
        for w in STARTS:
            i = t.find(w)
            if i >= 0: st = fr + 3*i; break
        if st is not None: break
    if st is None: why["core decapeptide not found"] += 1; continue
    best = None
    for cod in CODONS:
        for off in range(0, 3):
            pos = st + 3*cod + off
            if not SLIP.search(g[max(0,pos-6):pos+7]): continue
            cand = g[st:pos] + g[pos+1:]
            t = tr(cand); stop = t.find("*")
            if stop < 0: continue
            if not (LO <= stop <= HI): continue
            if set(cand[:3*(stop+1)]) - set("ACGT"): continue
            #  the product must still begin with core's own decapeptide; this
            #  is the check that catches an edit placed too early
            if t[:10] not in STARTS: continue
            best = cand[:3*(stop+1)]; break
        if best: break
    if not best: why["no clean +1 ORF at codon 10/11 on a 5-A run"] += 1; continue
    lens.append(len(best)//3 - 1)
    out.append((">F|%s Hepaciviridae F protein [%s | %s]"
                % (gid, names.get(gid, "Hepacivirus C"), gid), best.lower()))

p = os.path.join(HERE, "Transcript-Editing/Hepaciviridae/F.fasta")
with open(p, "a") as fh:
    for h, s in out:
        fh.write(h + "\n")
        for i in range(0, len(s), 60): fh.write(s[i:i+60] + "\n")
lens.sort()
print("added %d reference(s); set is now %d" % (len(out), len(existing) + len(out)))
if lens: print("  product length aa: min %d median %d max %d" % (lens[0], lens[len(lens)//2], lens[-1]))
for k, v in why.most_common(): print("  %-36s %d" % (k, v))
