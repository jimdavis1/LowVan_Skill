#!/usr/bin/env python3
"""Recover full-length nsP3 from the genome, between nsP2's end and nsP4's start.

WHY
nsP3's C-terminus is the nsP2 protease cleavage site, which lies SIX CODONS PAST
the opal. The nsP3 profiles cannot reach it in Chikungunya: the dominant profile
(NSP3.1, 514 columns) gives up at the opal with 7 columns unused, because those
columns sit in nsP3's hypervariable C-terminal domain and carry too little
information to pay tblastn's stop penalty. 90 of the 92 opal-carrying CHIKV
genome types are affected.

The interval between the two FLANKING mature peptides is nsP3 by definition, and
both flanks are called reliably in exactly these genomes. So take it from the
genome rather than from a profile that cannot reach the end.

THE STOP CHARACTER -- this is the trap
The recovered protein contains the readthrough stop. It cannot be written as
'*':
  * mafft silently rewrites '*' to '-', turning the residue into a GAP and
    shifting the column's meaning with no error;
  * psiblast refuses an MSA containing '*' outright
    ("CAlnReader::GetSeqEntry(): Seq_entry is not available until after Read()").
So the stop is written as 'X'. That is also the honest character for TRAINING
data: across the genus this position is UGA in opal strains and CGA/AGA (Arg),
CAA (Gln) or TGT (Cys) in others, so a low-information column is what the
profile should have there. It leaves the flanking conserved columns to carry the
extension, which is the whole point.

The emitted ANNOTATION still carries '*' -- that is the annotator's business and
is unaffected by this file.
"""
import collections, glob, os, sys

ANN = sys.argv[1] if len(sys.argv) > 1 else "coverage_eval_is/ann"
EX  = "coverage_eval/ex"
LO, HI = 410, 660                      # NSP3's declared bounds, widened a little

T = {}
_b = "TCAG"; _a = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
i = 0
for x in _b:
    for y in _b:
        for z in _b: T[x + y + z] = _a[i]; i += 1
def tr(s): return "".join(T.get(s[j:j + 3], "X") for j in range(0, len(s) - 2, 3))
def rc(s): return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))

meta = {}
for l in open("Toga.genome_meta.tsv"):
    p = l.rstrip("\n").split("\t")
    if len(p) >= 7: meta[p[0]] = dict(sp=(p[6].strip() or "?"), name=p[2])

stats = collections.Counter(); out = []; bysp = collections.Counter()
for t in sorted(glob.glob(os.path.join(ANN, "*.feature.tbl"))):
    if not os.path.getsize(t): continue
    g = os.path.basename(t).replace(".feature.tbl", "")
    f = collections.defaultdict(list)
    for line in open(t, errors="replace"):
        c = line.rstrip("\n").split("\t")
        if len(c) < 10: continue
        try: lo, hi = int(c[7]), int(c[8])
        except ValueError: continue
        f[c[6]].append((lo, hi, c[9]))
    if len(f.get("nsP2", [])) != 1 or len(f.get("nsP4", [])) != 1:
        stats["flanks not called exactly once"] += 1; continue
    (a1, a2, sa), (b1, b2, sb) = f["nsP2"][0], f["nsP4"][0]
    if sa != sb: stats["flanks on opposite strands"] += 1; continue
    p = os.path.join(EX, g + ".fna")
    if not os.path.exists(p): stats["no sequence"] += 1; continue
    s = "".join(x.strip() for x in open(p) if not x.startswith(">")).upper()
    if sa == "+":
        lo, hi = max(a1, a2) + 1, min(b1, b2) - 1
        nt = s[lo - 1:hi]
    else:
        lo, hi = max(b1, b2) + 1, min(a1, a2) - 1
        nt = rc(s[lo - 1:hi])
    if hi <= lo: stats["inverted interval"] += 1; continue
    if len(nt) % 3: stats["interval not a multiple of 3"] += 1; continue
    aa = tr(nt)
    if not (LO <= len(aa) <= HI): stats["length outside %d-%d" % (LO, HI)] += 1; continue
    ns = aa.count("*")
    if ns > 1: stats["more than one stop -- broken genome"] += 1; continue
    if aa.count("X") / len(aa) > 0.01: stats["too ambiguous"] += 1; continue
    if ns == 1:
        aa = aa.replace("*", "X"); stats["recovered WITH readthrough (stop -> X)"] += 1
    else:
        stats["recovered, no readthrough"] += 1
    bysp[meta.get(g, {}).get("sp", "?")] += 1
    out.append((g, meta.get(g, {}).get("name", "?"), meta.get(g, {}).get("sp", "?"), aa, ns))

os.makedirs("Extracted", exist_ok=True)
with open("Extracted/Togaviridae.NSP3.faa", "w") as fh:
    fh.write("; full-length nsP3 recovered from the nsP2-nsP4 interval.\n"
             "; The readthrough stop is written X -- mafft turns '*' into a gap\n"
             "; and psiblast refuses an MSA containing one. See extract_nsp3.py.\n")
    for g, nm, sp, aa, ns in out:
        fh.write(">extracted|%s %s | %s | %s\n%s\n"
                 % (g, nm, sp, "readthrough" if ns else "no readthrough", aa))
print("%d sequences written to Extracted/Togaviridae.NSP3.faa\n" % len(out))
for k, v in stats.most_common(): print("   %-44s %5d" % (k, v))
print("\n   top species recovered:", dict(bysp.most_common(6)))
L = sorted(len(a[3]) for a in out)
print("   lengths: min %d  median %d  max %d" % (L[0], L[len(L) // 2], L[-1]))
