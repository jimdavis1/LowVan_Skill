#!/usr/bin/env python3
"""Build Cardiovirus L* (overlapping ORF) and 2B* (-1 PRF transframe) per genome.

    special_cardio.py <ann_dir> <contig_dir> <out_dir>

L*  Kong & Roos 1991 (PMID 2033677); Chen et al. 1995 (PMID 7585219).
    An alternative AUG in the leader coding region, out of frame with the
    polyprotein, in the TO subgroup of TMEV (DA, BeAn) and absent from GDVII.
    Called here as the longest ORF that starts at an AUG within the first
    150 nt of the polyprotein, in a frame other than the polyprotein's, and
    runs >= 100 codons. Collinear, so it becomes an ordinary PSSM feature.

2B* Loughran et al. 2011 (PMID 22025686); Finch et al. 2015 (PMID 26063423).
    -1 programmed ribosomal frameshift at G_GUU_UUU just past the 2A/2B
    junction: the N-terminal 11-12 residues of 2B, then the -1 frame to its
    stop. EMCV: 128-129 aa. Theilovirus: the same site, a very short ORF.
    Slippery heptamer X_XXY_YYZ: the lone X is the THIRD base of a codon
    (phase 2). On the slip the P-site tRNA re-pairs one base back, so base h4
    is read twice: construct = genome[..h4] + h4 + genome[h4+1..], read in the
    0 frame from the 2B start to the first stop, stop codon kept (the product
    terminates there; trailing-* convention). N-terminus from the 2A/2B
    junction -> mat_peptide.
"""
import collections, glob, os, re, sys
from Bio.Seq import Seq

ANN, CON, OUT = sys.argv[1:4]
os.makedirs(OUT, exist_ok=True)

def contig(g):
    s = "".join(l.strip() for l in open(os.path.join(CON, g + ".fna")) if not l.startswith(">"))
    return s.upper()

def calls(t):
    d = {}
    for l in open(t):
        c = l.rstrip("\n").split("\t")
        #  feature.tbl: 6 symbol, 7 start, 8 stop, 9 strand
        if len(c) > 9 and c[9] == "+":
            d.setdefault(c[6], (int(c[7]), int(c[8])))
    return d

def tr(s):
    return str(Seq(s[: len(s) // 3 * 3]).translate())

lstar, bstar, log = [], [], collections.Counter()
for t in sorted(glob.glob(os.path.join(ANN, "*.feature.tbl"))):
    g = os.path.basename(t)[:-12]
    cl = calls(t)
    if "POLY" not in cl and "Ldr" not in cl: log["no POLY/Ldr call"] += 1; continue
    try: s = contig(g)
    except FileNotFoundError: log["no contig"] += 1; continue
    p0 = (cl.get("Ldr") or cl["POLY"])[0] - 1          # 0-based polyprotein start
    # ---- L*
    best = None
    for i in range(p0 + 1, min(p0 + 150, len(s) - 3)):
        if (i - p0) % 3 == 0 or s[i:i + 3] != "ATG": continue
        prot = tr(s[i:]).split("*")[0]
        if len(prot) >= 100 and (best is None or len(prot) > len(best[1])):
            best = (i, prot)
    if best:
        lstar.append((g, best[0] + 1, best[1])); log["L* ORF >=100 codons"] += 1
    else:
        log["no L* ORF"] += 1
    # ---- 2B*
    if "2B" not in cl: log["no 2B call"] += 1; continue
    b0 = cl["2B"][0] - 1
    hits = [m.start() for m in re.finditer(r"(?=GGTTTT[TC])", s[b0:b0 + 90])]
    #  the heptamer starts with the lone G at codon phase 2
    hits = [h for h in hits if h % 3 == 2]
    if not hits: log["no in-phase G_GUU_UUY in first 30 codons of 2B"] += 1; continue
    h1 = b0 + hits[0]; h4 = h1 + 3
    con = s[b0:h4 + 1] + s[h4] + s[h4 + 1:]
    prot = tr(con)
    if "*" not in prot: log["2B* no stop"] += 1; continue
    n = prot.index("*") + 1
    nt = con[: n * 3]
    if re.search("[^ACGT]", nt): log["2B* ambiguous base"] += 1; continue
    bstar.append((g, b0 + 1, hits[0] // 3, nt, prot[:n])); log["2B* built"] += 1

sp = {l.split("\t")[0]: l.rstrip("\n").split("\t")[4]
      for l in open(os.path.expanduser("~/lowvan_scratch/picorna/Picorna.full.tsv"))}
with open(os.path.join(OUT, "Lstar.faa"), "w") as f:
    for g, pos, p in lstar: f.write(">%s|Lstar|nt%d %s\n%s\n" % (g, pos, sp.get(g, "?"), p))
with open(os.path.join(OUT, "2Bstar_constructs.fna"), "w") as f, open(os.path.join(OUT, "2Bstar.faa"), "w") as fa:
    for g, pos, cod, nt, p in bstar:
        f.write(">%s|2Bstar|2Bstart=%d|slip_codon=%d %s\n%s\n" % (g, pos, cod, sp.get(g, "?"), nt.lower()))
        fa.write(">%s %s\n%s\n" % (g, sp.get(g, "?"), p))
print(dict(log))
def summ(name, rows, k):
    by = collections.defaultdict(list)
    for r in rows: by[sp.get(r[0], "?")].append(len(r[k]))
    print(name)
    for s_, L in sorted(by.items(), key=lambda x: -len(x[1])):
        L.sort(); print("   %-40s n=%3d  len %d-%d median %d" % (s_[:40], len(L), L[0], L[-1], L[len(L) // 2]))
summ("L* (aa)", lstar, 2)
summ("2B* (aa incl stop)", bstar, 4)
print("slip codon within 2B:", collections.Counter(r[2] for r in bstar).most_common(5))
