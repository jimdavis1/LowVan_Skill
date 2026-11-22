"""Junction chemistry and adjacency over an Anno_GTOs directory.

Rewritten for the GTO pipeline: the earlier version read the flat
`ann/*.feature.tbl` that evaluate_coverage's non-special path writes. The full
script annotates through New-annotate-viral-taxon.pl instead, so the calls live
in Anno_GTOs/<g>.anno.gto and the family_assignments carry the feature key.
"""
import json, glob, os, sys, collections

ANNO = sys.argv[1]

def key(f):
    for fa in f.get("family_assignments", []):
        if fa and fa[0] == "LOWVAN" and "." in fa[1]:
            return fa[1].split(".")[1]
    return None

per = {}
for p in sorted(glob.glob(os.path.join(ANNO, "*.anno.gto"))):
    g = os.path.basename(p).split(".anno.gto")[0]
    try:
        gto = json.load(open(p))
    except Exception:
        continue
    d = {}
    for f in gto.get("features", []):
        k = key(f)
        loc = f.get("location") or []
        if not k or not loc:
            continue
        c, start, strand, ln = loc[0][0], int(loc[0][1]), loc[0][2], int(loc[0][3])
        beg, end = (start, start + ln - 1) if strand == "+" else (start - ln + 1, start)
        prot = f.get("protein_translation", "")
        #  keep the LONGEST call per key; a duplicate is its own finding
        if k not in d or (end - beg) > (d[k][1] - d[k][0]):
            d[k] = (beg, end, prot)
    per[g] = d

def junction(a, b):
    flush = tot = 0
    hist = collections.Counter()
    for g, d in per.items():
        if a in d and b in d:
            tot += 1
            gap = d[b][0] - d[a][1] - 1
            hist[gap] += 1
            if gap == 0:
                flush += 1
    return flush, tot, hist

def endswith(k, aa):
    n = t = 0
    for g, d in per.items():
        if k in d and d[k][2]:
            t += 1
            if d[k][2].rstrip("*").endswith(aa):
                n += 1
    return n, t

print("  %d genome(s) with calls\n" % len(per))
print("  junction adjacency")
for a, b in [("NSP1A","NSP1B"),("NSP2","NSP3"),("NSP3","NSP4"),("NSP4","NSP5"),
             ("NSP5","NSP6"),("NSP6","NSP7A"),("NSP7A","NSP7B"),("NSP7B","NSP8"),
             ("NSP9","NSP10"),("NSP10","NSP11"),("NSP11","NSP12")]:
    fl, tot, hist = junction(a, b)
    if not tot:
        continue
    top = ", ".join("%+d:%d" % (k, v) for k, v in hist.most_common(4))
    print("    %-6s -> %-6s  %4d/%-4d flush (%5.1f%%)   %s"
          % (a, b, fl, tot, 100.0 * fl / tot, top))

print("\n  C-terminal residue (3CLSP cuts after Glu; PLP2 after Gly)")
for k, aa in [("NSP1A","M"),("NSP1B","G"),("NSP2","G"),("NSP3","E"),("NSP4","E"),
              ("NSP5","E"),("NSP6","E"),("NSP7A","E"),("NSP7B","E"),("NSP8","E")]:
    n, t = endswith(k, aa)
    if t:
        print("    %-6s ends %s  %4d/%-4d (%5.1f%%)" % (k, aa, n, t, 100.0 * n / t))
