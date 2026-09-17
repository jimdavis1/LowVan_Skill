#!/usr/bin/env python3
"""Collect the BV-BRC string -> annotation string aggregation.

Reads synonyms.tsv (emitted by build_collections.py) and joins it to the
annotation strings in <Taxon>_Viral_PSSM.json. Writes agg.json,
unbinned.json and typos.json for gen_collapse.py.

Run from the working directory that holds synonyms.tsv and the JSON.
"""
import csv, collections, json, sys, os, glob


def _pick_json(W):
    """Find the module JSON in this workdir.

    The family name used to be hardcoded as Rhabdoviridae_Viral_PSSM.json, so
    this script only ever worked for the family it was written for. Take
    whatever <Taxon>_Viral_PSSM.json is present, and say so if there are none
    or several.
    """
    c = sorted(glob.glob(os.path.join(W, "*_Viral_PSSM.json")))
    if not c:
        sys.exit("no *_Viral_PSSM.json in %s" % os.path.abspath(W))
    if len(c) > 1:
        sys.exit("several module JSONs in %s: %s -- keep one"
                 % (os.path.abspath(W), ", ".join(os.path.basename(x) for x in c)))
    return c[0]
if len(sys.argv) > 1 and sys.argv[1] in ("-h", "--help"):
    sys.exit("usage: collect_synmap.py [workdir]   # writes agg.json, unbinned.json, typos.json")
W = sys.argv[1] if len(sys.argv) > 1 else "."
rows = list(csv.DictReader(open(os.path.join(W, "synonyms.tsv")), delimiter="\t"))
J = json.load(open(_pick_json(W)))
key2anno, key2gene = {}, {}
for T, v in J.items():
    for k, e in v["features"].items():
        key2anno.setdefault(k, e["anno"]); key2gene.setdefault(k, e.get("gene_symbol", ""))

agg = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, set()]))
for r in rows:
    k = r["bins_to"]
    if not k: continue
    s = r["annotation"] or "(empty string)"
    agg[k][s][0] += int(r["count"]); agg[k][s][1].add(r["genus"] or "(no genus)")
out = []
for k, srcs in agg.items():
    out.append(dict(key=k, anno=key2anno.get(k, "(not in JSON)"), gene=key2gene.get(k, ""),
                    total=sum(v[0] for v in srcs.values()), nstrings=len(srcs),
                    srcs=sorted(([s, v[0], sorted(v[1])] for s, v in srcs.items()),
                                key=lambda x: -x[1])))
out.sort(key=lambda x: -x["total"])
json.dump(out, open("agg.json", "w"), indent=1)

un = collections.defaultdict(lambda: [0, set()])
for r in rows:
    if r["bins_to"]: continue
    s = r["annotation"] or "(empty string)"
    un[s][0] += int(r["count"]); un[s][1].add(r["genus"] or "(no genus)")
ul = sorted(([s, v[0], sorted(v[1])] for s, v in un.items()), key=lambda x: -x[1])
json.dump(ul, open("unbinned.json", "w"), indent=1)

def ed(a, b):
    if abs(len(a) - len(b)) > 3: return 9
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]
typos = []
for r in out:
    if len(r["srcs"]) < 2: continue
    top = r["srcs"][0][0]
    for s, c, g in r["srcs"][1:]:
        d = ed(s.lower(), top.lower())
        if 0 < d <= 2 and c >= 3: typos.append((r["key"], top, s, c, d, g))
typos.sort(key=lambda x: -x[3])
json.dump(typos, open("typos.json", "w"), indent=1)
print("agg.json %d targets, unbinned.json %d strings, typos.json %d variants"
      % (len(out), len(ul), len(typos)))
