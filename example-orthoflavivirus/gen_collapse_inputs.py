#!/usr/bin/env python3
"""Produce agg.json / unbinned.json / typos.json for the string-collapse artifact.

Reuses build_collections.py's own RULES and classify() rather than
reimplementing them: the point of the artifact is to show what the shipped
binning does, so a second copy of the rules that drifted would make the page
a description of something that was never run.

build_collections.py calls main() at import time, so the source is executed
with that trailing call stripped.

Where other modules report the GENERA a string came from, this reports
SPECIES: Orthoflavivirus is a single genus, so genus would be one value on
every row and say nothing.
"""
import collections, json, os, re, sys, difflib

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "build_collections.py")).read()
src = re.sub(r"(?m)^main\(\)\s*$", "", src)
ns = {"__name__": "_bc", "__file__": os.path.join(HERE, "build_collections.py")}
exec(compile(src, "build_collections.py", "exec"), ns)

RULES, NOT_MODELLED = ns["RULES"], ns["NOT_MODELLED"]
classify, compiled, T, MODULE = ns["classify"], ns["compiled"], ns["T"], ns["MODULE"]
decl = json.load(open(os.path.join(HERE, "%s_Viral_PSSM.json" % MODULE)))[MODULE]["features"]

species = {}
for line in open(os.path.join(HERE, "Flavi.id_name_gs")):
    p = line.rstrip("\n").split("\t")
    if len(p) > 1:
        species.setdefault(p[0].split(".")[0], p[1])
for line in open(os.path.join(HERE, "ortho_sample.meta.tsv")):
    p = line.rstrip("\n").split("\t")
    if len(p) > 4 and p[4]:
        species[p[0].split(".")[0]] = p[4]

rules = compiled(RULES)
notmod = [(p, w) for p, w in NOT_MODELLED]
md5s = [l.strip() for l in open(os.path.join(HERE, "%s.uniq.md5" % T))]
ann = [l.rstrip("\n").split("\t") for l in open(os.path.join(HERE, "%s.uniq.id_ann" % T))]
seq = [l.rstrip("\n").split("\t") for l in open(os.path.join(HERE, "%s.uniq.seq" % T))]

per_md5 = collections.Counter()
for l in open(os.path.join(HERE, "%s.id_md5" % T)):
    p = l.rstrip("\n").split("\t")
    if len(p) >= 2 and p[1].strip():
        per_md5[p[1].strip()] += 1


def sp_of(fid):
    m = re.match(r"fig\|(\d+)\.", fid or "")
    return species.get(m.group(1), "?") if m else "?"


bykey = collections.defaultdict(lambda: collections.Counter())
bykey_sp = collections.defaultdict(lambda: collections.defaultdict(set))
unb = collections.Counter()
unb_sp = collections.defaultdict(set)
notm = collections.Counter()
for i in range(len(md5s)):
    a = " ".join((ann[i][1] if len(ann[i]) > 1 else "").split())
    s = seq[i][1] if len(seq[i]) > 1 else ""
    if not s:
        continue
    fid = ann[i][0] if ann[i] else ""
    key, why = classify(a, rules, notmod)
    nfeat = per_md5.get(md5s[i], 1)
    sp = sp_of(fid)
    if key:
        bykey[key][a] += nfeat
        bykey_sp[key][a].add(sp)
    elif why:
        notm[(a, why)] += nfeat
    else:
        unb[a] += nfeat
        unb_sp[a].add(sp)

agg = []
for key in sorted(bykey, key=lambda k: -sum(bykey[k].values())):
    srcs = [[s, c, sorted(x for x in bykey_sp[key][s] if x != "?")]
            for s, c in bykey[key].most_common()]
    agg.append({"key": key,
                "anno": decl.get(key, {}).get("anno", key),
                "gene": decl.get(key, {}).get("gene_symbol", key),
                "total": sum(bykey[key].values()),
                "nstrings": len(bykey[key]),
                "srcs": srcs})

unbinned = [[s, c, sorted(x for x in unb_sp[s] if x != "?")] for s, c in unb.most_common()]

#  Near-duplicate strings inside one feature: the same protein named twice with
#  a typo or a punctuation difference. Reported against the most common form.
typos = []
for key in bykey:
    forms = [s for s, _ in bykey[key].most_common()]
    if len(forms) < 2:
        continue
    canon = forms[0]
    for s in forms[1:]:
        r = difflib.SequenceMatcher(None, canon.lower(), s.lower()).ratio()
        if r >= 0.85 and s.lower() != canon.lower():
            d = sum(1 for x in difflib.ndiff(canon.lower(), s.lower()) if x[0] in "+-")
            typos.append([key, canon, s, bykey[key][s], d,
                          sorted(x for x in bykey_sp[key][s] if x != "?")])
typos.sort(key=lambda r: -r[3])

for nm, obj in (("agg.json", agg), ("unbinned.json", unbinned), ("typos.json", typos)):
    json.dump(obj, open(os.path.join(HERE, nm), "w"), indent=1)
    print("wrote %-16s %d entries" % (nm, len(obj)))
print("\n%d features, %d distinct binned strings, %d unbinned strings, %d near-duplicates"
      % (len(agg), sum(a["nstrings"] for a in agg), len(unbinned), len(typos)))
print("binned features: %d" % sum(a["total"] for a in agg))
