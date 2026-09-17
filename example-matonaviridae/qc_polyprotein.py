#!/usr/bin/env python3
"""Cross-feature QC for a polyprotein taxon.

scripts/qc_cross_feature.py flags any two named features that share
near-identical sequence. In Matonaviridae that is the module's own architecture:
every mature peptide IS a substring of its precursor, so p150 is 100% identical
to the N-terminal 1301 residues of p200 and C to the first 300 of p110. Run
unfiltered it reports 13 MISLABELs, none of them mislabels, and a real one
would be invisible among them.

So: read its output, subtract the pairs the genome organisation predicts, and
report only what is left. Anything printed here is a genuine cross-feature
collision and needs a MISBINNED entry in build_collections.py.
"""
import collections, os, re, subprocess, sys

#  Every feature that is a substring of another, mapped to the family it sits
#  in. Two families, because the two rubivirus polyproteins never share
#  sequence -- they are separate ORFs, the second translated from a subgenomic
#  RNA:
#    nonstructural: p200 > p150, p90        (one cleavage, PMID 9557742)
#    structural:    p110 > C, E2, E1        (signalase, PMID 2273395)
FAMILY = {"NS":  {"P200", "P150", "P90"},
          "STR": {"P110", "C", "E2", "E1"}}
PARENT = {k: fam for fam, keys in FAMILY.items() for k in keys}
SKILL = os.path.join(os.environ.get("LOWVAN_KIT", ".."), "skill", "scripts")

out = subprocess.run([sys.executable, os.path.join(SKILL, "qc_cross_feature.py"),
                      "--workdir", "."], capture_output=True, text=True).stdout
rows = [m.groups() for m in
        (re.match(r"\s+Matonaviridae\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+([\d.]+)%\s+(\S+)", l)
         for l in out.splitlines()) if m]
expected, real = collections.Counter(), []
for binned, clus, master, ln, coll, pid, other in rows:
    base = coll.split(".")[0]
    if base == binned or (PARENT.get(binned) is not None
                          and PARENT.get(binned) == PARENT.get(base)):
        expected["%s / %s" % (binned, base)] += 1
    else:
        real.append((binned, clus, master, ln, coll, pid, other))

print("%d collision(s) reported by qc_cross_feature.py" % len(rows))
print("%d explained by the precursor/product architecture:" % sum(expected.values()))
for k, v in sorted(expected.items(), key=lambda x: -x[1]):
    print("   %4d  %s" % (v, k))
print("\n%d GENUINE cross-feature collision(s)%s"
      % (len(real), ":" if real else " -- clean"))
for r in real:
    print("   binned as %-8s cluster %-5s %s (%s aa) == %s at %s%% of %s" % r)
sys.exit(1 if real else 0)
