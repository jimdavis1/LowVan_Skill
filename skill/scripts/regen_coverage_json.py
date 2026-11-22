#!/usr/bin/env python3
"""Re-derive a coverage json from GTOs that are already on disk.

`evaluate_coverage.py` has no report-only mode, and the annotation is the
expensive half: 612 genomes took hours, and two defects in the REPORTING half
(the module taken from a transcript edit's feature key, and refs_used globbing
only *.feature.tbl) should not cost another run of it.  This imports the fixed
reporter and points it at the existing Anno_GTOs.

Usage: regen_coverage_json.py <work dir> <out.json>
"""
import importlib.util
import os
import sys

#  the sibling copy, so the kit is self-contained; LOWVAN_KIT still wins.
_here = os.path.dirname(os.path.abspath(__file__))
_ec = os.path.join(os.environ["LOWVAN_KIT"], "skill", "scripts", "evaluate_coverage.py") \
      if os.environ.get("LOWVAN_KIT") else os.path.join(_here, "evaluate_coverage.py")
spec = importlib.util.spec_from_file_location("ec", _ec)
ec = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ec)

work, out = sys.argv[1], sys.argv[2]
exd = os.path.join(work, "ex_gto")

ex, seqs = [], {}
for fn in sorted(os.listdir(exd)):
    g = fn.rsplit(".", 1)[0] if "." in fn else fn
    for suf in (".fna", ".fasta", ".dna"):
        if g.endswith(suf):
            g = g[: -len(suf)]
    s = []
    for line in open(os.path.join(exd, fn)):
        if line.startswith(">"):
            s.append([])
        elif s:
            s[-1].append(line.strip())
    ex.append(g)
    seqs[g] = ["".join(x) for x in s]
print("  %d exemplar(s) read from %s" % (len(ex), exd))
ec.report_gto(ex, seqs, work, out)
