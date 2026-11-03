#!/usr/bin/env python3
"""Emit check_literature_gaps.py-shaped probe files from probes.json.

check_literature_gaps.py wants a FLAT list of {name, pattern, note, source}
and scores it against collections/synonyms.tsv. probes.json is organised by
module and keyed by feature, which is the right shape for a curator to read
and the wrong shape for the checker. This converts one to the other rather
than keeping two hand-maintained copies that can drift apart.

A probe with no `pattern` cannot be scored, and a probe with no `source` is
reported UNSOURCED and fails the run -- both are errors here, not warnings.
"""
import json, os, sys

src = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "probes.json"))
bad = []
for mod, blk in src.items():
    if not isinstance(blk, dict) or "probes" not in blk:
        continue
    out = []
    for p in blk["probes"]:
        pm = p.get("pmids") or {}
        pmids = sorted({x for v in pm.values() for x in v})
        if not p.get("pattern"):
            bad.append("%s/%s: no pattern" % (mod, p["name"]))
            continue
        if not pmids:
            bad.append("%s/%s: no PMID" % (mod, p["name"]))
            continue
        out.append({"name": p["name"],
                    "pattern": p["pattern"],
                    "note": p.get("note", ""),
                    "source": "PMID_claude_generated " +
                              ", ".join(str(x) for x in pmids)})
    os.makedirs(mod, exist_ok=True)
    with open(os.path.join(mod, "probes.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("%-18s %d probes" % (mod, len(out)))
if bad:
    print("\nNOT EMITTED -- every probe needs a pattern and a PMID:")
    for b in bad:
        print("  " + b)
    sys.exit(1)
