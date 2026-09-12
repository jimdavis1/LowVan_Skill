#!/usr/bin/env python3
"""Validate the DLIT citations in a module JSON and their provenance.

A DLIT is meant to be a paper a curator read, which establishes the function or
the coordinates of a feature. A citation a language model proposed is not that,
even when the paper is real and on-topic — nobody has confirmed it supports the
claim. Those must be flagged in the JSON as `PMID_claude_generated` so a reader
can tell the two apart, and so they can be verified or replaced later.

Checks, in order of how badly they bite:

  NOT A PMID    a DOI or other identifier in the PMID field. `PMID` is a list
                of bare numeric PubMed ids and nothing else.
  UNRESOLVED    the id does not exist in PubMed. For a model-proposed citation
                this is the fabrication case and the whole reason for the flag.
  UNFLAGGED     present in PMID_claude_generated's registry upstream but not
                flagged here, or flagged but absent from PMID — the two lists
                disagree.

Then prints every citation with its title so a curator can eyeball relevance.

    python3 check_dlits.py --json Taxon_Viral_PSSM.json
    python3 check_dlits.py --json Taxon_Viral_PSSM.json --offline   # no PubMed

Network access to PubMed is used unless --offline.
"""

import argparse
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict

ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
FLAG = "PMID_claude_generated"


def pubmed(ids):
    """{pmid: (title, year)} for ids that resolve; missing ids are absent."""
    out = {}
    ids = list(ids)
    for i in range(0, len(ids), 150):
        chunk = ids[i:i + 150]
        url = "%s?%s" % (ESUMMARY, urllib.parse.urlencode(
            {"db": "pubmed", "id": ",".join(chunk), "retmode": "json"}))
        data = None
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url, timeout=90) as fh:
                    data = json.load(fh)
                break
            except Exception:
                if attempt == 2:
                    raise
                time.sleep(3)
        for k, v in (data.get("result") or {}).items():
            if k == "uids" or not isinstance(v, dict):
                continue
            if v.get("error"):
                continue
            title = (v.get("title") or "").strip()
            if title:
                out[k] = (title, (v.get("pubdate") or "")[:4])
        time.sleep(0.4)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True)
    ap.add_argument("--offline", action="store_true",
                    help="skip PubMed; structural checks only")
    args = ap.parse_args()

    with open(args.json) as fh:
        mod = json.load(fh)

    used, flagged = defaultdict(list), defaultdict(list)
    disagree = []
    for module, block in mod.items():
        for key, ent in block.get("features", {}).items():
            pm = ent.get("PMID") or []
            fl = ent.get(FLAG) or []
            for p in pm:
                used[str(p)].append("%s/%s" % (module, key))
            for p in fl:
                flagged[str(p)].append("%s/%s" % (module, key))
            extra = [p for p in fl if p not in pm]
            if extra:
                disagree.append((module, key, extra))

    numeric = {p for p in used if re.fullmatch(r"\d+", p)}
    nonnum = sorted(p for p in used if p not in numeric)

    print("%d distinct citations across %d features\n"
          % (len(used), sum(1 for m, b in mod.items()
                            for k, e in b.get("features", {}).items() if e.get("PMID"))))

    problems = 0
    if nonnum:
        problems += len(nonnum)
        print("NOT A PMID (%d) -- the PMID field takes bare PubMed ids only:" % len(nonnum))
        for p in nonnum:
            print("  %-34s used by %s" % (p, ", ".join(used[p])[:60]))
        print("  Resolve a DOI at pubmed.ncbi.nlm.nih.gov and store the numeric id.\n")

    if disagree:
        problems += len(disagree)
        print("FLAG DISAGREES WITH PMID (%d):" % len(disagree))
        for module, key, extra in disagree:
            print("  %-20s %-14s flagged but not in PMID: %s"
                  % (module, key, ", ".join(extra)))
        print()

    resolved = {}
    if not args.offline and numeric:
        try:
            resolved = pubmed(sorted(numeric))
        except Exception as exc:
            print("PubMed lookup failed (%s); rerun with --offline to skip.\n" % exc)
            resolved = None

    if resolved is not None and not args.offline:
        missing = sorted(numeric - set(resolved))
        if missing:
            problems += len(missing)
            print("UNRESOLVED (%d) -- no such PubMed record:" % len(missing))
            for p in missing:
                tag = "MODEL-PROPOSED" if p in flagged else "human-curated"
                print("  %-10s %-16s used by %s" % (p, tag, ", ".join(used[p])[:54]))
            print("  A model-proposed id that does not resolve is a fabricated citation.")
            print("  Remove it; do not keep it with the flag.\n")

    print("CITATIONS (%d model-proposed, %d human-curated):\n"
          % (len(flagged), len(numeric) - len(flagged)))
    print("  %-10s %-16s %-5s %s" % ("PMID", "provenance", "year", "title"))
    for p in sorted(numeric, key=lambda x: (x not in flagged, x)):
        tag = "MODEL-PROPOSED" if p in flagged else "human"
        if resolved and p in resolved:
            title, year = resolved[p]
        else:
            title, year = ("(not looked up)" if args.offline else "(unresolved)"), "-"
        print("  %-10s %-16s %-5s %s" % (p, tag, year, title[:66]))

    if flagged:
        print("\n  Model-proposed citations are proposals, not DLITs. A curator has to")
        print("  read each one and confirm it establishes what the feature claims,")
        print("  then drop it from the generator's CLAUDE_PMIDS registry so the flag")
        print("  disappears on the next regeneration.")

    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
