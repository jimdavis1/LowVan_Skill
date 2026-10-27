#!/usr/bin/env python3
"""What does the literature describe that BV-BRC does not call?

Collections are built from deposited annotation strings. A protein that is real
biology but was never annotated in the genomes you dumped is therefore INVISIBLE
to the build: it does not show up as a gap, as UNASSIGNED, or as a low count --
it simply never appears, and the module ships without it while every number on
the page looks clean.

Betaretrovirus Rem is the case that forced this check into the skill. It is a
described MMTV protein, a Rev-like nuclear export factor, and ZERO genomes in
that module carry a call for it. Nothing in the build would ever have said so.

So the check cannot be derived from the data -- it has to be driven by a list of
what the group is KNOWN to encode, written by whoever is building the module,
with a citation per entry. That list is the deliverable; this script only scores
it. Writing it is the point, and it is why this is a required step rather than a
report you can skim.

PROBE FILE (JSON), one entry per protein the literature describes:

    [
     {"name": "Rem", "pattern": "\\\\brem\\\\b|\\\\brej\\\\b",
      "note": "Rev-like nuclear export factor, doubly spliced",
      "source": "PMID 1656...."},
     {"name": "glyco-Gag", "pattern": "glyco-?gag|gpr80gag", ...}
    ]

Every entry needs a `source`. An entry without one is reported as UNSOURCED and
fails the run: "I think this protein exists" is not a literature check.

VERDICTS

    PRESENT       deposits found, and they bind to a declared feature
    RECOGNISED    deposits found, bound to NOT_MODELLED -- a deliberate call
    UNBOUND       deposits found but still UNASSIGNED -- a rule is missing
    BELOW_FLOOR   found, but fewer than --min-seqs, so no profile is possible
    ABSENT        no string anywhere in this module's vocabulary

ABSENT and UNBOUND are failures. ABSENT means the module is shipping without a
protein the group is known to encode; UNBOUND means your triage rules missed
vocabulary that is right there in the dump. Both need either a rule or a written
decision in the gaps file, which --write-gaps emits.
"""
import argparse, collections, json, os, re, sys

RESERVED = ("UNASSIGNED", "DISCARDED", "NOT_MODELLED")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--module", required=True)
    ap.add_argument("--probes", required=True,
                    help="JSON list of {name, pattern, note, source}")
    ap.add_argument("--min-seqs", type=int, default=3,
                    help="below this a collection cannot support a profile")
    ap.add_argument("--write-gaps", default=None,
                    help="write a markdown gaps record to this path")
    ap.add_argument("--verify-pmids", action="store_true",
                    help="resolve every PMID in the probe file against PubMed and "
                         "print its title. A remembered PMID is usually WRONG: of "
                         "three written from memory on 27 Oct 2026, two resolved to "
                         "unrelated papers -- an Italian cardiology trial and an SV40 "
                         "T-antigen study -- and both looked perfectly plausible as "
                         "bare numbers. Verify, or leave the field empty and fix it later.")
    ap.add_argument("--allow-absent", default="",
                    help="comma-separated probe names that are known-absent and "
                         "already decided; they report ABSENT but do not fail")
    a = ap.parse_args()

    syn = os.path.join(a.workdir, "collections", "synonyms.tsv")
    if not os.path.exists(syn):
        sys.exit("no %s -- build collections first" % syn)
    rows = []
    with open(syn) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        if head[:4] != ["count", "annotation", "genus", "bins_to"]:
            sys.exit("synonyms.tsv has the wrong schema: %s" % "/".join(head))
        for l in fh:
            f = l.rstrip("\n").split("\t")
            if len(f) >= 4:
                rows.append((int(f[0]), f[1], f[3]))

    probes = json.load(open(a.probes))
    allow = set(x.strip() for x in a.allow_absent.split(",") if x.strip())
    unsourced = [p["name"] for p in probes if not (p.get("source") or "").strip()]

    print("  %s: %d probe(s) against %d distinct strings\n" % (a.module, len(probes), len(rows)))
    results, fail = [], 0
    for p in probes:
        rx = re.compile(p["pattern"], re.I)
        hits = [(n, s, k) for n, s, k in rows if rx.search(s)]
        tot = sum(n for n, _, _ in hits)
        dests = collections.Counter(k for _, _, k in hits)
        if not hits:
            verdict = "ABSENT"
        elif tot < a.min_seqs:
            verdict = "BELOW_FLOOR"
        elif all(k == "UNASSIGNED" for k in dests):
            verdict = "UNBOUND"
        elif set(dests) <= {"NOT_MODELLED", "DISCARDED"}:
            verdict = "RECOGNISED"
        else:
            verdict = "PRESENT"
        bad = (verdict in ("ABSENT", "UNBOUND")) and p["name"] not in allow
        fail += bad
        results.append(dict(p, verdict=verdict, deposits=tot,
                            dests=dict(dests), failed=bad))
        print("  %-10s %-26s %4d deposit(s)  %s%s"
              % (verdict, p["name"][:26], tot,
                 ", ".join("%s(%d)" % (k, v) for k, v in dests.most_common(3)) or "-",
                 "   <== decide" if bad else ""))

    if a.verify_pmids:
        import urllib.request
        ids = []
        for p in probes:
            ids += re.findall(r"\b(\d{6,8})\b", p.get("source") or "")
        if ids:
            print("\n  resolving %d PMID(s) against PubMed:" % len(set(ids)))
            url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                   "?db=pubmed&retmode=json&id=%s" % ",".join(sorted(set(ids))))
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
                d = json.loads(urllib.request.urlopen(req, timeout=30).read().decode())
                res = d.get("result", {})
                for k in sorted(set(ids)):
                    e = res.get(k) or {}
                    if not e or e.get("error"):
                        print("    PMID %-9s *** NOT FOUND ***" % k); fail += 1
                    else:
                        print("    PMID %-9s %s. %s %s" % (k, e.get("sortfirstauthor") or "?",
                                                           e.get("source") or "", e.get("pubdate") or ""))
                        print("               %s" % (e.get("title") or "")[:96])
            except Exception as ex:
                #  BV-BRC and NCBI both reject the default urllib User-Agent, so a
                #  failure here is usually the UA or the network, not a bad PMID.
                #  Do not let that silently pass the run as if verified.
                print("    could not reach PubMed (%s) -- PMIDs NOT verified" % ex)
                fail += 1

    if unsourced:
        print("\n  UNSOURCED probe(s), no citation given: %s" % ", ".join(unsourced))
        fail += len(unsourced)

    if a.write_gaps:
        with open(a.write_gaps, "w") as fh:
            fh.write("# Literature proteins vs BV-BRC calls — %s\n\n" % a.module)
            fh.write("Collections are built from deposited strings, so a protein that is\n"
                     "real biology but unannotated in these genomes is invisible to the\n"
                     "build. This records what was probed for and what was found.\n\n")
            fh.write("| protein | deposits | verdict | binds to | source | note |\n")
            fh.write("|---|---|---|---|---|---|\n")
            for r in results:
                fh.write("| %s | %d | %s | %s | %s | %s |\n"
                         % (r["name"], r["deposits"], r["verdict"],
                            ", ".join(r["dests"]) or "—",
                            r.get("source") or "**UNSOURCED**",
                            (r.get("note") or "").replace("|", "/")))
        print("\n  wrote %s" % a.write_gaps)

    if fail:
        print("\n  %d probe(s) need a rule or a written decision." % fail)
        sys.exit(1)
    print("\n  every probed protein is accounted for.")


main()
