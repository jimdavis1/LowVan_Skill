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
    python3 check_dlits.py --json Taxon_Viral_PSSM.json --verify
    python3 check_dlits.py --json Taxon_Viral_PSSM.json --offline   # no PubMed

**--verify cannot confirm that a paper supports a claim.** It catches the gross
mismatch -- a citation that is about something else entirely -- and nothing
finer. A paper that is on-topic and still does not establish what the feature
asserts passes it. Only a curator reading the paper settles that.

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
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
FLAG = "PMID_claude_generated"
#  What the proposed paper SUPPORTS. A citation for what a protein does and one
#  for where its boundaries are are different claims with different failure
#  modes: a wrong function mislabels a feature, a wrong coordinate silently
#  moves a cut site and nothing downstream catches it. So the flag is a mapping
#  from category to ids, not a bare list.
CATEGORIES = ("function", "coordinates")

def malformed(ent):
    """Why this feature's flag field cannot be read, or None if it can.

    The field is {category: [pmid, ...]}. The inverted shape {pmid: category}
    looks just as reasonable when you are writing one by hand, and it fails
    silently and absurdly: `for ids in v.values()` then iterates the STRING
    "function" and extends the id list one character at a time, so the module
    reports citations named 'f', 'u', 'n', 'c'. Written out, that is obviously
    broken; inside a loop it just produces 12 NOT A PMID lines that look like
    somebody pasted DOIs. Check the shape instead of trusting it.
    """
    v = ent.get(FLAG)
    if v is None or isinstance(v, list): return None
    if not isinstance(v, dict): return "expected an object, found %s" % type(v).__name__
    for cat, ids in v.items():
        if isinstance(ids, str):
            return ("%r maps to the string %r -- the field is "
                    "{category: [pmid, ...]}, not {pmid: category}" % (cat, ids))
        if not isinstance(ids, list):
            return "%r maps to %s, expected a list of ids" % (cat, type(ids).__name__)
    return None


def flagged_ids(ent):
    """Every proposed id, whatever shape the field is in."""
    v = ent.get(FLAG)
    if not v or malformed(ent): return []
    if isinstance(v, list): return list(v)          # pre-2 October form
    out = []
    for cat, ids in v.items(): out.extend(ids or [])
    return out

def uncategorised(ent):
    """Proposed ids filed under neither 'function' nor 'coordinates'."""
    v = ent.get(FLAG)
    if not v or malformed(ent): return []
    if isinstance(v, list): return list(v)
    out = []
    for cat, ids in v.items():
        if cat not in CATEGORIES: out.extend(ids or [])
    return out



#  Words that appear in so many virology titles that sharing one says nothing.
STOP = set("""a an and are as at be by for from in into of on or the to with
virus viral viruses protein proteins gene genes genome genomes genomic rna dna
sequence sequences analysis study studies new novel characterization
identification molecular structure structural function functional role roles
complete type types strain strains isolate isolates expression evidence
during via using between among their its two three
""".split())

WORD = re.compile(r"[a-z0-9]+")

def content(text):
    """Lowercased content words of a string."""
    return {w for w in WORD.findall((text or "").lower())
            if len(w) > 2 and w not in STOP}


def feature_terms(module, key, ent):
    """The words that would have to appear in a relevant paper's title."""
    t = content(ent.get("anno", "")) | content(ent.get("gene_symbol", "")) | content(key)
    #  a bare positional key like "3D" or "2A" survives WORD but means nothing
    #  on its own in a title, so it is not evidence of relevance
    return {w for w in t if not re.fullmatch(r"\d+[a-z]?", w)}


RANK_SUFFIX = ("virinae", "viridae", "virales", "virus", "viroidae")
RANK_PREFIX = ("alpha", "beta", "gamma", "delta", "epsilon", "ortho", "para",
               "novi", "pseudo", "meta")

def name_stems(module):
    """Searchable stems for a taxon name.

    A module is named at whatever rank the partition sits at, and papers are
    not. "Alpharhabdovirinae"[:6] is "alphar", which does not occur in
    "ICTV Virus Taxonomy Profile: Rhabdoviridae" -- so a naive prefix test
    called 27 correct Rhabdoviridae citations fabrications. Strip the rank
    suffix and then the rank prefix, keeping both forms:
    Alpharhabdovirinae -> alpharhabdo -> rhabdo, which matches Rhabdoviridae,
    rhabdovirus and rhabdoviral alike.
    """
    out, n = set(), module.lower()
    for suf in RANK_SUFFIX:
        if n.endswith(suf): n = n[:-len(suf)]; break
    if len(n) >= 4: out.add(n)
    for pre in RANK_PREFIX:
        if n.startswith(pre) and len(n) - len(pre) >= 4:
            out.add(n[len(pre):])
    return {x for x in out if len(x) >= 4}


def module_vocab(module, block):
    """Every word a paper about this module could reasonably use.

    The module name alone is not enough -- nothing in a rubella paper says
    "Matonaviridae" and nothing in a hepatitis C paper says "Hepaciviridae".
    But the module already carries a list of the organisms it covers, in
    close_genomes, and those are exactly the names the literature uses.
    """
    v = set()
    for d in (block.get("close_genomes") or {}).values():
        v |= content(d.get("genome_name", ""))
    for ent in (block.get("features") or {}).values():
        v |= content(ent.get("anno", ""))
    return v


def relevant(text, module, block, extra, vocab):
    """Does this paper look like it is about this taxon at all?"""
    tl = (text or "").lower()
    if any(st in tl for st in name_stems(module)): return True
    if any(e.lower() in tl for e in extra if len(e) > 3): return True
    return bool(vocab & content(text))


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
                au = v.get("sortfirstauthor") or ""
                out[k] = (title, (v.get("pubdate") or "")[:4], au,
                          (v.get("source") or ""))
        time.sleep(0.4)
    return out


def abstracts(ids):
    """{pmid: abstract text}. Title alone is too thin to judge relevance --
    Strauss 1994 "The alphaviruses: gene expression, replication, and
    evolution" shares no word with "Sindbis virus", which is what the module
    is named after, and was called a fabrication on the title alone."""
    out, ids = {}, list(ids)
    for i in range(0, len(ids), 100):
        chunk = ids[i:i + 100]
        url = "%s?%s" % (EFETCH, urllib.parse.urlencode(
            {"db": "pubmed", "id": ",".join(chunk), "rettype": "abstract",
             "retmode": "text"}))
        try:
            with urllib.request.urlopen(url, timeout=120) as fh:
                txt = fh.read().decode("utf8", "replace")
        except Exception:
            continue
        #  records are separated by blank lines and carry "PMID: nnnn" at the
        #  end; split on that rather than guessing the record boundary
        cur = []
        for block in txt.split("\n\n"):
            cur.append(block)
            m = re.search(r"PMID:\s*(\d+)", block)
            if m:
                out[m.group(1)] = "\n".join(cur)
                cur = []
        time.sleep(0.4)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True, nargs="+",
                    help="one or more module JSONs. Pass them ALL when "
                         "verifying: a citation cleared in one module is "
                         "accepted in its siblings, which is what lets a "
                         "family-level paper support a genus-level module.")
    ap.add_argument("--offline", action="store_true",
                    help="skip PubMed; structural checks only")
    ap.add_argument("--verify", action="store_true",
                    help="resolve every citation and flag the ones whose paper "
                         "has nothing to do with the feature citing it")
    ap.add_argument("--taxon-terms", default="",
                    help="comma-separated extra names for this taxon that a "
                         "relevant title might use instead of the module name "
                         "(e.g. poliovirus,rhinovirus,coxsackie for Enterovirus)")
    args = ap.parse_args()

    mod = {}
    for path in args.json:
        with open(path) as fh:
            for k, v in json.load(fh).items():
                mod[k] = v

    used, flagged = defaultdict(list), defaultdict(list)
    disagree, uncat, badshape = [], [], []
    ctx = {}                      # pmid -> [(module, key, ent), ...]
    for module, block in mod.items():
        for key, ent in block.get("features", {}).items():
            pm = [str(x) for x in (ent.get("PMID") or [])]
            #  flagged_ids() handles both shapes. main() used to do
            #  `for p in ent.get(FLAG)`, which on the mapping form iterates the
            #  CATEGORY NAMES -- so every migrated module reported 0 proposals
            #  across 0 features and every check below silently had nothing to
            #  look at. 72 features were in that state.
            fl = [str(x) for x in flagged_ids(ent)]
            u = [str(x) for x in uncategorised(ent)]
            bad = malformed(ent)
            if bad: badshape.append((module, key, bad))
            if u: uncat.append((module, key, u))
            for p in pm:
                used[p].append("%s/%s" % (module, key))
            for p in fl:
                flagged[p].append("%s/%s" % (module, key))
            for p in set(pm) | set(fl):
                ctx.setdefault(p, []).append((module, key, ent))
            #  PMID and PMID_claude_generated are DISJOINT, not nested.
            #  PMID holds citations a curator supplied and read -- real DLITs.
            #  PMID_claude_generated holds model proposals, which are not DLITs
            #  however real the paper turns out to be. Listing a proposal in
            #  both duplicates it and blurs the one distinction the second field
            #  exists to draw: a consumer reading PMID would take proposals as
            #  curated. So the error is an id appearing in BOTH, not one
            #  appearing only in the flag.
            both = [p for p in fl if p in pm]
            if both:
                disagree.append((module, key, both))

    allids = dict(used)
    for p, w in flagged.items(): allids.setdefault(p, []).extend(w)
    numeric = {p for p in allids if re.fullmatch(r"\d+", p)}
    nonnum = sorted(p for p in allids if p not in numeric)

    nfeat = sum(1 for m, b in mod.items() for k, e in b.get("features", {}).items()
                if e.get("PMID") or e.get(FLAG))
    print("%d curated DLIT(s) and %d model proposal(s) across %d features\n"
          % (len(used), len(flagged), nfeat))

    problems = 0
    if nonnum:
        problems += len(nonnum)
        print("NOT A PMID (%d) -- the PMID field takes bare PubMed ids only:" % len(nonnum))
        for p in nonnum:
            print("  %-34s used by %s" % (p, ", ".join(allids[p])[:60]))
        print("  Resolve a DOI at pubmed.ncbi.nlm.nih.gov and store the numeric id.\n")

    if badshape:
        problems += len(badshape)
        print("MALFORMED %s (%d) -- cannot be read, so nothing below sees these:"
              % (FLAG, len(badshape)))
        for module, key, why in badshape:
            print("  %-20s %-14s %s" % (module, key, why))
        print()

    if uncat:
        problems += len(uncat)
        print("UNCATEGORISED (%d) -- every proposal must say what it supports,"
              " 'function' or 'coordinates':" % len(uncat))
        for module, key, ids in uncat:
            print("  %-20s %-14s %s" % (module, key, ", ".join(ids)))
        print()

    if disagree:
        problems += len(disagree)
        print("CITATION IN BOTH FIELDS (%d) -- they must be disjoint:" % len(disagree))
        for module, key, extra in disagree:
            print("  %-20s %-14s in PMID and in %s: %s"
                  % (module, key, FLAG, ", ".join(extra)))
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
                print("  %-10s %-16s used by %s" % (p, tag, ", ".join(allids[p])[:54]))
            print("  A model-proposed id that does not resolve is a fabricated citation.")
            print("  Remove it; do not keep it with the flag.\n")

    if args.verify and resolved:
        extra = [t.strip() for t in args.taxon_terms.split(",") if t.strip()]
        abst = abstracts(sorted(set(resolved)))
        vocab = {}
        for module, block in mod.items():
            vocab[module] = module_vocab(module, block)
        #  A citation is judged per (pmid, module) pair, and a module is named
        #  at whatever rank its partition sits at. So the ICTV Rhabdoviridae
        #  profile is accepted for Alpharhabdovirinae -- the name stems to
        #  "rhabdo" -- and rejected for Dichorhavirus and Merhavirus, which are
        #  rhabdovirid genera whose names do not contain the family stem. The
        #  paper is the same paper. If an id is relevant to ANY module in this
        #  run it is not a fabrication, so clear it everywhere; that is the
        #  whole reason --json takes several files.
        cleared = set()
        suspect, weak, ok = [], [], 0
        for p in sorted(numeric):
            if p not in resolved: continue
            title, year, au, src = resolved[p]
            body = title + "\n" + abst.get(p, "")
            for module, key, ent in ctx.get(p, []):
                terms = feature_terms(module, key, ent)
                if terms & content(body):
                    ok += 1; cleared.add(p)
                elif relevant(body, module, mod[module], extra, vocab.get(module, set())):
                    weak.append((p, module, key, ent.get("anno", ""), title))
                    cleared.add(p)
                else:
                    suspect.append((p, module, key, ent.get("anno", ""), title, au, year,
                                    p in abst))
        rescued = [r for r in suspect if r[0] in cleared]
        suspect = [r for r in suspect if r[0] not in cleared]
        for r in rescued:
            weak.append((r[0], r[1], r[2], r[3], r[4]))
        if rescued:
            print("  %d pair(s) cleared by a sibling module: the same paper is "
                  "relevant to\n  another module in this run, so it is not a "
                  "fabrication.\n" % len(rescued))
        print("VERIFY: %d of %d citation-feature pair(s) name their feature; "
              "%d are about\n        the taxon but not visibly the feature; "
              "%d are about neither.\n"
              % (ok, ok + len(weak) + len(suspect), len(weak), len(suspect)))
        if suspect:
            problems += len(suspect)
            print("SUSPECT (%d) -- resolves, but to a paper about something else:"
                  % len(suspect))
            for p, module, key, anno, title, au, year, hadabs in suspect:
                print("  %-9s %s/%s%s" % (p, module, key,
                      "" if hadabs else "   [no abstract -- judged on title only]"))
                print("            cited for : %s" % anno[:64])
                print("            resolves  : %s %s, %s" % (au, year, title[:58]))
            print("  A recalled PMID is not evidence. Look each of these up and "
                  "either\n  replace it with the right id or remove it.\n")
        if weak:
            print("NOT FEATURE-SPECIFIC (%d) -- right taxon, but the paper never "
                  "names this\n  feature. Normal for a genome or taxonomy paper "
                  "cited for coordinates;\n  worth a look where it is cited for "
                  "'function':" % len(weak))
            for p, module, key, anno, title in weak[:40]:
                print("  %-9s %-28s %-26s %s"
                      % (p, "%s/%s" % (module, key), anno[:26], title[:34]))
            if len(weak) > 40: print("  ... and %d more" % (len(weak) - 40))
            print()

    print("CITATIONS (%d model-proposed, %d human-curated):\n"
          % (len(flagged), len(numeric) - len(flagged)))
    print("  %-10s %-16s %-5s %s" % ("PMID", "provenance", "year", "title"))
    for p in sorted(numeric, key=lambda x: (x not in flagged, x)):
        tag = "MODEL-PROPOSED" if p in flagged else "human"
        if resolved and p in resolved:
            title, year = resolved[p][0], resolved[p][1]
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
