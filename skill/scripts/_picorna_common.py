# -*- coding: utf-8 -*-
"""Shared machinery for picornavirus genus triage.

The slot names 1A-3D are genuine homology groups INSIDE a genus and not
across one: enterovirus 2A is a chymotrypsin-like protease, cardiovirus 2A is
not one, aphthovirus 2A is an 18 aa StopGo peptide, and hepatovirus 2A is part
of the VP1-2A (pX) precursor. Each genus therefore gets its own rule list and
its own measured windows; only the matching mechanics are shared.
"""
import re

AMBIGUOUS = re.compile(
    r"^(protease|proteinase|capsid protein|capsid|capsid proteins?|unnamed|"
    r"hypothetical protein|polyprotein|genome polyprotein|viral protein|"
    r"P1|P2|P3|polypeptide|polyprotein precursor|mature peptide|)$", re.I)


def make(rules, prefer=None):
    """-> a triage(s) closure. `prefer` maps a regex to the key that wins when
    a string matches several, for the cases where one protein's NAME contains
    another's (a leader that is itself a protease, say)."""
    comp = [(k, re.compile(p, re.I)) for k, p in rules]
    pref = [(re.compile(p, re.I), k) for p, k in (prefer or [])]

    def triage(s):
        s = (s or "").strip()
        if AMBIGUOUS.match(s):
            return None, "ambiguous" if s else "blank"
        hits = sorted({k for k, rx in comp if rx.search(s)})
        if len(hits) > 1:
            for rx, k in pref:
                if k in hits and rx.search(s):
                    hits = [k]; break
        if len(hits) == 1: return hits[0], "rule"
        if not hits: return None, "unmatched"
        return None, "conflict:" + "+".join(hits)
    return triage


def slots(**extra):
    """The common picornavirus slots. Pass VP0=True for the genera whose
    capsid stays as VP0 rather than cleaving to VP4+VP2."""
    r = []
    #  Ldr, not L: in this project L is the Large protein (the mononegavirus
    #  RdRp, and the gene_symbol the closterovirid RDRP carries).
    if extra.get("L"): r.append(("Ldr", r"\bL\b|\bLab?\b|\bLb\b|\bLpro\b|leader"))
    if extra.get("VP0"):
        r.append(("VP0", r"\bVP0\b|\b1AB\b"))
    else:
        r += [("VP4", r"\b1A\b|\bVP[-\s]?4[ab]?\b|\bP1A\b"),
              ("VP2", r"\b1B\b|\bVP[-\s]?2\b|\bP1B\b")]
    r += [("VP3", r"\b1C\b|\bVP[-\s]?3\b|\bP1C\b"),
          ("VP1", r"\b1D\b|\bVP[-\s]?1\b|\bP1D\b"),
          ("2A",  r"\b2A\b|\bP2[-\s]?A\b|picornain\s*2A"),
          ("2B",  r"\b2B\b|\bP2[-\s]?B\b|viroporin"),
          ("2C",  r"\b2C\b|\bP2[-\s]?C\b|ntpase|atpase|helicase"),
          ("3A",  r"\b3A\b|\bP3[-\s]?A\b"),
          ("VPG", r"\b3B\b|\bvpg\b|genome[-\s]?linked"),
          ("3C",  r"\b3C\b|\bP3[-\s]?C\b|picornain|3Cpro"),
          ("3D",  r"\b3D\b|\bP3[-\s]?D\b|polymerase|\brdrp\b|3Dpol|\bpol\b")]
    return r