#!/usr/bin/env python3
"""synonyms.tsv -- one schema, one writer, one reader.

This file exists for the same reason `annotation_rules.py` does: the table was
being written in three places and read in two, and the copies drifted.

    collection_engine.py        count / annotation / genus / bins_to
    hiv2/build_collections.py   feature / annotation / count     (3 columns)
    genus/build_cols.py         did not write it at all

The consequences were not loud. `collect_synmap.py` read the 3-column file and
died with `KeyError: 'bins_to'`, which names a dict key rather than the real
problem. The genus modules simply had no file, so the first instinct was to
reconstruct the mapping from the built collections -- which measures something
else ("the sequence reached this collection", post length-window) and reported
EIAV at 61.0% against a true rule-bound 87.0%.

The rules were never missing. Only the emission was. So the emission lives
here, and every builder calls it.

WHAT THE TABLE MEANS. One row per (annotation string, genus, destination),
counting FEATURES not unique sequences. `bins_to` is the destination a RULE
chose, recorded before the length, ambiguity and cleavage filters run -- so a
string whose sequence was later dropped as too short still appears under its
feature. That is deliberate: the table is a record of VOCABULARY, not of
collection membership, and the two differ by exactly the filtered sequences.

RESERVED DESTINATIONS. These are not features and must never be counted as
collapse targets:

    UNASSIGNED     no rule matched
    DISCARDED      the source called it partial/truncated (annotation_rules)
    NOT_MODELLED   a real product the module deliberately does not declare
    SPLICE         a rule matched, but the feature is called by the splice
                   caller rather than a PSSM, so it has no profile

`UNBOUND` is what a collapse report must exclude from its numerator; the other
two reserved names are bound by a rule and belong in the denominator.
"""
import csv
import os

HEADER = ("count", "annotation", "genus", "bins_to")

#  Destinations that are not a feature. A collapse report divides by the number
#  of targets, so counting these inflates the target count and understates the
#  collapse -- that is how Tobamovirus once read 279 strings into 5 for 55.8x
#  when the truth is 267 into 4 for 66.8x.
UNBOUND = ("", "UNASSIGNED")

#  Destinations that are NOT a declared feature and must never appear as a
#  collapse target. They still belong in synonyms.tsv -- the table is a record
#  of vocabulary, and "the source called this partial" is part of that record --
#  but a collapse report counts features, so it skips these.
#
#  There is deliberately no "SPLICE" destination. A spliced feature is a
#  DECLARED feature that happens to be called by get_splice_variant_features.pl
#  rather than a PSSM, so its strings bin to its own key (TAT, REV) and it
#  resolves against the module JSON like any other. Writing "SPLICE" instead
#  loses that and renders on the page as "not modelled", which is wrong: the
#  module models it, just not with a profile.
NON_TARGET = ("UNASSIGNED", "DISCARDED", "NOT_MODELLED")


def write(path, rows):
    """rows: iterable of (count, annotation, genus, bins_to). Sorted by count."""
    rows = [(int(c), a, g or "", k or "UNASSIGNED") for c, a, g, k in rows]
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(HEADER) + "\n")
        for c, a, g, k in sorted(rows, key=lambda x: -x[0]):
            fh.write("%d\t%s\t%s\t%s\n" % (c, a, g, k))
    total = sum(r[0] for r in rows)
    bound = sum(r[0] for r in rows if r[3] not in ("UNASSIGNED", "DISCARDED"))
    return {"rows": len(rows), "features": total, "bound": bound,
            "pct": (100.0 * bound / total) if total else 0.0}


def read(path):
    """Rows as dicts. Refuses a file that is not this schema, and says why."""
    with open(path) as fh:
        first = fh.readline().rstrip("\n").split("\t")
    if tuple(first) != HEADER:
        #  The 3-column file predates this schema. Say so in terms of the fix,
        #  not in terms of a missing dict key.
        hint = ""
        if len(first) == 3:
            hint = ("\n  This looks like the old 3-column schema "
                    "(feature/annotation/count). Re-emit it with the module's "
                    "own builder: `python3 build_collections.py --synonyms-only`"
                    " -- that reuses the same rules and touches no FASTA.")
        raise SystemExit(
            "synonyms.tsv has the wrong schema: %s\n  found:    %s\n  expected: %s%s"
            % (os.path.abspath(path), "/".join(first) or "(empty)",
               "/".join(HEADER), hint))
    return list(csv.DictReader(open(path), delimiter="\t"))
