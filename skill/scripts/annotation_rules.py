#!/usr/bin/env python3
"""The two annotation-string tests, in one place.

These were copied into three files -- `collection_engine.py`, which bins;
`rescue_unassigned.py`, which must apply the same discard or it hands every
fragment straight back; and `check_annotations.py`, which asserts the rule
cannot delete a published subunit name. Widening the rule in one of them left
the other two disagreeing, which is how Sapovirus `Partial Capsid Protein` came
to be discarded by the binner and then offered to the rescue in the same run.
One definition, three importers.
"""
import re

#  PROJECT RULE (8 October 2026): a record the source annotates as partial or
#  incomplete is DISCARDED -- not binned, not offered to the homology rescue,
#  not parked in <FEAT>.outliers.fasta. The submitter is saying the sequence is
#  not the whole protein, and a profile has no business training on it.
#
#  This is a different test from the length gate and it runs first. The length
#  gate asks whether a sequence is the wrong SIZE; this asks whether the source
#  says it is a FRAGMENT. No length window substitutes: rotavirus VP7 is 326 aa
#  in RVA and 244-251 in RVB/RVG/RVJ, so the window must reach down to ~195,
#  and a TRUNCATED RVA VP7 at 276 then passes it. There were 10,105 of those
#  against 5,461 full-length RVA records.
#
#  TWO SHAPES, both discarded:
#
#    comma suffix   `putative VP1 polyprotein, fragment`   Sapovirus, 762
#                   `putative Capsid protein, N-terminal`  Sapovirus, 371
#    prefix         `truncated VP7`                        Rotavirus, 25
#                   `Partial Capsid Protein`               Sapovirus, 210
#
#  The prefix form is the one that keeps being missed, because the suffix rule
#  reads as complete. Sapovirus `Partial Capsid Protein` has a median length of
#  142 against a true VP1 of 556.
#
#  ANCHORING MATTERS AND IS NOT NEGOTIABLE. The qualifier must be
#  comma-introduced at end-of-string, or one of the three partiality words at
#  the very start. It must NEVER match mid-string: `, N-terminal` and
#  `, C-terminal` appear adjectivally inside 32 legitimate published names for
#  CLEAVED SUBUNITS -- Gn/Gc across five bunyavirus families and both
#  filoviruses, HA1/HA2, HEF1/HEF2, spike S1/S2, fusion F1/F2. Deleting those
#  would silently remove a declared feature from every module carrying one.
#  `check_annotations.guard_partial_vs_subunit` asserts exactly that, against
#  the published vocabulary, using this same object.
PARTIAL = re.compile(
    r",\s*(fragment|partial|incomplete|truncated|[NC]-?terminal)\s*$"
    r"|^\s*(truncated|partial|incomplete)\b", re.I)

#  `putative` is NOT a partiality claim -- it is the submitter saying they are
#  unsure WHICH protein this is. So it is kept, but it has to be EARNED: the
#  string is not taken at face value, the record goes to unassigned, and the
#  homology rescue adopts it only at >=80% identity over >=60% of the query
#  against a named collection. Binning it by name instead would let one
#  uncertain submitter's guess train the profile that then propagates it.
PUTATIVE = re.compile(r"^\s*(putative|probable|possible|predicted)\b", re.I)


def is_partial(anno):
    """True if the source calls this record a fragment, in either shape."""
    return bool(PARTIAL.search(" ".join((anno or "").split())))


def is_putative(anno):
    """True if the submitter is hedging about WHICH protein this is."""
    return bool(PUTATIVE.search(" ".join((anno or "").split())))
