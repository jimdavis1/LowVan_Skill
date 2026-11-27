#!/usr/bin/env python3
"""The one place that parses a sequence identifier. Import it; do not re-roll it.

WHY THIS FILE EXISTS. `fig|` IS PART OF THE IDENTIFIER, not a decoration to be
stripped by splitting on `|` and taking field 0. Doing that returns the
literal string `fig` for every sequence in the dump, and the failure is
silent in the worst possible way: nothing errors, nothing is empty, and every
record simply collapses into one bucket.

    fig|687384.4.peg.3   ->  split("|")[0]  ->  "fig"
    fig|687384.9.peg.3   ->  split("|")[0]  ->  "fig"

In `audit_collection_lengths.py` that turned the per-SPECIES length comparison
the tool exists for -- "Compare per species, never per module or genus", and
four paragraphs of Fowl adenovirus fiber-1 argument for why -- into the
module-wide comparison it forbids, in every engine-built module since
`collection_engine.py` was written, because the engine writes every collection
header as `>fig|<genome>.peg.<n>`. Found in Anelloviridae on 27 Nov 2026 only
because a report read `59 aa vs fig median 126`. Correcting it changed real
answers: a Gyrovirus_galga1 VP1 SPLIT finding disappeared (the second length
class was a different species) and genuine within-species outliers appeared in
Lambdatorquevirus and the Anelloviridae minor pool that the module-wide median
had hidden.

A FASTA HEADER IS SPLIT ON WHITESPACE FIRST. The id is the token before the
first space; everything after it is free text that may contain anything,
including `|`. Several call sites split on `|` without splitting on whitespace
at all, which silently folds the annotation into the id.

THE SHAPES THAT OCCUR, all of which these functions accept:

    fig|687384.4.peg.3                       collection_engine.py output
    fig|1128118.16.CDS.2                     ditto, CDS-numbered
    fig|687384.4.peg.2|Eta|2107c579ca11...   ditto with extra fields appended
    172297.3|VP_1                            capsid_length_audit.py input
    10524.1311 long fiber protein            plain dump id plus annotation
    >fig|687384.4.peg.3                      a header with its '>' still on

Anything that does not look like a BV-BRC id is returned whole rather than
mangled, so an unexpected format shows up as itself in a report instead of
quietly becoming a shared key.
"""

__all__ = ["header_id", "genome_of", "taxon_of"]


def header_id(header):
    """The identifier: the token before the first whitespace, '>' removed.

    Split on whitespace BEFORE anything else. `10524.1311 long fiber protein`
    is an id and an annotation, and the annotation is not part of the id.
    """
    if header is None:
        return ""
    t = header.strip()
    if t.startswith(">"):
        t = t[1:]
    return t.split()[0] if t.split() else ""


def genome_of(header):
    """The BV-BRC genome id: 'fig|687384.4.peg.3' -> '687384.4'.

    `fig|` is a prefix OF the id, so it is removed as a prefix, never by
    splitting the id on '|'. What remains is <taxon>.<version>[.peg.<n>],
    and the genome is its first two dot-separated fields.
    """
    t = header_id(header)
    if t.startswith("fig|"):
        t = t[4:]
    #  Formats like '172297.3|VP_1' put the genome first and a feature name
    #  after; formats like 'fig|...|Eta|<md5>' append extra fields. Either
    #  way the genome is in the first '|'-delimited field of what is left.
    t = t.split("|")[0]
    parts = t.split(".")
    if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
        return ".".join(parts[:2])
    return t


def taxon_of(header):
    """The NCBI taxon id: 'fig|687384.4.peg.3' -> '687384'.

    The genome id's prefix IS the taxon id, so the species key is exact and
    free. Never key on the genome NAME: "Human adenovirus 41 Hannover_2022"
    and "Human adenovirus 41 DOM03186" are one species under two strings.
    """
    return genome_of(header).split(".")[0]


if __name__ == "__main__":
    CASES = [
        ("fig|687384.4.peg.3", "687384.4", "687384"),
        ("fig|1128118.16.CDS.2", "1128118.16", "1128118"),
        ("fig|687384.4.peg.2|Etatorquevirus|2107c579ca1148d5", "687384.4", "687384"),
        (">fig|68887.6495.peg.3", "68887.6495", "68887"),
        ("172297.3|VP_1", "172297.3", "172297"),
        ("10524.1311 long fiber protein", "10524.1311", "10524"),
        ("10524.1311", "10524.1311", "10524"),
        ("fig|687384.4.peg.3 ORF3", "687384.4", "687384"),
        ("NC_001427", "NC_001427", "NC_001427"),
        ("", "", ""),
    ]
    bad = 0
    for h, g, t in CASES:
        gg, tt = genome_of(h), taxon_of(h)
        ok = (gg == g and tt == t)
        bad += not ok
        print("%-4s %-50r -> %-14r %-12r" % ("ok" if ok else "FAIL", h, gg, tt))
    #  The regression that motivated the file: distinct genomes must not
    #  share a key.
    keys = {genome_of(x) for x in ("fig|687384.4.peg.3", "fig|687384.9.peg.3",
                                   "fig|1138485.4.peg.3")}
    if len(keys) != 3:
        print("FAIL distinct genomes collapsed to %r" % keys)
        bad += 1
    else:
        print("ok   three distinct genomes give three distinct keys")
    raise SystemExit(1 if bad else 0)
