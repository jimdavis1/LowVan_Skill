#!/usr/bin/env python3
"""Check a module's annotation strings against the LowVan published vocabulary.

Every string a module emits becomes a permanent annotation in BV-BRC, so it has
to match S1 Table (the manuscript's supplemental vocabulary) where a string for
that protein already exists, and follow its conventions where one does not.
Reusing an existing string is always better than coining a near-duplicate: the
whole point of LowVan is to shrink the vocabulary, and "Nucleoprotein" next to
S1 Table's "Nucleocapsid protein" makes it grow.

The governing rule, from the manuscript: *the functional description is recorded
as the annotation string, while the commonly used symbolic name is recorded as
the gene symbol*. The measles polymerase is "RNA-dependent RNA polymerase" with
gene symbol "L" -- not "L protein". S1 Table carries the symbol in its own
column, and the module JSON carries it in `gene_symbol`.

Reports:

  REUSED    string already in S1 Table -- nothing to do, this is the good case
  VARIANT   same string but for case, punctuation or spacing. Always a defect.
            Deliberately narrow: this vocabulary distinguishes proteins by a
            single character (P5 / P6 / P7), so any edit-distance tolerance
            produces only false positives.
  NEW       not in S1 Table -- gets added, and is style-checked here

Also cross-checks `gene_symbol` against S1 Table's Symbol column wherever the
annotation string already exists there. S1 Table namespaces a symbol with a
lineage prefix when one module collapses several subgenera that each own a
same-numbered ORF (`Embeco_NS2a`, `Merbeco_ORF4a`, `G_COV_ORF4a` -- 41 of its
442 symbols), so `Tupa_SH` against a bare `SH` is that convention, not an
error, and is reported separately from a genuine disagreement.

    python3 check_annotations.py --json Taxon_Viral_PSSM.json
"""

import argparse
import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
VOCAB = os.path.join(os.path.dirname(HERE), "assets", "annotation-vocabulary.tsv")

GENERIC = {"protein", "hypothetical protein", "unknown", "orf", "putative protein",
           "uncharacterized protein"}
# qualifiers the manuscript explicitly refuses, to stop the vocabulary growing
BANNED = re.compile(r"\b(incomplete|truncated|partial|fragment|putative|probable)\b", re.I)


def squash(s):
    """Case, punctuation and spacing removed -- what a formatting check compares."""
    return re.sub(r"[^a-z0-9]", "", s.lower())


def style_notes(anno, ftypes, symbols):
    out = []
    if not anno:
        return ["empty"]
    if anno[0].islower():
        out.append("starts lower-case; S1 Table strings are sentence case")
    if anno.lower() in GENERIC:
        out.append('too generic; S1 Table uses "Uncharacterized lineage-specific '
                   '<symbol> protein" for unnamed ORFs')
    if anno.isupper() and len(anno) > 4:
        out.append("all-caps; spell the name out")
    m = BANNED.search(anno)
    if m:
        out.append('drop "%s" -- the manuscript excludes qualifying words from '
                   'annotation strings to stop the vocabulary proliferating' % m.group(0))
    if "mat_peptide" in ftypes and not re.match(r"^(Mature|Signal peptide)", anno):
        out.append('mat_peptide strings begin "Mature ..." or "Signal peptide of <symbol>"')
    if "CDS" in ftypes and anno.startswith("Mature "):
        out.append('"Mature ..." is the mat_peptide convention, not CDS')
    # the family/genus prefix was REMOVED from S1 Table in the current revision
    if re.match(r"^[A-Z][a-z]+(viridae|virinae|virus)\b", anno):
        out.append("drop the taxon prefix; S1 Table removed these "
                   '("Paramyxoviridae C protein" is now just "C protein")')
    # symbol belongs in gene_symbol, not embedded as the whole description
    for sym in symbols:
        if sym and anno.strip().lower() in (sym.lower(), sym.lower() + " protein"):
            out.append('the annotation is just the symbol "%s"; give a functional '
                       "description and keep the symbol in gene_symbol" % sym)
            break
    if len(anno) > 90:
        out.append("very long (%d chars)" % len(anno))
    return out


def symbol_collisions(vocab_rows, ours):
    """Same SYMBOL, different protein, in another taxon.

    check_annotations already asks "our symbol vs S1 Table's symbol for the same
    annotation string". It never asked the reverse, and the reverse is where the
    damage is: a symbol that already means something else elsewhere in the
    vocabulary. Togaviridae shipped `SP` for a 1,248 aa structural polyprotein
    when SP means "signal peptide" -- a mat_peptide of a few dozen residues --
    in ten other taxa, and nothing caught it.

    Not every shared symbol is wrong. `C` is "C protein" in six paramyxoviruses
    and the capsid here, and both are the community's own name, which the
    manuscript explicitly allows. So this reports rather than fails, and says
    what the other meaning is so the call can be made deliberately.
    """
    import collections as _c
    by = _c.defaultdict(set)
    for r in vocab_rows:
        if len(r) >= 3 and r[2]:
            by[r[2]].add((r[0], r[1]))
    out = []
    for taxon, key, sym, anno in ours:
        others = sorted({(t, a) for t, a in by.get(sym, set())
                         if t != taxon and a != anno})
        if others:
            out.append((taxon, key, sym, anno, others))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True, help="module JSON")
    ap.add_argument("--vocab", default=VOCAB, help="S1 Table TSV")
    args = ap.parse_args()

    known, taxa_of, symbol_of = set(), {}, {}
    with open(args.vocab) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            a = (row.get("Annotation") or "").strip()
            if not a:
                continue
            known.add(a)
            taxa_of.setdefault(a, set()).add((row.get("Taxon") or "").strip())
            s = (row.get("Symbol") or "").strip()
            if s:
                symbol_of.setdefault(a, set()).add(s)
    by_squash = {}
    for a in known:
        by_squash.setdefault(squash(a), a)

    #  raw vocabulary rows, for the reverse symbol check below
    with open(args.vocab) as fh:
        vocab_rows = [r for r in csv.reader(fh, delimiter="\t")][1:]

    with open(args.json) as fh:
        mod = json.load(fh)

    seen = {}
    quality_keyed = {}
    for module, block in mod.items():
        for key, ent in block.get("features", {}).items():
            a = (ent.get("anno") or "").strip()
            seen.setdefault(a, []).append(
                (module, key, ent.get("feature_type", ""), ent.get("gene_symbol", "")))
            #  viral_genome_quality.pl keys its %essential hash by the
            #  annotation string, for exactly these features.  Track them
            #  separately so the collision check below matches its rule.
            if ent.get("copy_num") and re.search(
                    r"(CDS|mat_peptide)", str(ent.get("feature_type", ""))):
                quality_keyed.setdefault((module, a), []).append((key, ent))

    #  Two features in one module that share an annotation string collide in
    #  viral_genome_quality.pl: %essential is keyed by the string, so the
    #  second feature read silently overwrites the first one's min_len,
    #  max_len and copy_num, and %anno_count sums both features' calls
    #  against the single surviving copy_num.  The symptoms are a genome
    #  flagged "too many HSPs for: <string>" plus length flags that flip
    #  between "too short" and "too long" from genome to genome, because
    #  Perl randomises hash order and the script runs once per genome.
    #  Tobamovirus shipped REP126 and REP183 both as "Nonstructural
    #  polyprotein" and scored 0% good on 89 genomes for this reason alone.
    dup = {k: v for k, v in quality_keyed.items() if len(v) > 1}
    if dup:
        print("ANNOTATION STRING COLLISION (%d) -- MUST FIX BEFORE INSTALL:" % len(dup))
        for (module, a), lst in sorted(dup.items()):
            print("  %s: %d features share %r" % (module, len(lst), a))
            for key, ent in sorted(lst):
                print("      %-14s copy_num=%-3s %s-%s"
                      % (key, ent.get("copy_num"),
                         ent.get("min_len"), ent.get("max_len")))
        print("  Give each feature its own string.  A readthrough or nested",
              "pair is still two proteins:\n  Togaviridae names them",
              "'Nonstructural polyprotein' (nsP1234) and",
              "'Nonstructural polyprotein P123' (nsP123) for this reason.\n")

    reused, variant, new = [], [], []
    for a, uses in sorted(seen.items()):
        if a in known:
            reused.append((a, uses))
        elif squash(a) in by_squash:
            variant.append((a, uses, by_squash[squash(a)]))
        else:
            new.append((a, uses))

    print("%d distinct annotation strings in %s\n"
          % (len(seen), os.path.basename(args.json)))

    print("REUSED from S1 Table (%d) -- these keep the vocabulary small:" % len(reused))
    for a, uses in reused:
        print("  %-52s %s" % (a, ", ".join(sorted({u[0] for u in uses}))[:56]))

    # gene_symbol agreement, only where S1 Table has a symbol for that string.
    # A "<Lineage>_<SYM>" prefix is S1 Table's own way of disambiguating a
    # symbol shared by several subgenera inside one module, so strip it before
    # comparing and report those separately.
    mism, nspaced = [], []
    for a, uses in reused:
        want = symbol_of.get(a)
        if not want:
            continue
        wl = {w.lower() for w in want}
        for module, key, _ft, sym in uses:
            if not sym or sym.lower() in wl:
                continue
            base = sym.split("_", 1)[1] if "_" in sym else None
            if base and base.lower() in wl:
                nspaced.append((module, key, a, sym, base))
            else:
                mism.append((module, key, a, sym, sorted(want)))
    if nspaced:
        print("\nGENE SYMBOL namespaced (%d) -- matches S1 Table practice, no action:" % len(nspaced))
        for module, key, a, sym, base in nspaced:
            print("  %-20s %-14s %-38s %-14s (base %s)"
                  % (module, key, a[:38], sym, base))
    if mism:
        print("\nGENE SYMBOL differs from S1 Table (%d) -- decide which is right:" % len(mism))
        for module, key, a, sym, want in mism:
            print("  %-20s %-14s %-38s ours %-8s S1 Table %s"
                  % (module, key, a[:38], sym, "/".join(want)))

    if variant:
        print("\nVARIANT of an existing string (%d) -- use the S1 Table spelling:" % len(variant))
        for a, uses, k in variant:
            print("  %-52s\n      S1 Table has:  %-42s (%s)"
                  % (a, k, ", ".join(sorted(taxa_of[k]))[:40]))

    print("\nNEW to the vocabulary (%d) -- add these rows to\n"
          "  assets/annotation-vocabulary.tsv:" % len(new))
    for a, uses in new:
        ft = {u[2] for u in uses}
        syms = {u[3] for u in uses if u[3]}
        print("  %-52s %-12s %s" % (a, "/".join(sorted(ft)),
                                    ", ".join(sorted({u[0] for u in uses}))[:34]))
        for note in style_notes(a, ft, syms):
            print("        style: %s" % note)

    #  The reverse of the symbol check above: not "is our symbol what S1 Table
    #  uses for this string", but "does this symbol already mean something else".
    ours = [(m, k, sym, a) for a, lst in seen.items()
            for (m, k, _ft, sym) in lst if sym]
    coll = symbol_collisions(vocab_rows, ours)
    if coll:
        print("\nSYMBOL COLLISION (%d) -- this symbol already names a different\n"
              "  protein in another taxon. Sometimes correct (a community name the\n"
              "  taxon genuinely owns), sometimes not. Decide it deliberately:" % len(coll))
        for taxon, key, sym, anno, others in coll:
            print("  %-10s %-10s ours: %s" % (sym, key, anno))
            for t, a in others[:6]:
                print("  %-10s %-10s ALSO %-22s = %s" % ("", "", t, a))
            if len(others) > 6:
                print("  %-10s %-10s ... and %d more" % ("", "", len(others) - 6))

    print("\n  columns: Taxon | Annotation | Symbol | Feature Type |")
    print("                         Segment | Used for Genome Quality | PubMed IDs*")
    return 1 if (variant or mism or dup) else 0


if __name__ == "__main__":
    sys.exit(main())
