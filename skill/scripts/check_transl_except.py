#!/usr/bin/env python3
"""Check that a feature's transl_except window is unambiguous.

`transl_except` declares that within `offset` residues of a named terminus, a
residue translating to one of `replace_residue` is really `insert_residue`.
The annotator substitutes only when EXACTLY ONE residue in that window
qualifies; when two do it substitutes nothing, because choosing between them
would be guessing.

That rule is safe but it is not free: every member on which it declines is a
genome that loses the feature. So the window has to be narrow enough that
declining is rare, and the only way to know is to measure it against the
collection the profiles are built from.

This is the gate that makes the declaration reusable. Run it whenever a
feature gains or changes a `transl_except`, and read the AMBIGUOUS count as
the fraction of genomes on which the substitution will refuse to fire.

    python3 check_transl_except.py --workdir W --module M --features F.json

It also reports what the window would cost at one residue wider on each side,
which is the comparison that tells you whether a range was chosen from the
data or from habit. On Deltavirus L-HDAg the answer is stark: offsets 19-20
leave 0 of 1,078 members ambiguous, 18-21 leaves 8, and 15-20 leaves 24.
"""
import argparse
import json
import os
import sys


def read_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def window_indices(n, lo, hi, frm):
    for k in range(lo, hi + 1):
        i = (n - k) if frm == "c-term" else (k - 1)
        if 0 <= i < n:
            yield i


def count(seqs, replace, lo, hi, frm):
    none = one = many = 0
    examples = []
    for sid, s in seqs:
        q = [i for i in window_indices(len(s), lo, hi, frm) if s[i] in replace]
        if not q:
            none += 1
        elif len(q) == 1:
            one += 1
        else:
            many += 1
            if len(examples) < 5:
                examples.append((sid, [len(s) - i for i in q]))
    return none, one, many, examples


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--module", required=True)
    ap.add_argument("--features", required=True)
    ap.add_argument("--collections")
    a = ap.parse_args()

    cdir = a.collections or os.path.join(a.workdir, "collections", a.module)
    if not os.path.isdir(cdir):
        cdir = os.path.join(a.workdir, a.module, "collections")
    feats = json.load(open(a.features))

    declared = {k: v["transl_except"] for k, v in feats.items()
                if isinstance(v, dict) and v.get("transl_except")}
    if not declared:
        print("no feature declares transl_except -- nothing to check")
        return 0

    print("# transl_except window check -- %s\n" % a.module)
    bad = 0
    for k, te in sorted(declared.items()):
        path = os.path.join(cdir, k + ".fasta")
        if not os.path.exists(path):
            print("  %-10s NO COLLECTION at %s" % (k, path))
            bad += 1
            continue
        seqs = list(read_fasta(path))
        replace = set(te.get("replace_residue") or [])
        ins = te.get("insert_residue")
        off = te.get("offset")
        lo, hi = (off, off) if isinstance(off, int) else (off[0], off[-1])
        frm = (te.get("from") or "n-term").lower()

        #  The collection is normalised: members already carry insert_residue
        #  at the site. Count it as a qualifying position too, or a normalised
        #  collection reports every member as "none" and the check says
        #  nothing at all.
        probe = replace | ({ins} if ins else set())
        none, one, many, ex = count(seqs, probe, lo, hi, frm)
        n = len(seqs)
        print("  %s  (%d members, offsets %d-%d from the %s, %s -> %s)"
              % (k, n, lo, hi, frm, "/".join(sorted(replace)), ins))
        print("     exactly one in the window : %d" % one)
        print("     none in the window        : %d" % none)
        print("     AMBIGUOUS (two or more)   : %d%s"
              % (many, "   <-- these lose the feature" if many else ""))
        for sid, offs in ex:
            print("        %s at offsets %s" % (sid, offs))
        if none:
            print("     note: %d members hold neither the replaced nor the "
                  "inserted residue in the window." % none)
        for d in (1, 2):
            _, _, m2, _ = count(seqs, probe, max(1, lo - d), hi + d, frm)
            print("     widening by %d on each side would make %d ambiguous"
                  % (d, m2))
        if many:
            bad += 1
        print()

    if bad:
        print("**%d feature(s) have an ambiguous window.** Narrow the offset "
              "range, or accept that the substitution will decline on those "
              "members." % bad)
        return 1
    print("**Unresolved: 0.** Every member resolves to a single position.")
    return 0


sys.exit(main())
