#!/usr/bin/env python3
"""Normalize `psiblast -out_pssm` ASN.1 text to the LowVan pipeline's format.

psiblast writes an intermediateData block (freqRatios) the pipeline's PSSMs do
not carry, renders a numeric query id as `local id N` rather than
`local str "N"`, and indents differently. This rewrites all three so a
recomputed PSSM is byte-identical to one the pipeline produced.

It also drops the `descr` block. Some pipeline runs emit

    descr {
      title "1 Untitled PSSM 1"
    },

and some do not; it is a label, it holds no scores, and psiblast omits it on a
rebuild. Leaving it in made rebuild_pssms.py report 399 of 808 Hepaciviridae
PSSMs as stale when every score was identical -- the block was the entire
difference, and it correlated perfectly: all 399 flagged carried it, all 409
clean ones did not, and nothing was genuinely out of date. A staleness check
that is wrong half the time is worse than none, because the next real hit gets
ignored with the rest.
"""
import re, sys

def strip_block(t, key):
    i = t.find(key + " {")
    if i < 0:
        return t
    j = t.index("{", i)
    d = 0
    for k in range(j, len(t)):
        if t[k] == "{":
            d += 1
        elif t[k] == "}":
            d -= 1
            if d == 0:
                e = k + 1
                while e < len(t) and t[e] == ",":
                    e += 1
                while e < len(t) and t[e] == "\n":
                    e += 1
                return t[:i] + t[e:]
    return t

def norm(t):
    t = strip_block(t, "intermediateData")
    t = strip_block(t, "descr")
    t = re.sub(r'local id (\d+)', r'local str "\1"', t)
    out, d, inq = [], 0, False
    for line in t.split("\n"):
        s = line.strip()
        if not s:
            continue
        if inq:                       # continuation of a quoted string
            out.append(s)
            if s.count('"') % 2:
                inq = False
            continue
        if s.startswith("}"):
            d -= 1
        out.append("  " * d + s)
        if s.count('"') % 2:
            inq = True
            continue
        d += s.count("{") - s.count("}")
        if s.startswith("}"):
            d += 1
    return "\n".join(out) + "\n"

if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] in ("-h", "--help"):
        sys.exit("usage: norm_pssm.py <raw.pssm>   # normalised PSSM on stdout")
    sys.stdout.write(norm(open(sys.argv[1]).read()))
