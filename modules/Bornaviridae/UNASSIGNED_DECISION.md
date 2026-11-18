# The unassigned pool, and why no UNC feature is declared

After `rescue_unassigned.py` (27 of 35 unbound sequences adopted at >=80 %
identity over >=85 % of the query), **8 unique sequences remain unbound** and
two strings reach the member floor of 5 by deposit count:
`unnamed protein product; p24 ORF` and `unnamed protein product; p10 ORF`.

**Neither is a protein this module fails to call. Both are fragment sets, and
all of them come from one submitter's BoDV-1 strain series** (genome ids
`1714621.67` – `1714621.87`), which deposited partial PCR amplicons under
positional "ORF" names:

| string | lengths deposited | full-length protein |
|---|---|---|
| `p24 ORF` | 147, 147, 147, 47, 47 aa | P is 201 aa |
| `p10 ORF` | 63, 63, 63, 26, 26 aa | X is 87 aa |
| `p40 ORF` | 176, 176, 13, 13 aa | N is 370 aa |

The structural test, not the length one: the 147 aa `p24 ORF` entries begin
`LPKDLRKNPSMISDPDQRTGREQLSNDELIKKLVTELAEN...` and the 26 aa `p10 ORF` entries
begin `LPKDLRKNPSMISDPDQRTGREQLSN` — **the same sequence, with no initiator
Met**, so they are internal pieces of one protein deposited twice under two
names. A deposit that starts mid-protein is a fragment whatever it is called.

They are therefore **tabled, not modelled**: declaring a feature from them
would build a profile out of fragments of P and X, which is the `E3_CR1D`
failure recorded in the skill ("a partial that became its own feature is still
a partial"). The full-length P and X of these same BoDV-1 strains are already
in the P and X collections.

`p40 ORF` does not reach the floor and is listed only for completeness.
