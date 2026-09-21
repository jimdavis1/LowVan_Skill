# Orthoflavivirus — when the coordinates are wrong and the annotation string is right

The fourth teaching case. `example-matonaviridae/` covers a source whose
annotation *strings* are wrong. This one covers the opposite: strings that say
exactly the right thing attached to coordinates that point somewhere else.

A 15-product polyprotein has 16 termini, and none of them is marked by a codon
— every one is where a named enzyme cuts. That makes the boundaries checkable,
and it makes a mis-drawn boundary invisible to every rule the kit had.

## What the wording and length rules could not catch

`fig|11060.18332` was filed as 2K carrying the annotation string for 2K. Its
sequence is `RGLRTLILAPTRVVAAEM`, which occurs **1,212 times inside the NS3
collection, always at residue 214 of 618** — a conserved motif in the NS3
helicase, 560 residues away from where 2K lives. It passed every gate: the
string genuinely says 2K, and 18 aa sits inside 2K's declared 15–31 window.

Five West Nile records filed as 2K begin `EKQR`, which is NS4A's own
C-terminal tail — their N-terminal cut is drawn four residues too early. Two
others end `...VAA` **+ `NEMG`** and eat into NS4B.

187 sequences filed as mature C are anchored capsid. C and ancC overlap in
length (C p90 114, ancC p10 110), so no length rule separates them, and the
source calls both "capsid protein".

## The rule, and where it comes from

* **PMID 8445732** — Lin, Amberg, Chambers & Rice, *J Virol* 1993;67(4):2327–2335.
  The NS4A/2K site is cleaved by NS2B-NS3 and sits exactly 23 residues upstream
  of the 2K/NS4B signalase site.
* **PMID 33494395** — Wahaab et al., *Pathogens* 2021;10(2):102. NS2B-NS3
  requires two basic residues (K-R, R-R, R-K, occasionally Q-R) at P2–P1 and a
  small residue (G, S, A) at P1'.
* **PMID 2174669** — Chambers, Hahn, Galler & Rice, *Annu Rev Microbiol*
  1990;44:649–688. The polyprotein cleavage map.

So an NS2B-NS3-generated C-terminus ends in R or K, a signalase-generated one
ends in a small residue, and an NS2B-NS3-generated N-terminus begins with one.
`JUNCTION` in `build_collections.py` is that table, one row per feature.

The check that it is right, rather than merely self-consistent: after
filtering, **2K measures 23 residues at p10, median and p90**. The filter
knows only residue classes at the two ends and is never told a length. It
reproduces Lin 1993 from sequence alone.

## Order matters

Re-binning runs **before** the length and composition filters. A record's
length only means something relative to the feature it actually belongs to: a
3,391-residue polyprotein filed as NS5 is 4× overlong for NS5 and exactly
right for POLY. Running the filters first silently discarded 36 of the 108
anchored-capsid records and both polyprotein records.

## Identifying what a rejected record is

`reassign.py` searches each junction failure against the clean collections and
moves it only where the match covers ≥80% **both ways** at ≥50% identity and
the record satisfies the target feature's junction rule.

Two traps, both of which produced wrong answers before they were fixed:

1. **Filter on coverage before ranking by score.** Every mature peptide is a
   substring of the polyprotein, so POLY hits align at ~100% identity over the
   full query while covering 3–4% of the target. Ranked by bits, POLY wins
   every contest; the coverage gate then discards the *record* instead of the
   hit, and a complete ancC looks like a partial polyprotein.
2. **Search mature peptides and polyproteins separately.** With 14,095
   polyproteins in the target set, they fill any `--max-seqs` window before the
   true assignment appears.

Of 1,422 junction failures, 113 could be placed. The other 1,309 break down as
1,094 whose best hit is their *own* feature — the right protein with broken
ends, which has no bin to move to — 188 matching nothing at ≥80% coverage, and
27 failing the target's own junction rule.

## Result

| | before | after |
|---|---|---|
| sequences | 68,202 | 55,119 |
| termini matching published chemistry | 84–98% | **100% on all 15 features** |
| 2K length p10 / median / p90 | 17 / 23 / 27 | **23 / 23 / 23** |
| PSSMs | 394 | 337 |

Fewer profiles, each trained on sequence that is the protein it claims to be.

## Files

| | |
|---|---|
| `build_collections.py` | binning, the `JUNCTION` table, ambiguity and length filters, re-bin application |
| `build_orthoflavivirus_json.py` | the module block; `min_len`/`max_len` are measured, not hand-drawn |
| `features.json` | annotation string and symbol per feature |
| `reassign.py` | homology identification of junction failures |

NS1′ is deliberately absent from `features.json`. It is real — 533 features,
median 404 aa against NS1's 352, the 52-residue extension of the −1 ribosomal
frameshift, and carried only by the JEV serogroup — but it is non-collinear
and needs `special: transcript_edit` with a hand-curated nucleotide reference
set, not a PSSM. The collection is kept so that work can start from it.
