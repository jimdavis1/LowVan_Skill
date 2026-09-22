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

## NS1′, and the limit of sequence-based identification

`NS1P` is built, as a `transcript_edit` feature with 104 nucleotide
references. `build_ns1p_refs.py` constructs them, and it is worth reading for
one reason: the construction is driven by the **known protein**, not by the
slippery-heptamer motif.

Motif-driven construction fails two ways this kit's own notes warn about —
testing the heptamer at the wrong codon phase rejects every genome and looks
like a finding, and a window anchored on the called feature misses the site
whenever an N-run truncated the call. Instead the script locates NS1's
N-terminus in all three frames, walks the 0-frame until the genome stops
agreeing with the protein, tries every duplication in a window, and keeps only
one that reproduces the protein byte-identically. All 104 donors round-trip
exactly, and no false positive appears outside the serogroup — 0 of 72
dengue, Zika and yellow fever genomes get a call or even a `partial_cds`.

It collects **every** valid duplication rather than the first. Taking the first
gave 109 references, five of which put the frameshift one residue before the
stop: those reproduce their own genome's protein, because nearly the whole
construct is copied verbatim, while placing the slip nowhere near the real
site. They were exactly the five with no heptamer at the chosen position.
Preferring motif-supported solutions drops them and leaves a set with no spread
in either diagnostic — 104 of 104 with the heptamer at phase 2 of NS2A codon
6, and 44 residues after the slip in every one. (44 + the 8 in-frame NS2A
codons before the slip is the 52-residue extension in the literature.)

### The limitation, which is intrinsic

Every JEV-serogroup genome carries that heptamer, intact and in frame — 70 of
79 readable JEV genomes do. What decides whether NS1′ is actually made is a
pseudoknot about 60 bases downstream, and the SA14-14-2 vaccine lineage carries
a single substitution there (`A` at NS2A nucleotide 66) that abolishes the
frameshift. A sequence match sees a near-perfect reference and fills the gap.

| NS2A nt66 | annotator | n | |
|---|---|---|---|
| `G` wild type | NS1′ called | 11 | correct |
| `G` wild type | `partial_cds` | 7 | no close reference; honest |
| `A` vaccine allele | NS1′ called | **16** | **false positive** |
| `A` | withheld | 0 | |

Wild-type JEV, West Nile, Usutu, Murray Valley and Cacipacore are all called
correctly. The error is confined to the attenuated lineage and to six other
JEV genomes carrying the same allele, whose status — unlabelled vaccine
derivative or genuine `A66` field isolate — cannot be settled from sequence.

Checking the heptamer before accepting a correction, the obvious guard, catches
none of them. **This is the limit of sequence-based identification for this
feature**, and it is why the feature ships with the limitation documented
rather than with a threshold tuned around it. Discriminating the vaccine
lineage would mean reading the pseudoknot, which is a change to
`get_transcript_edited_features.pl`, not to the reference set.

### Two things to know before measuring it

`run_gto_eval.py` does **not** run the transcript-edit pass, for any module, so
no quality figure on the coverage audit includes NS1′. Measure it by invoking
`get_transcript_edited_features.pl` on the annotated GTOs.

The reference set under-covers JEV: 97 of 104 references are West Nile and 4
are JEV, because that is the ratio in which BV-BRC annotates the feature. A
denser JEV set is the one improvement that would help.

### Provenance of the original exclusion

NS1′ is deliberately absent from `features.json`. It is real — 533 features,
median 404 aa against NS1's 352, the 52-residue extension of the −1 ribosomal
frameshift, and carried only by the JEV serogroup — but it is non-collinear
and needs `special: transcript_edit` with a hand-curated nucleotide reference
set, not a PSSM. The collection is kept so that work can start from it.
