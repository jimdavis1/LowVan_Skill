# Amnoonviridae and Roniviridae — build notes

Built 22 November 2026 as **three modules**, one per genome organisation.
Working trees `~/lowvan_scratch/amnoon` and `~/lowvan_scratch/roni`; test repos
`VA_amnoon` and `VA_roni`.

| module | family | genome | features | PSSMs |
|---|---|---|---|---|
| `Tilapinevirus` | Amnoonviridae | 10 segments, −ssRNA, orthomyxo-like | 11 | 13 |
| `Okavirus` | Roniviridae | 26.2–26.8 kb +ssRNA, 5–6 ORFs | 9 | 12 |
| `Nimanivirus` | Roniviridae | 28.9–29.1 kb +ssRNA, 4 ORFs | 4 | 3 |

Every protein declared here is backed by a paper whose PMID was resolved
against PubMed and whose title was read against the claim it supports. All
**19 PMIDs** resolve; none was written from memory.

---

## 1. The partition, and why it is three modules

### Amnoonviridae → one module

ICTV ratifies one genus, *Tilapinevirus*, one species, *Tilapinevirus
tilapiae* (tilapia lake virus, TiLV). 607 of 619 BV-BRC records are that
species. The other 12 are discussed in §7.

### Roniviridae → two modules, and the split is measured not asserted

BV-BRC and ICTV both place two genera here: *Okavirus* (gill-associated virus
GAV, yellow head virus YHV-1, YHV-8) and *Nimanivirus* (*Macrobrachium
rosenbergii* Golda virus, MrGV).

The ICTV profile states that MrGV **lacks the triple-membrane-spanning p25
protein** at the N-terminus of the okavirus glycoprotein precursor, and
**lacks ORFX**. A module declares a protein *set*, so that alone separates
them.

That is a claim about alignments, and it was re-measured. MrGV ORF3 (1503 aa)
blastp'd against YHV pp3 (1666 aa) returns **no meaningful alignment at all**:
best HSP **20.4 bits, 23.7% identity over 38 residues**, which is noise. The
two precursors are not alignable, so no okavirus cleavage site can be
transferred to MrGV and no PSSM could serve both.

Routing confirms it independently (§2).

---

## 2. Routing — measured on the shipped reference set

16 rep contigs in total (8 + 5 + 3), far inside the 25-per-module ceiling.
Measured at the **annotator's own BLASTn settings** (`-evalue 0.5 -reward 2
-penalty -3 -word_size 11 -soft_masking false`, 150-bit floor), with all three
modules' references in **one** database so that misrouting is visible — a
module measured alone has nowhere to send a genome and absorbs it.

| family | contigs | routed | misrouted |
|---|---|---|---|
| Amnoonviridae (in scope) | 546 | **546 (100%)** | 0 |
| Roniviridae | 234 | **234 (100%)** | 0 |

The Roniviridae split came out 190 Okavirus / 44 Nimanivirus, matching the
metadata genus counts exactly.

The 12 Amnoonviridae contigs that do **not** route are precisely the 12
out-of-scope records of §7 — the correct outcome: a group this module declines
rather than silently absorbs.

---

## 3. Tilapinevirus: binning without usable annotation strings

**363 of 399 distinct deposits are annotated `hypothetical protein`**, and the
only named protein in the whole dump is PB1. A rules table over the strings can
separate exactly one protein of eleven.

The ten segments are non-homologous to one another, so **homology** separates
them perfectly. The reference set is Til-4-2011 — the Bacharach 2016 isolate —
whose ten segments are deposited as **KU751814…KU751823 in segment order** at
exactly the published lengths (1641/1471/1371/1250/1099/1044/777/657/548/465
nt). TiLV-Israel-HK (OQ437054…OQ437063) reproduces the same ten protein
lengths independently, so the mapping is measured twice:

```
S1 519   S2 457   S3 419   S4 356   S5 343
S6 317   S7 195   S8 174   S9 118   S10 113   aa
```

**Validated against the submitters' own labels.** 429 of 558 contigs state
their segment in the FASTA description; blastn against the reference set agreed
with the stated segment in **429 of 429 (100%)**, and additionally assigned 117
contigs carrying no label at all.

### Segment 3 is the diagnostic RT-PCR target

151 of 546 in-scope contigs are segment 3, and **71 of its protein deposits are
amplicon translations**. They were excluded by the length window, and the
exclusion was then checked structurally rather than assumed: **0 of 71 begin at
reference residue 1**, and **58 of 71 begin at reference residues 100–124** —
one primer site — at 98.7% median identity. Internal fragments, not a short
form.

---

## 4. Two proteins the databases do not hold

### S9-F3, the eleventh TiLV protein — zero deposits

PMID 41513822 (Pankaew et al., *npj Viruses* 2026) reports an eleventh protein
translated from a **105-codon ORF in frame 3 of segment 9**, overlapping the
116-codon primary S9 ORF, from an AUG in a strong Kozak context. It was
confirmed by mass spectrometry — **4 unique peptides, 49% sequence coverage** —
and carries a CRM1-dependent nuclear export signal. The authors attribute it to
**leaky ribosomal scanning, not splicing**.

BV-BRC holds **none of it**. It was extracted from each genome's own segment 9:

- **32 extractions, every one exactly 105 aa**
- 82.9–100% pairwise identity, median 96.2%, 20 distinct sequences
- **all 32 share the identical N-terminus `MDRKYRFCVSNLDR`**

The two isolates PMID 41513822 names as exceptions (OM469318, premature stop at
codon 50; MG873154, no S9-F3 AUG) are not in the segment-9 anchor set and were
therefore **not tested either way** — stated here because it would be easy, and
wrong, to present their absence as confirmation.

### px / ORFX — zero deposits, and a published coordinate corrected

The ICTV profile states ORFX "commences 3 nucleotides downstream of the pp3
initiation codon". **A +3 nt offset is the same reading frame as pp3 and cannot
encode a different protein.** Measured directly in all three okavirus reference
genomes, the start is **+5 nt, in the +2 frame**:

```
YHV  NC_043505   pp3 starts nt 20985   ORFX nt 20990 (+5)  86 aa
GAV  NC_010306   pp3 starts nt 20675   ORFX nt 20680 (+5)  86 aa
YHV8 KX947267    pp3 starts nt 20979   ORFX nt 20984 (+5)  86 aa
```

15 extractions across the genus, **every one exactly 86 aa** (14 at +5, 1 at
+11) — matching the 86 aa the ICTV profile itself gives for px.

---

## 5. The okavirus glycoprotein cleavages — published, then re-measured

PMID 12655087 placed the two pp3 cuts between **Ala228/Thr229** and
**Ala1127/Leu1128** of the 1666 aa YHV precursor. Checked against every
full-length pp3 in the dump, both verify exactly in all five 1666 aa isolates:

```
site 1   ...LCFFITPAFA | TILSGI...     A-F-A, cut after A228
site 2   ...SIISPASA   | LAPRQ...      A-S-A, cut after A1127
```

Both are the **[Ala-X-Ala] signal-peptidase-type-1 motif** the ICTV profile
describes, with the cut after the second Ala.

**The residue numbers do not travel; the motifs do.** A 162 nt indel in ORF3
separates the YHV subtypes, so pp3 runs 1640–1666 aa and site 2 sits at 1101 in
GAV, 1104 in GAV genotype 3, 1117 in YHV-8 and 1127 in YHV-1. Cutting on a
column number would mis-cut every genotype but one.

Derivation refuses any precursor where a site does not match exactly once.
9 of 11 cut; 2 refused.

```
precursor           len   p25  gp116  gp64   sum-precursor
YHV-1 (x5)         1666   228    899   539        0
YHV-8              1656   228    889   539        0
GAV                1640   228    873   539        0
GAV genotype 3     1642   228    876   538        0
```

**p25 is 228 aa in all nine**, the published value; gp64 is 539 in eight and
538 in the ninth; gp116 carries all the variation, as the indel predicts.

---

## 6. The frameshift — two heptamers, read off the genomes

The ICTV profile gives the family's slippery sequence as "AAAUUUU or GGGUUUU".
Measured at the pp1a/pp1b junction, the two genera use **different** ones, and
the Nodaviridae lesson applied: a single assumed heptamer would have covered
half the work.

```
Okavirus     AAATTTT   GAV, YHV-1, YHV-8 — one hit each near the junction
Nimanivirus  GGGTTTT   both MrGV variants
```

For MW590703 the measured `GGGTTTT` starts at **nt 14925**, exactly the
coordinate PMID 35011135 publishes (5'-G14925GGUUUU14931-3').

Reference sets: Okavirus **7** (6673–6692 aa, matching the literature),
Nimanivirus **3** (all 7797 aa). Nimanivirus is thin because only one genome
deposits pp1ab; the other two were built from pp1a + ORF1b pairs, and 19 further
genomes were refused for ambiguous bases.

---

## 7. What these modules do NOT cover

### Eleven divergent amnoonvirid species

Flavolineata, Piscibus, Lipesoluc and the Catfish / Cod / Flatfish / Loach /
Lizard ×2 / Pufferfish amnoonviruses deposit **one segment each — PB1/RdRp only**
— and none has a second segment anywhere. Seven of them fall inside a ±25%
window around the TiLV PB1 length and would have trained the PB1 profile on
viruses this module does not route; they are excluded **on the species**, which
is what the scope decision is about, not on length.

Measured: **0 of their 12 contigs route to Tilapinevirus** at the 150-bit floor.
They cannot support a module of their own (one genome each, one segment each)
and are recorded as a gap.

The Loach amnoonvirus PB1 (821 aa) is so divergent it does not even blastp to
the TiLV reference at e<1e-3, though it carries the SDD polymerase motif.

### No mature peptides for either ronivirus replicase

The ICTV profile states the number and composition of the nsps released from
pp1ab **is not yet known**. BV-BRC holds deposits annotated `3C-like protease`
(335–343 aa), `helicase` (660–846 aa) and `RNA polymerase` (92–286 aa), but
these are domain fragments with no published junctions. Only the 3CLpro
*specificity* VxHE↓(L/V) is characterised (PMID 12502857) — and a specificity is
not a coordinate. Declaring a mature peptide would mean inventing its
boundaries, so none is declared.

### MrGV ORF3 ships whole

ICTV asserts it encodes gp116 and gp64 homologues. No cleavage site is published
for MrGV and the okavirus sites cannot be transferred (§1). The precursor is
declared whole and the products are a named gap, not an invented cut.

### ORF4 is GAV-only, and its size was settled by measurement

The ICTV profile gives 98 aa; the Walker & Sittidilokratna review gives 83 aa
(252 nt). **Both GAV deposits in this dump are 83 aa**, so 83 is used and the
disagreement is recorded rather than resolved silently. No `copy_num`: it is
genus-partial by construction, so marking it essential would flag every yellow
head virus genome.

### No splicing, no PA-X, no matrix protein in TiLV

PMID 41513822 searched for all three and reports none; S9-F3 arises by leaky
ribosomal scanning. The probe file tests for them explicitly, and ABSENT is the
expected and correct verdict — it confirms no orthomyxovirus splicing rule has
leaked into this module.

---

## 8. Results

### Held-out panels

References were excluded by name. Tilapinevirus: 21 held-out whole genomes of
29 (8 are references). Roniviridae: 36 held-out complete genomes of 44.

| module | genomes | good | poor | exact duplicate calls |
|---|---|---|---|---|
| Tilapinevirus | 21 | 20 (95.2%) | 1 | **0** |
| Okavirus + Nimanivirus | 36 | 33 (91.7%) | 3 | **0** |

Per-feature call rate, counted by the module the annotator itself routed each
genome to:

```
Tilapinevirus (21)   PB1 PA NP S5 S6 S7 S8 S9 S9-F3 S10  100%
                     PB2                                  95.2%

Okavirus (6)         PP1A N PP3 P25 GP116 GP64 PX        100%
                     PP1AB                                83.3%
                     ORF4                                 16.7%  (GAV-only: 1 of 6)

Nimanivirus (30)     N ORF3 PP1A                          100%
                     PP1AB                                 6.7%
```

**S9-F3 — the protein no database holds — is called on 21 of 21 held-out
genomes.**

### Nimanivirus PP1AB is genome quality, not the module

```
genomes with NO IUPAC ambiguity codes :  1, PP1AB called in 1
genomes WITH ambiguity codes          : 28, PP1AB called in 0
```

28 of 29 held-out Nimanivirus genomes are TPA assemblies carrying IUPAC
ambiguity codes, and the frameshift correction needs an unambiguous codon at
the slip site. The three "poor" genomes are flagged `Contig has too many
ambiguous bases` for the same reason.

### Cleaved-end tiling — the gate for a multi-product precursor

Sum of the three pp3 products minus the precursor, per held-out genome:

```
5 of 6 genomes:  sum − precursor = 0, every junction gap 0,
                 gp64 ending exactly at the precursor's stop
1 of 6 (YHV-8):  +1, a one-codon overlap at the p25|gp116 junction
```

The single drift is a divergent genotype where the GP116 profile extends one
residue past the motif cut. Recorded, not hidden.

### Overlapping calls — all of them expected architecture

0 exact duplicates in either family. Within-contig overlaps:

```
Tilapinevirus   S9 / S9-F3                     21   different frames, by design
Okavirus        PP1A / PP1AB                    7   the frameshift
                P25, GP116, GP64 inside PP3     6   cleavage products
                PX inside PP3 / P25             6   ORFX overlaps ORF3's 5' end
                GP116 / P25                     1   the YHV-8 codon drift above
```

### Whole-taxon coverage, on 95% exemplars

**Tilapinevirus — 39 exemplars of 546 in-scope contigs**

```
routed                                      39 / 39   100%
expected feature called                     33 / 39   84.6%
wrong feature called                              0
cross-segment leakage                          NONE
```

Every miss is explained exactly: **33 of 33 contigs encoding ≥65% of their
expected protein are called (100%)**, and all 6 misses sit at **0.23–0.33**
coverage — diagnostic RT-PCR amplicons, below the declared `coverage_cutoff` of
0.65. The module does not call a 419 aa protein from a 295 nt amplicon, and
should not.

**Roniviridae — 17 exemplars of all 234 contigs**

```
routed            17 / 17   100%,  misrouted 0
called something  11 / 17
```

All six that called nothing cover **only ORF1b or the ORF1a/ORF1b junction**,
which this module declines to model by design. KX883632 (5040 nt, GAV) is the
instructive one: it routes at 6216 bits and every profile returns `stop-free 0`,
because the contig spans the frameshift and no single reading frame is open
across it.

### Cutoffs — measured, not inherited

Generated cutoffs are floors (30–60 bits). Raised to half the lowest genuine
panel score, capped at 500, features with fewer than five panel calls left
alone:

```
Tilapinevirus  11 raised   S9-F3 101 … PB1 500
Okavirus        7 raised   PX 90 … PP1A/PP3/GP116/GP64 500
                ORF4 LEFT ALONE (2 calls — one low score is not a distribution)
Nimanivirus     3 raised   N 134, ORF3 500, PP1A 500
```

Re-running both panels after the raise gave **identical** good/poor counts and
call rates — no regression.

### copy_num — set only after measurement

```
Tilapinevirus   PB1, PB2, PA, NP
Okavirus        N, PP3
Nimanivirus     N, ORF3, PP1A
```

**S9-F3 was withheld despite measuring 100%.** `apply_copy_num.py` proposed it,
correctly noting every call is *out of frame* with the S9 ORF containing it and
so is a separate gene rather than a cleavage product. The reason for withholding
is biological: PMID 41513822 names isolates that genuinely lack the protein, so
marking it essential would flag correctly-sequenced genomes as defective. The
100% is of a panel that happens to exclude the known exceptions.

**GP64 was withheld** for consistency with P25 and GP116, which the tool blocked
as lying wholly inside PP3. All three are products of one precursor and are
called exactly when PP3 is; PP3 carries the essentiality.

The uncharacterized features (S5–S10, PX, ORF4) are never marked, per the
project rule.

---

## 9. A kit defect found and fixed

`check_literature_gaps.py` documents a `RECOGNISED` verdict — "deposits found,
bound to NOT_MODELLED, a deliberate call" — and tests for it with
`dests <= {NOT_MODELLED, DISCARDED}`. But `collection_engine.py` wrote
**`UNASSIGNED`** into `synonyms.tsv` for every NOT_MODELLED string, collapsing a
curator's decision into "nobody wrote a rule for this". **The RECOGNISED verdict
was unreachable**, and every deliberate decision came back `UNBOUND` — "a defect
in your rules" — to be waved through with `--allow-absent`, a flag that means
something else entirely.

Found on Okavirus, where the nsp domain fragments are declared not-modelled on
the ICTV's own statement that the nsp composition is unknown.

Fixed in two places: the engine now writes `NOT_MODELLED` (and `DISCARDED` for
partials), and `synonyms_table.py` excludes `NOT_MODELLED` from the *bound*
count, since those deposits are not called by the module and counting them would
overstate the collapse figure the manuscript reports. The Okavirus binned figure
is unchanged at 84.3% (188 of 223), confirming the fix is behaviour-preserving
for the headline number.

---

## 10. Build parameters and known limits

- **`m` lowered 5 → 2** on every feature. These are small, homogeneous
  collections; at the documented default most features built nothing. The
  identity floor was **not** lowered — it stays at 0.80, and identity still does
  the quality work.
- **Single-sequence profiles** via the leftover path at `--min-seqs 1` for
  Okavirus ORF4 (1 unique sequence) and Nimanivirus PP1A and N. Their cutoffs
  are measured, never inherited — except ORF4, which has too few panel calls to
  measure and keeps its generated floor.
- **Panels are small**: Okavirus has 6 held-out genomes because only 9 complete
  okavirus genomes exist in BV-BRC. Nimanivirus has 30, but 28 carry ambiguity
  codes.
- **Segment length floors are deliberately permissive.** `segments` feeds the
  annotator's contig-length gate, which runs *before* routing; a floor set from
  a segment's full length would reject every partial CDS deposit — the majority
  of this family — silently, having never offered them to a module. That is the
  Adenoviridae failure the skill records. Completeness is judged by which
  features are called, not by contig length.
- **`PMID` is empty on every feature.** All 19 citations sit in
  `PMID_claude_generated` awaiting a human curator's read.
