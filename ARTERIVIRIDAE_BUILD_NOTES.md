# Arteriviridae — build notes

Built 19 November 2026 with the LowVan kit at `~/lowvan_runtime/kit-verify`.
One module, `Arteriviridae`, plus the pre-existing `Simarterivirinae`.

Local build tree: `~/lowvan_scratch/arteri` (local disk; Box is synced at the
end, per the work-local rule). Runtime repo used for every measurement here:
`~/lowvan_scratch/arteri/VA_arteri`.

> **Read this first if you are re-running anything.** `_infrastructure/env.sh`
> sets `LOWVAN_DATA_DIR` to the **Box** checkout, and both
> `annotate_by_viral_pssm.pl` and `annotate_by_viral_pssm-GTO.pl` resolve their
> rep contigs and PSSMs from that variable. A hand-rolled
> `New-annotate-viral-taxon.pl` call therefore routes every genome against
> *every* module's reference contigs, read over Box: measured 8 genomes in 25
> minutes at **0.0% CPU**, and the processes could not be killed afterwards
> (SIGKILL included) until the Box reads timed out. `env.sh` in this tree now
> ends with `export LOWVAN_DATA_DIR=$LOWVAN_TESTREPO`. Prefer
> `evaluate_coverage.py`, which sets it itself from `--repo`.

---

## 1. The dump

| | |
|---|---|
| genome records in family Arteriviridae | 5,517 |
| `id_md5` rows | 30,408 |
| unique protein sequences | 17,689 |

Genera: Betaarterivirus 5,102, Alphaarterivirus 362, Gammaarterivirus 10,
Kappaarterivirus 9, Nuarterivirus 2, Muarterivirus 2, Lambdaarterivirus 2,
28 unassigned.

**Taxon 28344 is a legacy bin holding both PRRSV-1 and PRRSV-2.** Any
per-species comparison that groups on that taxon id mixes the two lineages;
apparent "low conservation" inside it was traced to exactly that and is not a
contamination finding.

## 2. Features

27 declared. Three groups:

* **Structural**, PSSM-called: E, GP2, GP3, GP4, ORF5A, GP5, M, N
* **Replicase and its mature peptides**, PSSM-called: ORF1A, NSP1, NSP1A,
  NSP1B, NSP2, NSP3, NSP4, NSP5, NSP6, NSP7A, NSP7B, NSP8, NSP10, NSP11, NSP12
* **Non-collinear, `special: transcript_edit`**, nucleotide references and no
  PSSM: ORF1AB, ORF1ATF, ORF1AN, NSP9

Gene order: `ORF1a-ORF1b-GP2-GP3-GP4-ORF5a-GP5-M-N`, with E overlapping GP2.

## 3. The two defects this build fixed

### 3.1 M and GP4 shared a bin

The M rule matched `(un)?glycosylated membran[e]?`, which also matched GP4's
own deposits. 35 GP4 deposits sat in the M bin — **2% of the bin**, which
misannotated **44% of genomes**: 268 of 612 called M at GP4's locus.

Fixed by tightening the rule to `unglycosylated…` and letting the homology
rescue arbitrate; all 36 bare-string deposits went to GP4 by sequence. M
1,686 → 1,650, GP4 1,741 → 1,777. Profiles M.4 and M.11 retired.

Verification is positional, not an absence of findings:

| | before | after |
|---|---|---|
| M called within 30 nt of GP4's start | 268 / 612 | **0** |
| M between GP5 and N, its true locus | — | **604 / 604 = 100%** |
| `qc_duplicate_calls.py` EXACT | GP4 + M | **none** |

This is the standing lesson that **contamination damage is not proportional to
the contaminated fraction**. Count foreign members, never the fraction.

### 3.2 The E profile fired a second time inside nsp2

`copy_num` blocked E at 98.7% raw against 2.8% standalone. E was being called
twice in 87 of 611 genomes: once correctly, once ~8.2 kb upstream inside nsp2.
Genuine E calls are 68–85 aa; the spurious ones were 150–199 aa.

Scoring every E call against all nine E profiles separates them completely:

| | n | best bits |
|---|---|---|
| spurious, inside ORF1a | 87 | 13.3 – **27.6** |
| genuine | 601 | **30.3** – 132 (median 99) |

`bit_cutoff` was **25**, below the spurious maximum. Raised to **30**: it loses
none of the 601 and removes all 87. Simarterivirinae's E was already at 30.

**A kit-level finding fell out of this.** `max_len` is declared on every PSSM
feature and is **never enforced** on a PSSM call — `annotate_by_viral_pssm.pl`
uses it only for genome size and for `non_pssm_partner` features. E declares
`max_len: 89` and was emitting 199 aa calls. `min_len` is used, but only as a
sanity check on read-through termination.

## 4. Transcript-edited features

Four products are reached by ribosomal frameshifting and cannot be called by a
collinear profile:

| feature | mechanism | site |
|---|---|---|
| ORF1AB | −1 PRF | `UUUAAAC` (PRRSV/LDV), `GUUAAAC` (EAV) at the ORF1a/1b junction |
| ORF1ATF | −2 PRF | `GGUUUUU` + downstream `CCCANCUCC` within 34 nt, inside nsp2 |
| ORF1AN | −1 PRF | same nsp2 site |
| NSP9 | derived from ORF1AB | |

### 4.1 The ORF1ab junction convention was wrong in the shipped set

The existing references used five different insertion offsets at a uniform
phase, producing **different junction residues**. UniProt **P19811** (EAV
replicase, reviewed) annotates `VAR_SEQ 1728..3175 "Missing (in isoform
Replicase polyprotein 1a)"`, so residue 1727 is the last shared with pp1a and
1728 is the first made in the −1 frame. That fixes the convention:

```
+6  -> L-K-L   WRONG, drops the Asn
+7  -> L-N-L   correct
+8  -> L-N-L   correct (what EAV uses; matches P19811)
+9  -> L-N-W   WRONG, inserts a Trp
```

The mechanism is tandem slippage on X_XXY_YYZ: the A-site tRNA is already
charged and selected by YYZ and **keeps its amino acid** after slipping onto
YYY, so Asn is incorporated before the frame changes. I had argued from first
principles that the A site re-decodes, which would give Lys; P19811 settles it
the other way.

**1,132 of 1,194 existing references were wrong.** The rebuild corrected 1,138
by exactly one residue: **W→L 804, K→N 338**. EAV re-verified against P19811 at
3,167 / 3,175 with an identical junction (`NQLNLRAP`).

### 4.2 References built from genomes, not only from deposits

`build_transcript_edit_refs.py` can only build a reference for a genome that
deposits the product, because it validates by reproducing that protein. For the
nsp2 PRF products that was the binding constraint: 540 ORF1aTF references for
2,374 complete genomes, and the mechanism only calls a product when a reference
**≥99% identical** exists. The cliff is sharp — ≥99.6% identity gives exactly
one gap run, ≤97.7% gives 3 to 69, and a mixed-direction set is rejected.
Failure tallies before the rebuild: ORF1AB 45,820 / 46,256 rejected as mixed
(99.1%), ORF1AN 96.8%, ORF1ATF 95.1%, NSP9 51.7%.

A genome does not have to deposit the protein to testify to it; it has to be a
good genome, and that is checkable structurally. `build_te_refs_from_genomes.py`
accepts a genome only with a canonical `GGTTTTT` carrying the conserved
downstream `CCCANCTCC` inside the 34-nt window (PMID 24825891), exactly one
complete PRF signal, an ORF1a start matching a deposited N-terminus at ≥80%
over 20 residues, a clean terminator, no ambiguous base, and a product inside
the literature window. **What is written is the corrected transcript** — the
duplication already applied — because the annotator transfers the reference's
extra base in as a subject gap, and a raw genomic span carries no edit.

The insertion convention was measured, not guessed: all 540 existing ORF1ATF
references duplicate `TT` at slip + 7 (539 of 540; one at 8). Reproduction of
the existing set was **536 / 536 byte-identical**.

Start anchoring matters and fails quietly without it. "The first in-frame
stop-free ATG" gave 104 distinct start positions spanning 0–403 nt on taxon
28344, picking the real start only 34.5% of the time; `start 0` means the
contig begins mid-ORF. Anchoring on the 115 distinct deposited ORF1a
N-terminal 20mers rejected 29 genomes and caught one truncated reference
already shipping.

### 4.3 nsp2TF and nsp2N are lineage-restricted

**486 of 615 exemplars carry the nsp2 PRF signal.** EAV (0 of 10) and Wobbly
possum disease virus (0 of 5) genuinely lack it. So **486, not 612, is the
honest denominator for ORF1aTF and ORF1aN**, and any report that uses 612
understates them.

## 5. Positional QC

Absence of a duplicate-call finding is not evidence a feature is called in the
right place. Measured over the 612-genome evaluation:

**The mature-peptide cassette is in the correct 5′→3′ order in every genome.**
Every consecutive pair is 100.0%:

```
NSP1A->NSP1B 586/586   NSP4->NSP5  600/600   NSP7A->NSP7B 569/569   NSP10->NSP11 607/607
NSP2 ->NSP3  596/596   NSP5->NSP6  578/578   NSP7B->NSP8  571/571   NSP11->NSP12 577/577
NSP3 ->NSP4  609/609   NSP6->NSP7A 581/581   NSP8 ->NSP10 584/584
```

That is the evidence that **nsp6 is callable at all**. It is 15–23 aa from a
5-member collection, which is normally too short to profile; here it is correct
578/578 after nsp5 and 581/581 before nsp7α, because it sits in a fixed
cassette whose neighbours anchor it.

**The structural cassette** is 100% for GP2→GP3→GP4→ORF5a and GP5→M→N. Two
apparent exceptions, both resolved rather than assumed:

* **E before GP2 in 90 genomes** at an offset of **−8,198 nt**. Not an overlap
  — the spurious nsp2 calls of §3.2, found independently.
* **ORF5a vs GP5 in two discrete classes**, −10 nt (n=409) and +5 nt (n=170),
  with almost nothing between. The −10 class is the documented PRRSV
  arrangement. Both produce well-formed ORF5a — 43–51 aa, Met-initiated, with
  the basic `RKCLRRR` / `CQRQRQ` motif:

  ```
  -10  MLQHVGQVLDSCLLLLVAFFVVYRLVLFLCQRQRQQLNDAIHI
   +5  MFKQIGRVFDTILAIAILFCLVYSLYRKCLRRRRQLNQELDLHLNLV
  ```

  Two lineages with different start arrangements, not one protein called in two
  places. Recorded, not "fixed".

## 6. Nested calls are not duplicate calls in a polyprotein module

`qc_duplicate_calls.py` reports 21 NESTED rows and **none is a defect**. Its
verdicts are written for modules without polyproteins and are wrong here:

* 13 rows are mature peptides inside ORF1a, blastp 100% identity over 100%
  coverage — which is what a cleavage product *must* give. "SAME PROTEIN,
  merge" would delete the mature peptides the module exists to call.
* 5 are ORF1A and nsp9–nsp12 inside ORF1AB, the same relation one level up.
* ORF1ATF and ORF1AN inside ORF1A share their entire N-terminus by mechanism.
* **ORF5A inside GP5** (id 52%, cov 43%, 22 bits) and **E inside GP2** (id 31%,
  cov 60%, 20 bits) are genuine overlapping ORFs in different frames. The low
  scores are the proof: a promiscuous profile calling one protein twice gives
  high identity, as the nsp rows do.

`apply_copy_num.py` reaches the same conclusion mechanically for ORF5A, which
it marks essential despite containment because **0 of 602** calls are in GP5's
frame.

## 7. `Arteriviridae sp.` (taxon 2697055) is out of scope

54 of the 5,954 genomes in the dump carry the organism name
`Arteriviridae sp. <isolate>` -- taxon 2697055, a single BV-BRC catch-all bin
holding rodent-derived assemblies (`RtRe-ArtV/...`, `RtBi-ArtV/...`,
`GX/B30F0408/2019`, `MgAV1`). The skill's intake rule drops an `sp.` taxon
before collections are built: it has no assigned species, so its members are
not known to be homologous to one another.

**It never reached a profile.** Taxon 2697055 has **zero** rows in
`Art.id_md5` -- BV-BRC holds the assemblies but no protein features for them --
so nothing of it entered a collection, an alignment or a PSSM. There is no
contamination to undo.

**It did reach the coverage audit**, because that is built from `Contigs`
rather than from the protein dump. 35 of the 700 exemplars (5.0%) are taxon
2697055. They are excluded from the coverage denominator, which is therefore
**661 in-scope exemplars**, and the exclusion is applied at report time rather
than by deleting contigs, so the run stays reproducible.

The exclusion matters for more than the denominator: these 35 genomes produced
**6 of the 12** out-of-window calls among the `copy_num` features -- all three
GP2 over-length calls, the GP5 under-length call and both N under-length calls.
Scored as in-scope they would have argued for widening four windows to cover
assemblies that are not known to be arteriviruses of any species.

**No other unclassified taxon is present.** Scanning every exemplar's organism
name for ` sp.`, `unclassified` and `unidentified` returns taxon 2697055 and
nothing else.

## 8. Length windows were re-derived before `copy_num` was applied

`min_len`/`max_len` are **not** tested by `annotate_by_viral_pssm.pl`; a window
cannot suppress a call. They are consumed by `viral_genome_quality.pl`, which
compares the protein translation to the window and flags
`Feature is too short` / `Feature is too long`, marking the feature *and the
genome* Poor -- and it reads the window **only for features carrying
`copy_num`** (`$essential`, line ~111).

So the window has to be right **before** `copy_num` goes on, not after.
Applying `copy_num` to the seven earned features against the windows as they
stood would have converted correct calls into quality flags.

Each window was reset to the literature range in `length_refs.json`, extended
only where the observed in-scope distribution required it:

| module | feature | was | now | literature | in-scope calls |
|---|---|---|---|---|---|
| Arteriviridae | GP2 | 175-322 | **150-320** | 150-320 | 609 |
| Arteriviridae | GP3 | 153-285 | **130-280** | 130-280 | 609 |
| Arteriviridae | GP4 | 76-196 | **80-210** | 80-210 | 606 |
| Arteriviridae | GP5 | 151-270 | **160-285** | 160-285 | 607 |
| Arteriviridae | M | 148-190 | **150-190** | 155-190 | 614 |
| Arteriviridae | N | 94-142 | **95-140** | 100-140 | 615 |
| Arteriviridae | ORF5A | 40-63 | **40-80** | 40-80 | 605 |
| Simarterivirinae | GP2 | 148-249 | **150-320** | 150-320 | 45 |
| Simarterivirinae | GP3 | 112-286 | **120-280** | 130-280 | 42 |
| Simarterivirinae | GP4 | 142-192 | **80-210** | 80-210 | 45 |
| Simarterivirinae | GP5 | 216-294 | **160-285** | 160-285 | 46 |
| Simarterivirinae | M | 149-183 | **155-190** | 155-190 | 46 |
| Simarterivirinae | N | 90-139 | **100-140** | 100-140 | 46 |
| Simarterivirinae | ORF5A | 42-72 | **40-80** | 40-80 | 45 |

M and N are the only two widened past the literature, both at the minimum and
both by a few residues, to cover measured calls (M 150 aa, N 96 aa).
Simarterivirinae GP5's declared minimum had been 216 against a literature
minimum of 160 -- a window drawn from the collection rather than from the
papers, which is what made a literature-plausible 176 aa call look like a
defect.

**The audit must be keyed on `(module, feature)`.** A feature key is declared
in both modules with different windows, and the GTO's `family_assignments`
carries the *profile* id `Module.FEAT.N`. Keying on the annotation string, or
on the bare feature name, merges the two modules: the first pass of this
measurement did exactly that and reported **436 phantom GP2 violations** by
scoring Arteriviridae calls against the Simarterivirinae window.

### 8.1 The six calls that remain outside their window

Across 8,565 calls in 661 in-scope genomes, six fall outside. Each is a real
per-genome anomaly, and a Poor flag is the correct response to all six:

| feature | call | genome | what it is |
|---|---|---|---|
| GP3 | 290 aa | 28344.18027 | PRRSV; same start as its neighbours, ~111 nt further 3'. The sibling PRRSV genome 28344.37248 calls GP3 at 254 aa from the same start, so this is a C-terminal extension in one isolate, not a profile defect |
| GP3 | 297 aa | 1806636.5 | Rodent arterivirus RtMc-Arterivirus/Tibet2014, a divergent lineage |
| ORF5A | 35 aa | 28344.37248 | PRRSV; truncated against a 51 aa norm |
| ORF5A | 83 aa | 3084900.3 | Wufeng rodent arterivirus 1 |
| ORF5A | 83 aa | 3084900.4 | the same, second isolate |
| GP4 | 270 aa | 2682610.3 | **not a length problem** -- see section 9 |

The two Wufeng ORF5a calls sit exactly where ORF5a belongs: starting inside
GP4's last codons (11710 vs GP4 ending 11710) and overlapping the GP5 start at
11729, which is the arrangement Firth 2011 describes. They are well formed and
very likely genuine, and the window was still **not** widened to swallow them,
because 83 aa exceeds every length that paper reports. Widening a literature
window to cover an un-attested value would launder an observation into a fact;
two Poor flags on a divergent rodent lineage is the honest outcome.

## 9. Two genomes call a paralogue pair on one ORF

Simarterivirinae carries a duplicated minor-glycoprotein cassette, so GP2/GP2P,
GP3/GP3P and GP4/GP4P are paralogues. In two of the 661 in-scope genomes a
single interval is claimed by more than one of them:

```
2682610.3    GP3 + GP4 + GP4P   all at 12884..13696 (270 aa)
1885929.26   GP2 + GP2P         one interval
```

The 2682610.3 sequence is plainly GP3 -- signal peptide, `NYTTC` and `NITLHG`
glycosylation sites -- so GP4 and GP4P are the spurious callers, and the 270 aa
"GP4" in section 8.1 is this collision rather than a long GP4.

Two genomes in 661 is 0.3%, and no other feature set collides anywhere in the
panel. The fix is the standard one for paralogues -- position curates the call,
profile per orthogroup -- and it is **not** a cutoff change, because both
profiles score the same ORF honestly.

**Consequence for `copy_num`:** the P-paralogues are lineage-specific
accessories and must not carry it, which was already the rule; this measurement
is the evidence for it. A `copy_num` on GP4P would mark every Arteriviridae
genome outside Simarterivirinae as missing an essential feature.

## 10. NSP10 is clipped by ~70 residues outside PRRSV

**Measured 20 Nov 2026 on the 661 in-scope exemplars. 38 NSP10 calls fall
below the literature minimum of 420 aa, in a tight 377-406 band.**

Length alone cannot say whether a short call is a short protein or a truncated
one. Adjacency can, because nsp10 and nsp11 are consecutive cleavage products
of pp1ab and must abut:

```
normal  NSP10 -> NSP11   gap 0 nt   565 of 570
normal  NSP10 -> NSP11   gap 12 nt    5
short   NSP10 -> NSP11   gap 219 nt  10      <- EAV
short   NSP10 -> NSP11   gap 138-234 nt  25
short   NSP10 -> NSP11   gap 0 nt     3
```

The short calls leave a **gap**. They are clipped, not short.

EAV settles it arithmetically. All ten EAV genomes call NSP10 at exactly
399 aa with a 219 nt gap to NSP11; 219 nt is 73 codons, and 399 + 73 = 472
against the 467 aa the literature records for EAV nsp10 (PMID 9971783,
2371-2837 of pp1ab). The profile is stopping about seventy residues before the
nsp10/nsp11 junction.

**Cause: one profile, built from a collection that is almost entirely PRRSV.**

```
NSP10 collection: 65 members
   441 aa  x62     (PRRSV)
   440 aa  x1
   442 aa  x1
   467 aa  x1      (EAV -- one sequence)
NSP10 profiles:   1  (Arteriviridae.NSP10.1)
NSP10 config:     downstream_ext 0, upstream_ext 0
```

With `downstream_ext: 0` a call is exactly the profile's match extent, so a
divergent lineage is called across the conserved core and cut wherever
similarity falls away. The affected lineages are precisely the ones the
collection barely represents: EAV (10), Wobbly possum disease virus (5),
Lactate dehydrogenase-elevating virus (5), Lopma virus (2), Praja virus (2),
Wufeng rodent arterivirus (2), Olivier's shrew virus 2 and 3, and one PRRSV.

**This is a collection gap before it is a profile gap.** Splitting the existing
profile cannot help: there is only one non-PRRSV sequence in the bin to split
off. The fix is to add nsp10 members for the divergent lineages and rebuild
per orthogroup. The material to do it with is already in hand -- in each
affected genome NSP11 is called correctly, so the missing residues are exactly
the interval between the NSP10 call's end and the NSP11 call's start.

**The window was deliberately NOT widened to absorb these calls.** 420-490 is
the literature range and the calls are wrong, so they flag Poor, which is the
correct outcome until the profiles are fixed. Widening to 380 would make 38
truncated calls look clean and would hide the defect in exactly the way a
length window is most dangerous.

**This was repaired on 20 Nov 2026 — see section 12.** The collections were
rebuilt from projected cassettes covering 22 lineages instead of 2, and the
window moved only to 405 to admit Kappaarterivirus at a measured 413 aa, not
to 380 to hide the clipping.

**Contrast with NSP7A, section 11**, where the same symptom had the opposite
cause. Short length is not evidence of anything on its own; the neighbour's
coordinate is.

## 11. NSP7A at 117 aa is genuine, and its literature entry was wrong twice

19 complete PRRSV genomes call NSP7A at exactly 117 aa against a literature
PRRSV figure of 149 aa. By the section 10 test these are **not** clipped:

```
NSP6 ends  -> NSP7A begins   gap 0 nt   20 of 20
NSP7A ends -> NSP7B begins   gap 0 nt   13 of 20   (3, 6, 12 nt in the rest)
```

Flush on both sides, so the cassette is contiguous and 117 aa is a real
shorter nsp7alpha in this lineage.

`length_refs.json` had `NSP7A` at **130-170**, and its own note records
*"EAV nsp7alpha 123 aa"* -- the declared window excluded the reference virus
the entry cites. It was corrected to 120 for that reason and then to **115** to
admit the measured 117 aa class. Both corrections are recorded in the entry
itself with their evidence.

## 12. NSP10/11/12 were rebuilt from projected cassettes (the fix for section 10)

Section 10 diagnosed the clipping. This section records the repair.

### 12.1 The collections held two taxa out of twenty-two

| feature | members before | taxa before | members after | taxa after |
|---|---|---|---|---|
| `NSP10` | 65 | 28344, 11047 | 97 | 22 |
| `NSP11` | 9 | 28344, 11047 | 40 | 22 |
| `NSP12` | 8 | 28344, 11047 | 39 | 22 |

Every ORF1b nsp in the dump comes from PRRSV (28344) or EAV (11047): 62 of 65
NSP10 members are 441 aa PRRSV and exactly one is the 467 aa EAV. No other
lineage ever deposited an ORF1b mature peptide, so no amount of re-binning
could have fixed this — the members did not exist and had to be derived.

### 12.2 Why the kit's projector could not be used

`project_matpeptides.py` refuses to transfer a cut below `--min-id 0.80`, and
its own help says why: *"at 0.50 the alignment is no longer trustworthy
residue-by-residue and the cut lands wherever the gaps fell."* The affected
lineages are **0.34-0.51** identical to the only two references, so the tool
cannot reach them. That floor is correct **for a method in which the alignment
places the cut.**

### 12.3 Alignment localises, chemistry places

The cut was placed in two stages: a global alignment to the nearest reference
cassette to get within a dozen residues, then a snap onto the nearest nsp4
(3CLSP) cleavage site. Arterivirus nsp4 cuts after Glu, and after Gln at EAV's
nsp10/11 site (PMID 9971783).

The chemistry is nearly unique over that interval, which is what makes the
weak alignment sufficient. Measured on the **573 exemplars whose nsp10 call is
independently correct**:

```
E/Q-followed-by-small sites within +/-12 of the TRUE cut
    exactly one   562   (98.1%)
    two            10
    none            1
```

Validation, same 573 genomes, projection run blind against them:

```
nsp10 start reproduced exactly   570 / 573
nsp10 length reproduced exactly  562 / 573
nsp11 length reproduced exactly  565 / 570
P1 at the projected junction     E in 572 / 573
```

EAV lands on **467 / 219 / 119** — the literature values (PMID 9971783) — with
the Gln at P1 that the same paper reports as unique to EAV.

All 35 repaired cassettes land on canonical sites:

```
nsp10|nsp11    E|G  22     E|S  4     Q|S  9  (the nine are EAV)
```

### 12.4 The nsp9/nsp10 end needed a second, structural anchor

P1' at nsp9|nsp10 is G in 568 of 573 correct genomes, but **K** in
Gammaarterivirus, Praja, Wufeng rodent, Bamboo rat and Hedgehog arterivirus,
and no E/Q-followed-by-small alternative exists within +/-12 in any of them.
Chemistry alone could not settle whether those starts were right.

The zinc-binding domain settles it. Nsp10 is helicase + ZBD and the ZBD is its
N-terminus:

```
cysteines in the first 45 residues of nsp10
    573 of 573 correct genomes          6 or 7
    all 35 projected starts             6 or 7
```

and the conserved `HxH...CPV` motif is intact in every E|K lineage, e.g.
Gammaarterivirus `KKKCRTCAHCGAPSTLVSSCGLNLCDYHGHGHPHCPVVLPCGHAV`. So E|K is
genuine lineage variation at P1', not a misplaced start.

### 12.5 Gates applied before a cassette was written

A projected cassette was accepted only when all of:

* P1 of all three junctions is E or Q
* nsp10 405-495, nsp11 195-245, nsp12 105-175
* 6-7 cysteines in nsp10's first 45 residues

P1' is deliberately **not** gated, per 12.4. 596 of 608 cassettes passed; up to
two genomes per taxon were written to `Extracted/`, giving 35 sequences per
feature across 22 taxa in all six genera.

### 12.6 One window changed

`NSP10` **420-490 -> 405-495**. Kappaarterivirus (wobbly possum disease virus,
7 genomes) is **413 aa** with both junctions canonical E|G and the ZBD intact;
the old floor excluded it. Recorded in `length_refs.json` with the measurement.
No other window moved.

### 12.7 `merge_and_json.py` clobbered 33 fields and they were restored

The merge recomputes `min_len`/`max_len` as 0.8x/1.2x of the collection median
and `bit_cutoff` as `max(20, 0.45x)` for **every** key with a collection, not
only the changed ones. It overwrote 33 fields across 19 features, including the
measured `E` bit_cutoff of 37 — the one that removed 35 spurious E calls —
and invented a `bit_cutoff` of 308 for `NSP9`, a `transcript_edit` feature with
no PSSM.

All 33 were restored from `Arteriviridae_Viral_PSSM.json.bak_nspfix` and only
the 12.6 change was applied. **Any future merge needs the same diff-and-restore
pass**; this is recorded in `KIT_PATCHES_2026-11-19.md`.

### 12.8 What is not claimed

The repair rests on a cleavage-chemistry argument, not on deposited boundaries,
for 20 of the 22 lineages. The measurements in 12.3 and 12.4 are the evidence
that the argument holds, and they are reproducible from `nspfix/`. Where a
lineage later acquires a deposited ORF1b mature peptide, it should be preferred
over the projection and the projection re-measured against it.

One genome is deliberately left unrepaired: `28344.44501`, a PRRSV whose nsp10
is **388 aa** with E|G at both junctions, flush neighbours and an intact ZBD.
The projection agrees with the call exactly. It is a genuine short form below
the window, not a clipped call, and it flags Poor.

## 13. Open items

* `build_transcript_edit_refs.py` chooses the insertion "closest to the
  pre-shift stop". That rule produced the §4.1 junction error here; whether it
  has put the same error into Tobaniviridae and the coronavirus modules has not
  been checked.
* The annotator lowercases the first letter of the original product name when
  it builds the demoted-partial string: `Uncorrected oRF1ab polyprotein
  encoding region`, `Uncorrected rNA-dependent RNA polymerase encoding region`.
  Cosmetic, but it ships into every GTO.
* 34 `hypothetical protein` sequences remain uncalled by decision, recorded in
  `TABLED_UNASSIGNED.md`.

* **GP5 profile 7 puts spurious residues on two calls.** `check_met_start.py`
  reports its alignment starting 4 residues into its members. The cause is
  real: its two members have genuinely different N-termini — `MLSRLFDCGDMK…`
  and `MGHCHHLLSRLLDCGDMK…`, six residues apart — so mafft trimmed to the
  common region. With the profile matching downstream of the true start,
  `upstream_ext: 1` scans past the real Met and stops at whatever it finds:

  ```
  2012619.3   F KYHRSSQRWCLDTSYSFPGNDQ VGHCHHLLSRLFDCGDMK…   (24 spurious residues)
  3084632.3   I QGFLFFSFLPACHAVVRSHCVHLCILLRSDTMTCSNNWARP…
  ```

  **602 of 604** GP5 calls begin at Met, so the effect is 0.3% — but spurious
  residues are never acceptable, and the fix the tool names is right: split
  cluster 7 on its N-terminus into two single-sequence profiles, each with a
  **measured** cutoff (they cannot inherit one). Deferred only so that profiles
  are not changed underneath a coverage run already in flight.

* **The Simarterivirinae paralogue collision (section 9) is not yet fixed.**
  Two genomes, 0.3% of the panel. The fix is per-orthogroup profiles with
  position curating the call, not a cutoff change.
* The 54 taxon-2697055 contigs are still present in `Contigs`. They are
  excluded at report time (section 7) rather than deleted, so any re-run of
  the coverage audit must apply the same exclusion or the denominator moves.

* **NSP10 is clipped outside PRRSV (section 10), 38 of 661 genomes.** Not
  fixed. Needs nsp10 members for EAV/LDV/WPDV and the other divergent
  lineages added to the collection, then profiles rebuilt per orthogroup.
  This is the largest open defect in the module.
* **NSP8 emits 1494 and 1502 aa calls** in two PRRSV genomes against a 35-60 aa
  window -- the same shape as the E defect, a weak hit extended far past the
  protein. Two genomes, correctly flagged Poor, cause not yet investigated.
* **ORF1A emits 3829 and 3922 aa calls** in the same two genomes, against a
  2700 aa maximum. Likely the same underlying readthrough.

### 12.9 The bin audit after the rebuild

`audit_bins.py` re-run over all 24 collections on 20 Nov 2026 after the
NSP10/11/12 rebuild: **Unresolved: 0** (11 FRAGMENT findings, all 11 signed
off; 1 paralogue bin; 24 collections).

The eleven FRAGMENT findings are the mature peptides of pp1a/pp1ab reported as
"partial of ORF1A" on identity alone. They are expected: a cleavage product
*is* a sub-sequence of its polyprotein, so 100% identity over 1-6% of ORF1A is
the correct reading of a correct bin, not contamination. Each carries its
literature size and junction chemistry in the sign-off (PMID 11602729,
PMID 30353315).

One SPLIT is **reported and deliberately not retired**:

| bin | n | median aa | finding |
|---|---|---|---|
| `NSP12` | 39 | 143 | 9 groups, largest holds 44% of the bin |

nsp12 is the C-terminal product of pp1b and the least conserved protein in the
cassette; the length audit puts it at 118-164 aa across the family. Nine groups
in 39 sequences is lineage structure, not a second protein — the groups track
genus, and every member passed the E|[GSA] junction and ZBD-anchored projection
gates in section 12.5. Splitting into per-orthogroup profiles is what
`build_features.py` does with the bin anyway, so the finding is recorded rather
than acted on separately.

Post-rebuild length audit rows:

```
NSP10   97 members   413 / 441 / 467     (window 405-495)
NSP11   40 members   219 / 223 / 225
NSP12   39 members   118 / 143 / 164
```

### 12.10 Verification: the cassette after the rebuild

`coverage6`, 612 exemplars, 612 routed, 0 unrouted. Measured with the same
script on both runs, so the two columns are comparable.

```
                          BEFORE (coverage5)        AFTER (coverage6)
NSP10 length         648  301 / 441 / 443      610  348 / 441 / 467
NSP11 length         646  154 / 223 / 225      609  154 / 223 / 225
NSP12 length         583  118 / 154 / 165      610   92 / 153 / 165

NSP10 -> NSP11      569 / 645 flush  88.2%    605 / 608 flush  99.5%
NSP11 -> NSP12      410 / 581 flush  70.6%    601 / 608 flush  98.8%
gapped junctions           247                        10

NSP10 ends in E      571 / 648  88.1%          597 / 610  97.9%
         or Q (EAV)        --                  607 / 610  99.5%
```

**The single most telling number is NSP10's maximum: 443 -> 467.** EAV's nsp10
is 467 aa in the literature (PMID 9971783) and no call in the old run reached
it, because one profile built from 62 PRRSV sequences capped every lineage at
its own extent. EAV's nsp12 now comes out at 119 aa against the 118 of the
literature cassette.

**The NSP11 -> NSP12 discrepancy recorded as an open item is closed.** It was
71% flush and is now 98.8%; it had the same cause as the nsp10 clipping and was
fixed by the same rebuild, not separately.

#### The 10 junctions still gapped

Eight are single genomes rather than lineages, and seven of those are PRRSV
deposits (`28344.*`) that are aberrant in their own right:

```
28344.35931  NSP11->NSP12  324 nt   NSP11 154 aa ends R
28344.25698  NSP10->NSP11  281 nt   NSP10 348 aa ends S
28344.26426  NSP11->NSP12  169 nt   NSP11 167 aa ends K
28344.12616  NSP11->NSP12   38 nt   NSP11 209 aa ends H
1806636.5    NSP10->NSP11   12 nt   NSP10 434 aa ends A
1806636.5    NSP11->NSP12    9 nt
1806636.3    NSP11->NSP12    6 nt
28344.26604  NSP10->NSP11    3 nt   NSP10 440 aa ends L
1965067.23   NSP11->NSP12    3 nt
28344.38494  NSP11->NSP12    3 nt
```

None ends in E or Q where the gap is large, which is the signature of a
truncated or frameshifted deposit rather than a profile that stopped early --
a clipped profile call still ends on the conserved core's last aligned residue,
and these do not. They are recorded, not repaired.

#### Three NSP10 calls do not end in E or Q

`1806636.5` (A, 12 nt gap), `28344.25698` (S, 281 nt gap) and `28344.26604`
(L, 3 nt gap). Down from 77 such calls before the rebuild.

#### One new outlier to resolve

NSP12's minimum fell from 118 to 92: `28344.39812`, called by `NSP12.5` at
92 aa. Its neighbours are flush, so this is not clipping. Open.

## The UNASSIGNED `hypothetical protein` pool, answered rather than tabled (21 Nov 2026)

`install_module.py` reports, every time it runs:

```
UNASSIGNED POOL: 1 string(s) at >= 5 uncalled sequences were never clustered
    34    hypothetical protein
```

The skill requires that each group reaching the floor be declared or the
decision to table it recorded. Declaring it needs to know what the sequences
are, so they were pulled out of the dump (98 unique sequences whose md5 is in
no Arteriviridae collection) and identified. **None of them is a missing
feature.** They fall into three groups, and each group is settled by a
measurement rather than by a judgement call.

**(a) 61 of 98 are already modelled — in the sister module.** Blasted against
`Simarterivirinae/collections/`, they hit `ORF5A` at **100% identity over
their full length** (e.g. `fig|1885929.10.peg.11`, 68 aa, 2.9e-48). They are
Kibale red colobus virus 1 and 2 ORF5a, verbatim members of Simarterivirinae's
own ORF5A collection. They appear in Arteriviridae's unassigned pool only
because both modules are built from one family-wide dump — the pool is
computed per workdir and does not know the sister module exists. **This is an
artefact of the shared dump, not a gap**, and it accounts for nearly two
thirds of the pool.

**(b) 30 of the remaining 37 are minus-strand ORFs — impossible, and therefore
discarded.** Every one of them was located on its own contig by six-frame
search. All 30 are on the **minus strand of a positive-sense ssRNA genome**:

| sequence family | aa | genomes | strand / frame |
|---|---|---|---|
| `MDRRTLTHVGKQVFPHKAKD…` | 134 and 197 | 15 | − |
| `MLIVSLVHGRNVIVAGQGGI…` | 81 | 7 | − |
| `MKFWPGVGGRENSRYLIQSA…` | 100 | 4 | − |
| `MQIILIRHIGRDSDVDRCSY…` | 60 | 3 | − |
| others | 67–101 | 1 each | − |

All 30 come from one block of PRRSV-2 deposits, `28344.384xx`, which carry
12–13 pegs where a PRRSV genome has ten. Arteriviruses have no minus-strand
coding capacity, so these are ORFs an annotation pipeline called on the reverse
complement. They are pipeline artefacts and are not binned. **The strand test
is what settles this** — by length, by name, or by homology they are simply
"unknown", and three of the four families are above the five-sequence floor and
would otherwise have demanded a feature.

**(c) 7 are genuine plus-strand ORFs, in 5 genomes, and are tabled with their
reason.**

| genome | organism | pegs | aa |
|---|---|---|---|
| `3084631.3/.4/.5` | Jingmen rodent arterivirus 1 | peg.3 | 310–312 |
| `2116439.3` | Hainan hebius popei arterivirus | peg.3 | 287 |
| `2867482.3` | Pelodiscus sinensis respiratory and haemorrhagic syndrome virus | peg.2, peg.4, peg.5 | 168, 260, 385 |

These are divergent lineages — rodent, snake and turtle arteriviruses — and
several carry signal-peptide-like N-termini (`MTTLKVFTILILAIALPQ…`,
`MAKSHVIILLAVILSVTAAQA…`), so they are plausibly the minor glycoproteins. That
is a guess, and it stays one: **each genome deposits only 3–5 pegs in total,
every one of them annotated `hypothetical protein`**, so the source never named
a single protein in these genomes. No literature places their boundaries
either. Under the source-annotation ceiling there is nothing to bin from, and
the largest group is 3 sequences — below the floor — so they are tabled rather
than modelled. They are recorded here as the module's known lineage gap.

## Held-out panel: 590 genomes, 530 of 587 scored Good (90.3%) (21 Nov 2026)

The panel is the 612 coverage-9 exemplars minus the 25 genomes whose contigs
are the module's own rep contigs, so no genome is scored against a reference
built from itself. **It needs no second annotation run**: the full
`evaluate_coverage.py` already annotates every exemplar through the GTO
pipeline and writes `Quality_GTOs/`, so the panel is a filter over measurements
that already exist. 587 of the 590 produced a quality GTO.

**530 Good, 57 Poor.** The Poor genomes split cleanly:

* **20 of 57 are short deposits** (<14,000 nt; the shortest is 12,639). Every
  one of the 530 Good genomes is >= 14,262 nt. These are incomplete genomes and
  the flag is correct.
* **37 are full length.** Their failures are single-feature and scattered over
  fifteen different features -- NSP1A 7, NSP1B 4, NSP11 4, ORF1A 4, NSP3 3,
  then ones and twos. There is no second systematic defect hiding in them;
  almost all are short calls on divergent isolates.

**Three genomes produce runaway calls, and the cause is in the deposit.**

| genome | call | length | |
|---|---|---|---|
| `28344.39183` | ORF1A / NSP8 | 3923 / 1495 aa | |
| `28344.39739` | ORF1A / NSP8 | 3830 / 1503 aa | |
| `11047.866` | NSP12 | 1260 aa | |

For the two PRRSV genomes this was tested rather than assumed: reading the
ORF1a frame codon by codon, **the first stop codon is at nt 11955 and 11666**,
roughly 4,000 nt past where a PRRSV ORF1a ends. The deposits have no authentic
ORF1a stop, so `downstream_ext: 1` scans on into ORF1b and both ORF1A and NSP8,
which share that terminus, run to the same place. The annotator is behaving
correctly on a defective sequence. `11047.866` is a different fault -- the stop
is exactly where the call ends (nt 9746), so NSP12's **start** is ~1,140 codons
too early -- and it is one EAV genome.

These counts are measured before the NSP1A/NSP1B repair is installed. The
largest single flag group on the panel is `missing NSP1-beta` (18) and
`missing NSP1-alpha` (14), which is the PRRSV-1 profile-reach defect, so the
panel is re-scored after coverage 10 rather than quoted from here.

## Every EAV genome was Poor for features EAV does not have (21 Nov 2026)

All ten Equine arteritis virus exemplars scored **Poor**, with the same three
flags:

```
missing Papain-like protease involved in polyprotein cleavage NSP1-alpha
missing Papain-like protease and interferon antagonist NSP1-beta
missing NSP7-beta protein
```

Two different faults are hiding behind one symptom, and they need opposite
remedies.

### NSP1A and NSP1B must not be essential: EAV has one nsp1

**EAV has a single, uncleaved nsp1, not nsp1α/nsp1β** -- den Boon 1995, PMID
**7769711**: "The replicase of the related virus equine arteritis virus
contains only one active PCP in the corresponding region." The module already
models this correctly: it declares a separate `NSP1` feature, and the calls
follow the biology exactly --

| taxon | n | NSP1 | NSP1A | NSP1B |
|---|---|---|---|---|
| PRRSV-2 `28344` | 511 | 0 | 510 | 509 |
| PRRSV-2 `1965067` | 34 | 0 | 34 | 34 |
| PRRSV-1 `1965066` | 16 | 0 | 16 | 16 |
| **EAV `11047`** | **10** | **10** | **0** | **0** |

`NSP1` correctly carries no `copy_num`. But `NSP1A` and `NSP1B` both did, so
every EAV genome was flagged twice for having the right architecture. This is
the failure the `copy_num` rule exists to prevent -- a feature that one genus
genuinely lacks must not be declared essential -- and it was invisible in the
coverage numbers, which only ever reported that NSP1A called 590 of 612.

**`copy_num` removed from NSP1A and NSP1B.**

### NSP7B was a real miss in EAV, and the cause was NSP7A overrunning

The third flag is not an architecture difference: EAV *does* have nsp7α and
nsp7β. The published boundaries are pp1ab **1453-1575 (nsp7α, 123 aa)** and
**1576-1677 (nsp7β, 102 aa)**.

The module called EAV **NSP7A at 133 aa** -- ten residues long -- which ate the
start of nsp7β and left a 92 aa slot between NSP7A and NSP8. **92 + 10 = 102**,
the published nsp7β length, so the arithmetic closes.

The chemistry places the same boundary without reference to the literature:
in **all ten** genomes, residue 123 of the NSP7A call is **Glu** and residue
124 is **Ala** -- `E|A`, the nsp4 3CLSP site (EAV cuts after Glu, PMID
**9971783**). Re-derived on that site: **10 NSP7A at exactly 123 aa, all
ending Glu; 10 NSP7B at exactly 102 aa, all ending Glu.**

### The other lineages that lack NSP7B are left alone, deliberately

Twenty-four genomes have NSP7A and NSP8 called with nothing between them, and
19 of the 24 slots end in Glu, so the slot is real. But outside EAV the slots
are **104-152 aa** -- LDV 140, Rat arterivirus 139, Lopma 152, bamboo rat 149 --
and several sit outside NSP7B's 90-135 window. No publication places the
nsp7α/nsp7β junction in any of these lineages, and each contributes one to five
genomes. Under the source-annotation ceiling they are not projected. They are
the reason **NSP7B's own `copy_num` is re-examined** once the EAV members are
in: a feature that cannot be called in a complete genome of a divergent lineage
must not flag that genome as defective.

## The `copy_num` set, re-derived from the calls rather than inherited (21 Nov 2026)

`copy_num` makes a feature **essential**: `viral_genome_quality.pl` scores only
features that carry it, and flags any genome that does not call one.  That is
the right behaviour for a feature every member of the module really has.  It is
the wrong behaviour for a feature a lineage genuinely lacks, or one whose
boundary no publication places outside the characterised genera -- there it
flags a complete, correctly-annotated genome as defective, and the flag is a
statement about the module's reach, not about the genome.

Twenty-two of twenty-seven features carried `copy_num`.  `qc/copynum_audit.py`
asks the only question that settles it: **is there a complete genome, routed to
this module, in which this feature is never called?**  Not "what fraction of
genomes call it" -- a feature at 97% looks healthy while being uncallable in an
entire genus.

### Two filters, both of which changed the answer

The first version of the audit had neither, and both of its headline findings
were artefacts.

**A genome in the run may have routed to a different module.**  Its features
then carry *that* module's family assignments.  Tallying on the bare feature
name credits another module's calls to this one and counts its misses against
it.  `2916858.3`, *Eospalax fontanierii baileyi* arterivirus, routed to
**Simarterivirinae** and appeared as an Arteriviridae genome calling 5 of 22
features -- which made eighteen different features look non-universal.

**The source's own completeness flag.**  A 13.1 kb deposit is not a complete
arterivirus genome merely because it clears a nucleotide floor.  Taxon
`3084900` (Wufeng rodent arterivirus 1, `OQ715740`/`OQ715741`, 13,116 and
13,141 nt) is marked `Partial` by BV-BRC, and it was the *only* evidence
against `GP2` and `GP3`.  A feature uncalled in a partial genome is precisely
the case `copy_num` exists to flag, so GP2 and GP3 keep it.  Length alone
cannot see this: both genomes clear 12,000 nt comfortably.

With both filters applied: **537 complete genomes in 20 taxa**, 74 skipped as
`Partial`, 1 as another module's.

### Two taxa the module does not cover

| taxon | genomes | features called (of 22) |
|---|---|---|
| `2499677` WPDV-AU1 | 2 | 13 |
| `1118369` Wobbly possum disease virus | 5 | 13 |

These genomes **route correctly** -- three of the twenty-five reference contigs
are WPDV (`MN508192`, `NC_026811`, `MN635448`) -- and the module calls the full
structural set in them.  What it does not call is a block:

```
called   GP2 GP3 GP5 M N  ORF1A  NSP2 NSP3 NSP4  NSP10 NSP11 NSP12   (+E, NSP5 in 1 of 7)
never    NSP1 NSP1A NSP1B NSP6 NSP7A NSP7B NSP8   GP4  ORF5A
         NSP9 ORF1AB ORF1AN ORF1ATF
```

That is the **5' half of the replicase plus two accessories**, and it is at
least partly real biology rather than a profile gap: the 5' coding region of
arterivirus genomes is highly divergent and carries different numbers of
papain-like protease domains between lineages (PMID 28053107), and ORF5a was
characterised in PRRSV, EAV and LDV (PMID 21307222, 21307223) with no
publication placing it in Zealarterivirinae.

`E`, `GP4`, `NSP5`, `NSP7A` and `ORF5A` fail **only** here, and they keep
`copy_num` deliberately.  Removing it would make these seven genomes score
Good while thirteen features are missing, which hides the gap; keeping it makes
them Poor, which is the truthful answer -- *this module does not fully annotate
this genome*.

### The six features that lost it

Every other taxon in the panel annotates to 17-22 features -- genomes the
module genuinely covers.  In those, six features are still never called:

| feature | taxa | why it is not essential |
|---|---|---|
| `NSP7B` | 11 | the nsp7alpha/nsp7beta junction is published for PRRSV and EAV and **for nothing else**; LDV, rat, Lopma, Praja and bamboo-rat slots run 104-152 aa with no paper placing the cut |
| `NSP6` | 8 | 16-22 aa.  A peptide that short cannot be profiled across this much divergence |
| `NSP8` | 4 | the short pp1a C-terminal product; the ORF1a/1b region differs in the divergent lineages |
| `NSP2` | 4 | hypervariable by design -- the largest and least conserved arterivirus protein |
| `NSP1B` | 3 | beyond EAV, three further taxa cannot call it |
| `NSP1A` | 0 beyond EAV | EAV has a single uncleaved nsp1 (PMID 7769711), which the module already models correctly as its own `NSP1` feature |

`copy_num` removed from all six.  **Sixteen features keep it**, and they are
exactly the set called at least once in every taxon the module covers:

```
E GP2 GP3 GP4 GP5 M N NSP10 NSP11 NSP12 NSP3 NSP4 NSP5 NSP7A ORF1A ORF5A
```

### The criterion, stated so it can be argued with

An earlier draft of this section drew the line at a feature count -- keep
`copy_num` where the genome annotates to 17 of 22, drop it below.  That is a
proxy, and a bad one, because the threshold is arbitrary and the same genome
can sit either side of it after an unrelated fix.  The criterion that actually
decides these cases is:

> **Is the feature absent or unplaceable in that lineage, or is the whole
> genome under-annotated by this module?**

If the feature is absent or its boundary is unpublished -- EAV's single nsp1,
the nsp7alpha/nsp7beta junction outside PRRSV and EAV, a 16 aa peptide across
40% divergence -- then flagging the genome is false, and `copy_num` comes off.
That is `NSP1A`, `NSP1B`, `NSP2`, `NSP6`, `NSP7B` and `NSP8`, failing across
three to eleven otherwise well-annotated lineages.

If instead the module simply did not reach a lineage, the flag is true and must
stay, because the quality score exists to say whether the annotation is
complete.  That is Zealarterivirinae, and it is the whole reason `E`, `GP4`,
`NSP5`, `NSP7A` and `ORF5A` keep `copy_num` despite failing in seven complete
genomes.

### Simarterivirinae, audited the same way

The sister module was put through the identical test: 54 annotated genomes, 19
skipped as `Partial`, **35 complete genomes in 8 taxa**, every one annotated to
12 or 13 of its 14 features.  **Every one of its nine `copy_num` features is
called in every complete genome routed to it**, so its essential set stands
unchanged.  Its four paralogues -- `EP`, `GP2P`, `GP3P`, `GP4P` -- correctly
carry no `copy_num`, and neither does `ORF1AB`, which is a transcript edit.

## Every cited PMID verified against PubMed (21 Nov 2026)

`PMID_claude_generated` carries **46 distinct references in Arteriviridae and
17 in Simarterivirinae, 51 distinct in all**.  Each was resolved through the
NCBI esummary API and its title read, because the project has been burned
before: on 27 Oct 2026 three PMIDs written from memory all resolved to
unrelated papers.

**All 51 resolve, and 49 are arterivirus papers on the subject they are cited
for.**  The dates span 1988 to 2021 and the titles match the claims -- the
nsp4 3CLSP specificity on `19646449`, the nsp1alpha/nsp1beta auto-cleavage on
`20006994`, the -2 frameshift on `23043113`, ORF5a on `21307222`/`21307223`,
nsp8 as the RdRp N-terminus on `30353315`.

**Two are not arterivirus papers and are flagged here rather than quietly
kept.**  `ORF1AB` cites `15680415` ("Programmed ribosomal frameshifting in
decoding the SARS-CoV genome") and `2447506` ("Characterization of ribosomal
frameshifting in HIV-1 gag-pol expression").  Both are foundational **-1 PRF
mechanism** references and neither reports arterivirus data.  They are
legitimate as mechanism citations and are kept as such; the arterivirus
ORF1ab frameshift itself rests on `20696193`, the PRRSV replicase review,
which is cited alongside them.  A reader should not take either as evidence
about this family.

## The whole cleavage cascade measured, junction by junction (21 Nov 2026)

`check_cleaved_ends.py` had never been run for this module, because it reads a
flat `ann/*.feature.tbl` tree and `evaluate_coverage.py` writes GTOs.
`qc/gto2tbl.py` bridges them.

**Its first version wrote the feature key one column early**, putting the
symbol where the start coordinate belongs.  The check then read a coordinate
as a feature name and reported `no genomes carry both` for all eleven
junctions -- a clean, confident, entirely wrong answer with no error anywhere.
The column order is the annotator's and is load-bearing: `5` is the feature
id, `6` the symbol, `7` and `8` the coordinates.

### pp1a, over 611 annotated genomes

| junction | n | exactly adjacent |
|---|---|---|
| NSP1A start vs ORF1A start | 585 | **99%** |
| **NSP1B start vs NSP1A end** | 586 | **68%** |
| NSP2 start vs NSP1B end | 580 | 99% |
| NSP3 start vs NSP2 end | 602 | 94% |
| NSP4 start vs NSP3 end | 609 | 99% |
| NSP5 start vs NSP4 end | 600 | 97% |
| NSP6 start vs NSP5 end | 588 | **100%** |
| NSP7A start vs NSP6 end | 591 | **100%** |
| NSP7B start vs NSP7A end | 571 | 98% |
| NSP8 start vs NSP7B end | 573 | 99% |
| NSP8 end vs ORF1A end | 589 | 99% |

### pp1ab

| junction | n | exactly adjacent |
|---|---|---|
| NSP10 start vs NSP9 end | 518 | **100%** |
| NSP11 start vs NSP10 end | 608 | **100%** |
| NSP12 start vs NSP11 end | 608 | 99% |
| NSP12 end vs ORF1ab end | 484 | 99% |

Every junction in the module is crisp except one, and it is the one already
under repair: **NSP1A|NSP1B at 68%**, which is the PRRSV-1 PCP&alpha; site the
collections had no member for.  The 31% that are not adjacent sit at +72 nt and
its multiples -- the signature of a boundary that is wrong by a fixed amount
rather than drifting.

**One row in that table is not a measurement and should not be read as one.**
`NSP9 start vs ORF1AB start` reports a mode of **+6,981 nt**, because the tool
assumes the first product in the list begins at the precursor's own start.
nsp9 begins roughly 7 kb into pp1ab, after the whole of pp1a, so the offset is
the architecture rather than creep.  The junction that matters for nsp9 is its
downstream one, `NSP10 start vs NSP9 end`, and that is 100%.

## Simarterivirinae: three genomes where the duplicated cassette collapses (21 Nov 2026)

`qc_duplicate_calls.py` over all 54 annotated Simarterivirinae genomes, with
the family's declared overlaps passed as `--expected`, reports **no nested
call and three exact coordinate collisions**:

| genome | accession | collision |
|---|---|---|
| `1885929.21` | `KC787623` | GP2 + GP2P, both 10659 +699 |
| `1885929.22` | `KC787628` | GP2 + GP2P, both 10659 +699 |
| `2682610.4` | `NC_048209` | GP3 + GP4 + GP4P, all 12884 +813 |

The third is the instructive one.  The cassette is otherwise laid out exactly
as it should be -- GP2P 10806, GP3P 11354, GP4P 11878, GP2 12423, GP5 13732 --
and the annotator's own trace puts GP3 at **12884-13210** and GP4 at
**13229-13693**, which is correct.  The collision happens *after* the profile
match: both features carry `upstream_ext 1` and `downstream_ext 1`, so each
call scans out to a start and a stop, and ORF3 and ORF4 overlap in this family.
All three features therefore land on the same 271 aa ORF.

**Not repaired, by instruction**: the fix would be to stop scanning to a Met
start, and that behaviour is to be left alone.  Position-based disambiguation
of the paralogues is the real answer and is a larger change than this build.
Three genomes of 54, recorded rather than patched.

### The three nested calls are architecture, not defects

`E` inside `GP2` (100% of calls), `ORF1A` inside `ORF1AB` (98%) and `ORF5A`
inside `GP5` (46%) are **overlapping reading frames**, which is how
arteriviruses encode these products -- ORF2b inside ORF2a, ORF5a inside ORF5
(PMID 21307222, 21307223), and ORF1a inside ORF1ab by definition.  The tool's
default advice for a nested call is "promiscuous profile, raise bit_cutoff or
drop"; here that would be exactly wrong, and the low blastp similarity between
the pairs (E vs GP2: 41% identity over 27% coverage, 19 bits) confirms they are
different proteins sharing genome space rather than one profile firing twice.
They are declared through `--expected` so that a future run reports them as
EXPECTED and a real nesting is not lost among them.

## One sequence in 588 put E on the ORF5a reading frame (21 Nov 2026)

`qc_duplicate_calls.py` over the 612-genome coverage run, with the family's
overlapping reading frames declared through `--expected`, leaves exactly one
finding: **`E` and `ORF5A` called on byte-identical coordinates in three
genomes** -- `28344.37242`, `28344.39812` and `28344.39814`.

E is the one in the wrong place.  ORF2b sits inside ORF2a at about 12.0 kb;
these calls are at 13368, 13414 and 13432, on the ORF5a frame inside ORF5, a
kilobase downstream.  GP5 is called correctly ten bases later in all three, so
the genome is fine.

**One profile did it in all three: `Arteriviridae.E.19`, built from one
sequence.**

| | |
|---|---|
| the member | `fig|28344.22115.peg.9` |
| deposited as | `Envelope protein` |
| length | **51 aa** -- ORF5a's length; E is ~70 |
| vs the ORF5A collection | **98.0% identity over all 51 residues** (and 90-96% to 23 more) |
| vs the E collection | nothing but itself |

It is an ORF5a filed under E.  **Re-binned into `ORF5A`** with the evidence in
its header -- moved, not deleted, so the decision can be read back.  E 588 ->
587, ORF5A 374 -> 375.

### Why the bin audit could not catch it, and what now does

`audit_bins.py`'s foreign-subset test is the right test and it has a floor:
`SUBSET_MIN = 5`, the point at which a group of foreign members is numerous
enough to survive clustering and become a profile of its own.  That floor was
measured on the M/GP4 case, where 35 foreign members made `M.4`.

**But a single sequence can become a profile too.**
`build_leftover_pssms.py --min-seqs 1` duplicates the row and builds from one
sequence, which is correct and deliberate -- a lineage-specific protein with a
single deposit must not disappear silently, and Birnaviridae `X` is exactly
that case.  The cost is that a single *mislabelled* deposit gets the same
treatment, and the resulting profile carries **no conservation information at
all**: every column is that one sequence, so it matches its true family far
better than the feature it was filed under.

`qc/singleton_contaminants.py` closes the gap: for every member of every
collection, the best hit in every *other* collection, reported at >=90%
identity over >=85% of the query, **with no count floor**.

### The sweep, validated and then run on everything

A detector that has never caught its own motivating case is not a detector.
Run against the **pre-fix** E collection, `singleton_contaminants.py` reports
exactly `fig|28344.22115.peg.9` as the one E member of 588 whose best foreign
hit is ORF5A, and 21 ORF5A members pointing back at it.  Run against the
collections as they now stand, **no member of any of the 24 matches a foreign
feature**.  Simarterivirinae is clean over all 13 of its collections --
which also confirms the duplicated cassette holds genuine paralogues rather
than duplicated bins, since GP2/GP2P and their siblings fall well below 90%
identity.

**Two design points, both of which cost a wrong answer first.**

*Exclude the queried feature from its own subject database.*  The first
version built one database of everything and blasted everything against it.
That is quadratic in the wrong term -- GP5's 2,375 members all hit each other
-- and it ran 31 minutes of CPU without finishing.  Capping `-max_target_seqs`
would have been worse than slow: for a 2,375-member collection the top hits
are all same-feature, so the single foreign hit being looked for falls off the
end and the sweep reports CLEAN.  Per-feature exclusion asks the question
directly and finishes in about two minutes.

*Require the hit to cover BOTH sequences.*  Query coverage alone reports every
mature peptide in the module, because a mature peptide is a 100%-identical
substring of its own precursor: nsp7alpha covers 149 of its own 149 residues
against ORF1A and is not a contaminant, it is the definition.  Subject coverage
separates them with no hand-maintained pair list -- **a contaminant is the same
protein as something else; a product is part of something else.**  nsp7alpha
covers 6% of ORF1A; the mislabelled ORF5a covered 100% of the ORF5a it matched.

**That heuristic has a limit, found by running the sweep on Birnaviridae.**
pVP2 is VP2 plus the C-terminal structural peptides -- 512 aa against 441 in
IBDV -- so VP2 covers **86%** of pVP2 and clears an 85% floor.  All 398 VP2
members were reported.  Lowering the floor would hide real contaminants and
raising it would not fix this pair, so near-length maturation steps are
declared explicitly (`VP2:PVP2`).  The list stays short precisely because only
near-length pairs need it.

## Coverage 11 — the final measurement (21 Nov 2026)

One rebuild pass covering four keys (`NSP7A,NSP7B,E,ORF5A`), then install and
re-measure. 612 exemplars, **612 routed, 0 unrouted**, 13,289 LowVan Annotate
calls and 1,852 transcript-edited calls.

| | coverage 9 | coverage 10 | coverage 11 |
|---|---|---|---|
| genome quality, Good | — | 542 | **569** |
| genome quality, Poor | — | 70 | **43** |
| NSP1A → NSP1B flush | 68.4% | 96.3% | 96.3% |
| NSP7A ends E | 92.0% | 92.0% | **93.5%** |

### Both repairs verified on the annotated GTOs, not on the collections

**EAV nsp7α/nsp7β.** All ten *Alphaarterivirus equid* exemplars now call
`NSP7A` at **123 aa** and `NSP7B` at **102 aa** — exactly the pp1ab 1453–1575 /
1576–1677 boundaries the literature gives. Before the repair the EAV genomes
had no NSP7B call at all. Taxon 11048 (LDV) calls NSP7A at 126 aa and still no
NSP7B; that is a residual, recorded here rather than forced.

**The E/ORF5a contaminant.** Re-checked by walking every annotated GTO for two
features on byte-identical coordinates: **0 E/ORF5A collisions**, down from 3.
E is called on 604 of 612 and ORF5A on 603, so removing the one mislabelled
deposit did not cost either feature any calls.

### The quality panel is NOT held out, and that was worth measuring

The page said the panel was held out. It is not. Blasting the 25 reference
contigs against the 612 exemplar sequences gives **25 of 25 present at 100%
identity over their full length** — both sets are drawn from the same genome
collection by the same clustering, so the references are necessarily cluster
representatives too.

The first version of this test blasted against `genomes.fna`, which is the
**pre-clustering pool** (615 sequences), not the 612-genome panel in `ex_gto`.
Same answer here, but only by luck: the pool is a superset, so a hit there does
not establish a hit in the panel. Re-run against `ex_gto`.

Scoring with the references removed:

| | n | Good | Poor | |
|---|---|---|---|---|
| whole panel | 612 | 569 | 43 | 93.0% |
| **held out** | **587** | **553** | **34** | **94.2%** |
| the references alone | 25 | 16 | 9 | 64.0% |

**The held-out figure is higher**, which is the opposite of what a panel
contaminated with training data does. The reason is that references are chosen
for *breadth*: they are disproportionately the divergent lineages — rodent
arteriviruses, wobbly possum disease virus, the simian outliers — which are
precisely the genomes this module annotates least completely. Being a reference
bought them nothing, because routing to the module was never what they failed.

Simarterivirinae: 4 of 11 references in the 24-genome panel; 20/24 Good whole,
**17/20 Good held out**.

## derive_cutoffs_separation.py is the wrong tool for a module with mature peptides

Run as a dry run (no `--write-json`) over the rebuilt profiles. It reports
**16 of 27 features as INSEPARABLE**, and 14 of those 16 collide with `ORF1A`.

That is a tautology, not a finding. **Every nsp IS a sub-sequence of ORF1A**, so
a profile built from nsp5 necessarily scores against the pp1a precursor at least
as well as against its own training set — the tool reports `NSP5 self=360
foreign=378` and concludes no threshold separates them. The remaining two,
`N` and `ORF5A` against `NSP2`, are short features measured against the longest
profile in the module.

The proposed cutoffs were correspondingly destructive (`E` 27, `N` 21, `NSP12`
18, against floors that currently work), and **none were written**. The tool's
foreign set would have to exclude a feature's own precursor before its output
means anything here.

The real question it was run to answer — whether the new **single-sequence
profiles** (12 for E, 11 for ORF5A, 8 for NSP1B) call loosely — is answered by
the coverage run instead, and answered better: E and ORF5A are called on 604 and
603 of 612 genomes with **zero coordinate collisions** between them. A loose
single-sequence profile shows up as inflation and there is none. Cutoffs in this
module are per-feature, so the single-sequence profiles are governed by a floor
that was measured on the feature's full collection.

## Two silent defects in evaluate_coverage.py's GTO reporter, both fixed

Neither raised an error. Both were found by checking a number against an
independent measurement rather than by reading it off the page.

### 1. A transcript edit's feature KEY was read as its MODULE

`family_assignments` has **two shapes**, and only one puts the module inside
`fam[1]`:

```
['LOWVAN',       'Arteriviridae.ORF1A.1', 'ORF1a polyprotein',   'LowVan Annotate dev']
['Arteriviridae','ORF1ATF',               'ORF1a-TF polyprotein','LowVan Transcript Edited Feature']
```

The reporter did `key.partition(".")` on `fam[1]` for both. For a transcript
edit that key carries no dot, so the **feature key became the module name**.
`coverage11.json` filed 563 of 612 genomes under four modules that do not
exist — `ORF1AB`, `NSP9`, `ORF1AN`, `ORF1ATF` — and left 49 under
`Arteriviridae`. `gen_facts.py` selects `rows` by `module == G`, so the page
would have reported **49 routed genomes out of 612**.

This defect only appears on a module that declares transcript edits or splice
variants, which is why ten families' pages were unaffected by it.

Fixed by branching on `fam[0]`: anything other than the literal `LOWVAN` *is*
the module. After the fix: `per_module = {Arteriviridae: 611,
Simarterivirinae: 1}`, and that single Simarterivirinae genome is the same one
the copy_num audit independently reports as routed away.

### 2. refs_used was zero on every GTO run

`scan_refs()` globs `*.feature.tbl`. **The GTO path never writes one** — the
annotator's per-genome stdout *is* the feature table, saved as
`<genome>.stdout.txt` beside the GTO. So the glob matched nothing, `refs_used`
came back `{"module": null, "own": 0, "all": 0}`, and the routing sentence
would have published *"25 reference contigs, of which 0 were used"*.

This is the second time this field has gone silently to zero by a different
route; the function's own docstring records the first. Fixed by globbing
`*.stdout.txt` as well and adding `Anno_GTOs` to the directories searched.

**Cross-checked independently** before trusting the fix: BLASTn of the 25
reference contigs against the 612-genome exemplar panel, counting distinct
best-hit subjects, gives **25 of 25** — the same number the repaired
`scan_refs` now reports. Simarterivirinae: **11 of 11** both ways.

### Re-reporting without re-annotating

`evaluate_coverage.py` has no report-only mode, and the annotation is the
expensive half — 612 genomes is hours, the report is seconds. `qc/regen_coverage_json.py`
imports the fixed reporter and points it at the GTOs already on disk. Both
coverage files were regenerated that way; the annotation was not repeated.

## The pages, and two numbers that disagreed

Ten pages, all clean through `check_artifact_html.py`, both coverage audits
carrying the required two SVGs (genome map and rarefaction curve).

| page | Arteriviridae | Simarterivirinae |
|---|---|---|
| string collapse | 420 strings &rarr; 24, 17.5&times; | 79 &rarr; 1, 79&times; |
| PSSM registry | 27 features, 417 PSSMs | 14 features, 124 PSSMs |
| vocabulary saturation | 273 vs 23 at 5,187 proteins, **11.9&times;** | 56.9 vs 10.0 at 290, **5.7&times;** |
| coverage audit | 611/612 routed, 93.0% clean | 24/24 routed, 83.3% clean |
| genome map | `28344.92`, 26 features | `38143.25`, all 14 |

### "String collapse" meant two different things on two pages

`gen_collapse.py` divides **all distinct strings** by controlled targets:
420/24 = **17.5&times;**. `gen_facts.py` divides **string-to-feature pairs** by
the same targets: 362/24 = **15.1&times;**. Both were on the same module's
pages, under the same words, with no indication they were different
measurements. The gap is exactly the 7 strings that bind to more than one
feature. Both numbers are now named and defined in section 09 rather than left
to be read as a contradiction.

### A pair count reported as a string count, and published

`gen_collapse.py` printed *"16 string(s) bind to more than one feature"*. The
expression behind it is `nbound - len(bound_strings)`, which counts **extra
pairs**: a string bound to five features contributes four. The real answer is
**7 strings, 16 extra pairs**, and the 16 had already been written into the
page prose as a string count. Fixed in the generator, which now reports both.

The seven are worth listing, because they are the whole argument for the genus
column in every binning rule:

| string | features it is deposited on |
|---|---|
| glycosylated membrane protein | GP2, GP3, GP4, GP5, ORF5A |
| hypothetical protein | GP2, GP3, GP4, GP5, ORF1A |
| glycoprotein | GP2, GP3, GP5 |
| structural glycoprotein | GP2, GP3, GP5 |
| glycosylated minor envelope protein | GP2, GP3, GP4 |
| envelope glycoprotein | GP2, GP5 |
| viral protein | GP3, GP4 |

### The panel table was rebuilt rather than re-measured

`gen_facts.py` wants `panel_eval.tsv` from `run_gto_eval.py`, which
**re-annotates** to produce it. The coverage run had already annotated and
scored this exact panel, so `qc/panel_eval_from_quality.py` reads the six
columns out of `Quality_GTOs/` instead. Cross-check: it reports 569 genomes
with no flag of any kind, and the quality model independently reports 569 Good.

## Status at close -- 21 November 2026

Complete. Nothing is mid-flight; no background job is running.

| | |
|---|---|
| modules | `Arteriviridae` (27 features, 417 PSSMs, 16 essential), `Simarterivirinae` (14 / 124 / 9) |
| coverage | 612 exemplars, 612 routed, 0 unrouted, 569 Good / 43 Poor; simian 24 / 24 / 0, 20 Good |
| held out | 553 of 587 = 94.2%, higher than the 93.0% whole-panel number because only 16 of the 25 references score Good |
| repo | `jimdavis1/LowVan_Skill` `main` at `bf1b7ae5`, push verified by remote ref |
| pages | 12 live artifacts (4 each for Arteriviridae, Simarterivirinae, Birnaviridae), indexed in `ARTIFACTS.md` |

Build tree is `~/lowvan_scratch/arteri` and is NOT in the repo; it holds the
collections, alignments, eval runs and the `qc/` scripts. The shipped copy of
each module is `kit-verify/modules/<Module>/`.

Three scripts written during this build live in `qc/` and are also copied into
the kit: `regen_coverage_json.py` (re-derive a coverage json from GTOs on disk
without re-annotating), `panel_eval_from_quality.py`, `add_genome_map_figure.py`.

Two untracked backups are left in the kit on purpose:
`evaluate/evaluate_coverage.py.bak_seg` and `skill/scripts/gen_facts.py.bak_shadow`.

**The one open question, for whenever it is wanted:** the module-derivation
bug fixed in `evaluate_coverage.py` was silent on every module declaring
`special: transcript_edit` or `special: splice`, so those modules' published
coverage numbers understate routing. Arteriviridae is the first module here
with four such features, which is why it surfaced now. `regen_coverage_json.py`
answers it per module in minutes, without re-annotating.

**Deliberately unrepaired**, all argued above: the three Simarterivirinae
cassette collisions (`scan_to_met_start` left alone), `ORF1ATF`/`ORF1AN` at
~70%, Zealarterivirinae keeping `copy_num` so its gap stays visible, and
`28344.39812` dropped as an outlier.
