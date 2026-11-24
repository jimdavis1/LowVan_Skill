# Hepadnaviridae — LowVan build notes

Built November 2026.  Two modules, **Orthohepadnavirus** (7 features) and
**Avihepadnavirus** (5 features).  This is the first circular-genome family in
LowVan and it is the test case the circular-replicon support in kit commit
`63f0b1ce` was written for.

---

## 1. Why this family could not be built before

A circular genome reaches GenBank as the loop cut at one arbitrary point, so a
CDS lying across that cut arrives split: part at the end of the contig, part at
the start.  Searched as deposited it is two partial hits, neither of which
clears its `bit_cutoff`, and the feature is lost.

For hepatitis B that is not an edge case, it is the normal state of the data.
Measured on 400 near-full-length deposits (17 Oct 2026):

| protein | crosses the origin |
|---|---|
| P (polymerase) | **97.2%** |
| L (large surface) | **100%** |
| X | 4.2% |
| S, preC, C | never |

And there is no better place to cut.  Coding density is **98.1%** — 3,141 of
3,203 nt lie inside at least one CDS — the polymerase alone covers 78.7% of the
genome, and the longest non-coding window is **62 nt**, so a break avoids every
ORF only 1.9% of the time.  The field's own convention does not help: position 1
is nominally the EcoRI site, which sits inside both P and the envelope ORF, and
96.4% of deposits break inside P.

The reference genome says the same thing in its own annotation.  RefSeq
**NC_003977.2** (HBV genotype D, 3182 bp, circular) writes three of its seven
CDS as joins:

    P     join(2309..3182,1..1625)   832 aa
    L     join(2850..3182,1..837)    389 aa
    M     join(3174..3182,1..837)    281 aa
    S     157..837                   226 aa
    X     1376..1840                 154 aa
    preC  1816..2454                 212 aa
    C     1903..2454                 183 aa

Duck hepatitis B virus, RefSeq **NC_001344.1** (3027 bp, circular), wraps a
*different* pair — preC `join(2524..3027,1..414)` and C `join(2653..3027,1..414)`
— which is the point: the origin is wherever the submitter cut the circle, and
nothing about the biology puts it in the same place twice.

Genome references: Galibert *et al.* 1979, PMID **399327** (HBV, subtype ayw);
Mandart, Rieger & Galibert 1984, PMID **6699938** (DHBV); Galibert *et al.*
1982, PMID **7086958** (woodchuck).

---

## 2. Two modules, because the genome organisation differs

The partitioning rule is organisation, not phylogeny, and the deposits settle it
before any reference is consulted.  Across 636 unique Avihepadnavirus protein
sequences there is **exactly one** `X protein` string and **five** `middle S
protein` strings, against 11,393 and 7,390 in Orthohepadnavirus.  All six are
mislabelled.  Avihepadnaviruses carry neither the X gene nor the middle
(preS2+S) surface form, which NC_001344.1 confirms by annotating neither.

| module | genomes | features |
|---|---|---|
| Orthohepadnavirus | 142,997 | P, L, M, S, preC, C, X |
| Avihepadnavirus | 352 | P, L, S, preC, C |

Both clear the 25-reference budget comfortably: Orthohepadnavirus routes
**99.9%** of its 18,021 in-window contigs on 25 references, Avihepadnavirus
**100%** on 2.

### Intake

Dropped before any collection was built:

| dropped | n | why |
|---|---|---|
| `... sp.` records | 17 | no assigned species; LowVan does not build modules for unclassified viruses |
| named only `Hepadnaviridae <code>` | 7 | a family-rank label on a metagenomic assembly, not a species |
| fish/amphibian genera | 64 | tabled, see §3 |
| other unassigned fish/novel | 6 | tabled with the above |
| *Alewife hepatitis B virus* | 2 | a **fish** hepadnavirus that BV-BRC files under genus Orthohepadnavirus. It had been selected as `Orthohepadnavirus.23.dna`; a reference contig for a taxon the module cannot call is worse than none, because it converts a visible routing gap into a silent miscall |

Kept, and folded into Orthohepadnavirus because they carry X and the three
surface forms: bat hepatitis virus, bat hepadnavirus, ringtail hepadnavirus and
Tai Forest hepadnavirus, which BV-BRC leaves without a genus.

---

## 3. What is NOT built, and why

**Metahepadnavirus, Parahepadnavirus and Herpetohepadnavirus are tabled.**
Three genera, 64 genomes, and **80 unique protein sequences between them** —
fish and amphibian viruses that are deeply divergent from each other.  Per
feature that is 7 to 31 sequences spread across three genera, with ORFY and
ORFZ present as single deposits.  Routing is not the obstacle (33 contigs fall
into 4 clusters and 25 references would cover them twice over); the obstacle is
that there is not enough protein diversity per genus to build profiles that
would cover any of the three, and sharing profiles across genera that distant
would produce a module whose cutoffs cannot separate them.  Recorded as a
documented gap rather than forced.

**HBSP, the hepatitis B spliced protein, is declared and not built.**  It is
real: PMID **12586301**, "The expression of hepatitis B spliced protein (HBSP)
encoded by a spliced hepatitis B virus RNA is associated with viral replication
and liver fibrosis" (*J Hepatol* 2003).  Four deposits in this dump carry
`spliced protein` at exactly **111 aa**, which is HBSP's published length, and
four more carry `splice-generated protein` at 68, 76, 96 and 364 aa.  Eight
sequences is below the member floor, and building it properly means a
`special: splice` feature with literature-anchored donor and acceptor
coordinates on a genome that also wraps — the one interaction the new circular
support has not been exercised on.  Named here so it is a known gap and not an
invisible one.

---

## 4. The nested forms, which are the real hazard in this family

L, M and S are three in-frame starts on **one** surface ORF sharing **one**
stop, so S is an exact C-terminal substring of M and M of L.  preC and C are the
same arrangement on the core ORF.  This is the VP1/VP2 trap, and in this family
it has two separate consequences.

### 4a. Contamination between the forms

`fasta-cluster-pssm` drops a cluster whose consensus is a substring of
another's, so a mixed collection does not merely blur the forms — it deletes the
shorter one and ships its profile as the longer.  Measured against the
literature size classes before any profile was built:

| named | L-size | M-size | S-size | other |
|---|---|---|---|---|
| L | 9366 | 127 | **484** | 6528 |
| M | 0 | 5647 | **502** | 2182 |
| S | **47** | 10 | 9814 | 6889 |

About 1–5% cross-contamination, an order of magnitude better than
Parvoviridae's 66-of-89 but still enough to matter, and removed by the length
windows.

### 4b. The generic surface strings, which no name rule or homology can place

6,776 deposits carry a string that names the ORF and not the form — `surface
antigen`, `surface protein`, `HBsAg`, `envelope protein`, `S antigen`, `sAg`.
Homology cannot resolve them, because the forms are 100% identical over their
shared region: a 226 aa query matches the L canonical over 100% of itself
exactly as well as it matches S, so the rescue's ≥80%/≥60% test passes on all
three and the winner is arbitrary.

Length alone cannot resolve them either, and that is the part worth stating.  A
226 aa surface deposit is either a genuine S **or an L truncated after 226
residues**, and those are opposite things: one belongs in the S collection, the
other is a partial that must be discarded.  Binning on length alone builds the
S profile out of L's N-terminus.

`resolve_surface_forms.py` applies the structural test instead — **which end
does the sequence share?**  Every genuine form ends at the common stop, so every
one of them carries the S C-terminus; a truncated L carries L's N-terminus and
stops wherever the submission stopped.  The anchor is taken from the module's
own S collection (`GPSLYNILSPFLPLLPIFFCLWVYI` for Orthohepadnavirus) rather than
written down, so it stays right for the woodchuck and bat viruses.

| verdict | n |
|---|---|
| assigned → L / M / S | 947 / 65 / 1439 |
| discarded: does not share the C-terminus (truncated parent) | 18,564 |
| discarded: non-standard residue | 1,889 |
| discarded: shares the C-terminus but matches no form window | 1,767 |

The last row is 32–359 aa — N-terminally truncated S, plus a real population of
preS-deletion variants at 319–340 aa, which are genuine but defective and do not
belong in a profile's training set.

The same pass moved 28 Avihepadnavirus deposits filed as `core protein` that are
305 aa — those are precore — and 8 filed as surface strings that are the large
form.

**`build_collections.py` truncates the collection files, so re-running it alone
silently discards all of this.**  `build.sh` chains the two for that reason.

---

## 5. Windows, set per species

Orthohepadnavirus is not only HBV.  Woodchuck P is 883–884 aa against HBV's
832–845, its large surface protein is 431 aa against HBV's 389–400, and Tai
Forest core is 220 aa against HBV's 183.  A window cut to HBV alone deletes
every non-human orthohepadnavirus silently.

| feature | window | literature anchors |
|---|---|---|
| P | 790–900 | HBV 832/843/845, WHV 883/884, GSHV 877, WMHBV 835, TFHV 819 |
| L | 360–445 | HBV 389/400, WHV 431 |
| M | 255–300 | HBV 281, ringtail 271 |
| S | 210–250 | HBV 226, WHV 242, WMHBV 228, ringtail 225 |
| preC | 195–235 | HBV 212, WHV 224/225, ringtail 219 |
| C | 170–225 | HBV 183/185, WHV 188, GSHV 187, ringtail 191, TFHV 220 |
| X | 130–165 | HBV 154, WHV 141/142, GSHV 138, WMHBV 152, ringtail 139 |

P's window discards the most, and correctly: the median `polymerase` deposit in
this genus is **344 aa**, because drug-resistance surveillance sequences the
reverse-transcriptase region only.  3,557 deposits sit at exactly 344 aa.  They
are real partial polymerase and they are not the protein.

---

## 6. `upstream_ext` is 0 on M, S and C, and that is load-bearing

`scan_to_met_start` begins at `$start - 3` — one codon **upstream** of the
call — and walks left in 3-nt steps until it meets a Met or a stop.  Between
S's start and M's, and between M's and L's, there is no stop: they are the same
ORF.  So `upstream_ext: 1` on S does not refine S's start, it converts every S
call into M and then into L.  The same holds for C against preC.

    L        upstream_ext 1    longest surface form, nothing upstream in frame
    M, S     upstream_ext 0    a longer sibling start lies upstream, no stop between
    preC     upstream_ext 1    longest core form
    C        upstream_ext 0    precore's start lies upstream in the same frame
    P        upstream_ext 1    its own ORF
    X        upstream_ext 0    see section 9 -- HBwx's start lies 56 codons
                               upstream in the same frame, with no stop between
                               on the genomes that carry it
    HBwx     upstream_ext 0    longest form of the X ORF; nothing above it

`downstream_ext` is 1 on all of them: every feature ends at a genuine stop
codon, none is a protease product, and on a circular contig the scan runs
through the origin because the search is done on the doubled sequence.

---

## 7. Segment windows are the circular guard

See `SEGMENT_WINDOWS.md`.  In short: on a circular module `min_len`/`max_len`
are not only the quality check — the annotator consults them before doubling a
contig and refuses anything outside them, because a partial contig is not a
circle.  Orthohepadnavirus 3130–3400, Avihepadnavirus 2980–3120.  The lower
bound costs the wrap recovery on 992 near-complete HBV records (5.2%), which is
the intended trade; the upper bound is what refuses the head-to-tail 1.3-mer and
2-mer constructs HBV is commonly deposited as.

---

## 8. HBSP: spliced and origin-crossing at the same time

HBSP is `special: "splice"`, so it has no PSSM and is called from a reference
junction in `Splice-Variants/Orthohepadnavirus/HBSP.fasta`.  It is 93 aa: the
first 46 residues of the polymerase, then 47 residues of novel sequence in a
shifted frame.

The coordinates in circulation for it — donor 2447, acceptor 489 — do **not**
work on NC_003977.2.  They land on `AT` and `CC`, not `GT` and `AG`, because
they are in EcoRI numbering; this is the same +2 offset that puts HM153811's
P at 2307 against 2309 here.  Soussan 2000 (PMID 10619861) gives no nucleotide
coordinates at all and names no reference genome.  What it gives is the
product: *"This protein, with 93 amino acids, begins at the polymerase start
codon and includes the first 46 amino acids of the polymerase, followed by the
47 amino acids corresponding to a new HBV amino acid sequence generated by the
frameshift."*

So the junction was located by splice chemistry on NC_003977.2 and then checked
four ways, all of which had to agree:

    donor 2449      followed by GT
    acceptor 545    preceded by AG
    product         exactly 93 aa
    partition       46 + 47, first 46 byte-identical to unspliced P

giving `join(2309..2449, 545..685)`.

The reference junction spans the origin.  That works only because the contig is
still doubled when `get_splice_variant_features.pl` sees it — folding happens
afterwards, in `fold_circular_features.pl`.  HBSP is the first feature anywhere
in the kit that is both spliced and origin-crossing, so that ordering had never
been exercised before.

---

## 9. Whole-X (HBwx) is invisible from the reference genome

HBwx is 210 aa: an in-frame ATG 56 codons upstream of X's, reading through to
X's stop.  Pre-X 56 + X 154, a 630 bp gene.

The reference cannot make it.  NC_003977.2 is genotype D and carries **two**
in-frame stops in the pre-X region, so translating those 56 codons gives
`LGLGHGPSAHAWNLFGSSADPYCGTPSRLFCSQQVWSKHYRD**LCCPIPQIYIVS`.  A module built
by reading the reference and its RefSeq annotation would never have seen this
protein at all.  It was found by scanning every complete deposit for an ORF
that ends at X's stop and is longer than X:

    complete contigs scanned (3.1-3.4 kb)   18,021
    carrying X                               8,609
    with an open pre-X frame                 1,178   (13.7%)
    of those, exactly 210 aa                 1,170

which sits in the range the literature reports — 21% of healthy carriers and
64% of chronic liver disease patients in the pre-X survey (PMID 15948258).
HBwx and X are functionally distinct proteins (PMID 27255760).

Two parameters are deliberate and neither is the module default:

  * **no `copy_num`.**  Absence is biology here, not an annotation failure.  A
    `copy_num` of 1 would flag 86% of HBV genomes as defective.
  * **`coverage_cutoff` 0.85**, against 0.65 everywhere else in the module.  A
    154 aa X covers 73% of a 210-column profile.  At 0.65 the HBwx profile
    would fire on every X-only genome in the database.

And X itself moves to `upstream_ext: 0`.  With extension on, X's Met-walk
climbs out of X and into the pre-X region on precisely the genomes where that
region is open — the walk stops at the first Met or stop, and on an HBwx-positive
genome the first thing it meets going left is HBwx's start codon.  X would
silently become HBwx and the real X call would be lost.  This is the same trap
documented for M, S and C in section 6, with one more instance.

## 10. Pre-pre-S (whole-L): rejected on a search for the wrong word

**EVERY PMID IN THIS SECTION IS CLAUDE-GENERATED AND UNVERIFIED.**  None of it
may be written into a `PMID` field; if any of this is built, the citations go in
`PMID_claude_generated` and stay there until Jim has checked them.

42 Orthohepadnavirus deposits (35 unique sequences) carry the annotation
`putative prepreS protein` or `putative prepre-S protein`.  These were tabled
on 23 Nov on the ground that no experimental report could be found.  That was
wrong, and the reason is worth recording: the search was for *preS*, and the
term in the literature is **pre-pre-S**.  Exactly two PubMed records use it.

The substantive one is a review (PMID 16364485, Faure, *Virus Res* 2006;
117:185–201), which states that in more than 300 Asian genomes a **pre-X** or
**pre-pre-S** region is found just upstream of, respectively, the X and the
pre-S1 genes; that both are translated in frame with their own gene; that the
corresponding proteins have a demonstrated transactivating role; and that both
are genotype- and serotype-specific (C/adr).  That is one sentence covering
pre-X and pre-pre-S together — the same sentence that licensed HBwx in section
9.  Building one and rejecting the other was not a defensible split.

The deposits agree with the paper without having been asked to.  Measured
across all 35 unique sequences, 31 carry an in-frame N-terminal extension of
**exactly 45 residues** ahead of the L start:

```
MQLIITSKLGIIYILCGRLAFYIREKLHAAPHFVGHHILGNKSYS | MGGWSSKPRQGM...
└────────────── 45 aa = 135 nt ──────────────┘   └─── L (preS1) ───
```

A secondary source (PMC4504892 / PMID 15991289) attributes to Dong *et al.* an
ORF upstream of preS1, recovered by long-accurate PCR from serum-derived
genomes, **"about 135 bp, provisionally called pre-pre-S"**, with promoter
activity.  135 nt was measured here from the deposits before that number was
read anywhere, so the agreement is independent.  The Dong primary has been seen
only through that citation chain and has not been read.

The remaining 4 of 35 are the same architecture with L's own initiator mutated
away — `...KSYS·KPRQGM`, 439–440 aa against 445.  That is the HBwx failure mode
exactly (section 9): a nested isoform whose only distinguishing feature is an
upstream ATG.  So if this is built as `WHOLEL`, it takes `require_met_start: 1`
and L takes `upstream_ext: 0`, by the same argument as HBwx and X.

**Not built.**  The review clears the "seen experimentally" bar for the
transactivation claim; the primary has not been read; and the decision to add a
ninth-and-tenth nested surface form is Jim's, not mine.

## 11. p25 / p22 / p17 are mature peptides on PRECORE, and the cut is genotype-specific

**CLAUDE-GENERATED, UNVERIFIED.**  Same rule as section 10.

The precore product is processed twice:

| | what it is | made by |
|---|---|---|
| **p25** | the precore primary translation product | the precore AUG, 29 codons above C's |
| **p22** | p25 minus the 19 aa ER signal | signal peptidase, in the ER |
| **p17** | secreted HBeAg | **furin**, in the Golgi, removing the C-terminal arginine-rich domain |

The furin residue blocked this for a week because the reviews disagreed —
Arg-154 against Arg-167.  They are not two readings of one fact.  **The cut is
genotype-dependent** (PMID 19193799, Ito, Kim, Lok & Tong, *J Virol* 2009;
83(8):3507–17, which identifies furin using furin-deficient LoVo cells and a
convertase inhibitor):

* non-A genotypes cut after the **first** RXXR motif, `151-RRGR-154` → C-terminus at **Arg-154**
* genotype A carries a two-residue insertion (`151-RRDR-154`) and shifts to the
  **second** (`156-RSPR-159`) or **third** (`166-RRRR-169`) motif

**This is why they cannot ship as a fixed-coordinate mat_peptide trio.**  p25
and p22 have genotype-independent boundaries and could be built today; p17's
C-terminus cannot be written as one number.  Either p17 is built per genotype,
or it is built to the non-A site with the genotype-A error stated on the page,
or it is not built.  That is a decision, not a measurement, and it is Jim's.

One hard constraint if any of these is built: a mature peptide does **not**
begin with Met, so `require_met_start` must be absent — setting it would delete
the feature entirely.

## 12. Protected bin pairs, and why each one is protected

`audit_bins.py` asks of every collection whether it holds one protein and its
own protein. In a family built almost entirely out of **nested in-frame
isoforms sharing a stop**, that question has a false answer by construction:
every S really is 100% identical to part of its own genome's M, because it
*is* part of it. The audit cannot tell that from contamination, so the pairs
are declared in `collections/Orthohepadnavirus/protected_pairs.json` and
excluded from judgement.

Protection suppresses findings, so each pair needs a reason that does not come
from the audit itself.

### The envelope family: L / M / S, and WHOLEL

One ORF, three in-frame starts, one stop. L = preS1 + preS2 + S;
M = preS2 + S; S alone. WHOLEL adds the pre-pre-S extension upstream of L.

The independent evidence that no sequence is misfiled is **length**, measured
24 Nov 2026 over the built collections:

    S   n=10879   210 - 249 aa
    M   n= 6007   255 - 299 aa
    L   n=11082   360 - 431 aa

Three disjoint windows with no overlap at any point. A misfiled M in S would
be 255+ aa in a bin whose largest member is 249. The length audit is therefore
doing the work that the bin audit cannot, and it is clean.

(WHOLEL was split out of L on this same date: 26 whole-L sequences, 22 of them
at exactly 445 aa, were sitting in L, and L's upper bound had been widened to
445 to accommodate them. Removed, L's max returns to 431, which is the
woodchuck protein.)

### The core family: C / PRECORE, and the mature peptides

preC is C plus the 29-residue precore peptide, same stop. p25 is the PRECORE
CDS primary product; p22 is p25 minus the 19-residue signal peptide; p17
(HBeAg) is p22 minus the C-terminal arginine-rich domain. Every one of these
is a sub-sequence of the one above it **by definition of the processing
pathway**, which is the thing being annotated.

Here the length windows do overlap — C runs 170-220 and PRECORE 195-234 —
because the module spans several host species and a 220 aa C belongs to a
genome whose preC is 249. The overlap is across species, not within a genome.

Note that the audit's `KEPT` row for C is correct and should stay readable:
C and PRECORE are a mutual-duplicate pair, and retiring both would delete the
protein. Protection is the right instrument precisely because neither bin is
wrong.

### What protection does NOT cover

A pair not listed is still judged. Adding one is a verdict-level change and,
since 24 Nov 2026, costs nothing to re-run: `audit_bins.py --cache` stores the
measurement, keyed on the content of every collection plus the measure
thresholds, so editing the protections re-judges in seconds while editing a
collection correctly invalidates it.
