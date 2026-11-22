# Nodaviridae — build notes

Built 21–22 November 2026 as **four modules**, one per genome organisation.
Working tree `~/lowvan_scratch/noda`, test repo `VA_noda`.

Every protein, cleavage product and frameshift product declared here is backed
by a paper whose PMID was resolved against PubMed and whose title was read
against the claim. The full survey, with one table per lineage and an
evidence level on every row, is `LITERATURE_SURVEY.md`.

---

## 1. The partition, and why it is four modules

ICTV MSL41 (2025) ratifies **two** genera in this family, Alphanodavirus and
Betanodavirus (report PMID 30431412). The dump holds four distinct genome
organisations, and a module declares a protein *set*, so each one needs its
own:

| module | RNA1 | RNA2 | genomes |
|---|---|---|---|
| `Alphanodavirus` | A + B1 + B2 | alpha, **cleaved** to beta + gamma | 84 |
| `Betanodavirus` | A + B1 + B2 | CP, **not cleaved** | 1,456 |
| `Macrobrachium_rosenbergii_nodavirus` | A + B2, **no B1** | CP, not cleaved | 153 |
| `Le_Blanc_nodavirus` | A only, **no B1 or B2** | CP + CP-delta (**−1 frameshift**) + delta | 101 |

Alphanodavirus and Betanodavirus carry the same gene *count* and still cannot
share a module: the alphanodavirus capsid cleaves itself and the betanodavirus
capsid does not (PMID 9672600). A combined module would declare beta and gamma
for 1,456 genomes that never make them.

**Naming.** The crustacean and nematode groups have no genus at any rank —
NCBI files all seven species under "unclassified Nodaviridae" — so each module
is named for its most-sampled species, the convention already used for
`Orthopneumovirus_muris` (curator's rule, 22 Nov 2026). Membership is recorded
in each module's docstring because the name understates it. The margin is thin
for the nematode module and worth stating: Le Blanc 49 genomes, Santeuil 44,
Orsay 4, Melnik 4 — **Orsay virus carries nearly all the literature** but a
fraction of the sequence.

*Gammanodavirus* was not used: PMID 30346944 asks "A new genus within the
Nodaviridae?" in its title and proposes no ratified name.

### The routing budget, measured before building

25 references per module is the project ceiling. Every module reaches **100%**
of its own genomes well inside it:

| module | references used | genomes routed |
|---|---|---|
| Alphanodavirus | 10 | 31/31 |
| Betanodavirus | 3 | 603/603 |
| Macrobrachium_rosenbergii_nodavirus | 5 | 63/63 |
| Le_Blanc_nodavirus | 7 | 101/101 |

Re-measured independently at the annotator's own BLASTn settings
(`-evalue 0.5 -reward 2 -penalty -3 -word_size 11 -soft_masking false`, floor
150 bits) over all 788 in-scope genomes ≥800 nt: **788 routed, 0 unrouted,
0 misrouted.** `check_rep_contigs.py` was not used for this number because it
calls bare `blastn`, which is megablast and under-reports below ~85% identity.

### No other module in the repository can capture a nodavirus

The skill warns against a module-only routing repo because it absorbs
misrouting genomes. That was tested rather than assumed: every in-scope genome
was BLASTed against **all 1,466 reference files in the runtime**, at the
annotator's settings. The best hit to any non-nodavirus module is **45 bits**
(Kobuvirus), against a routing floor of 150. There are no neighbours, so a
Nodaviridae-only repo is the family repo.

---

## 2. What was dropped at intake, and why

416 records, each logged with its reason in `SCOPE_DECISIONS.tsv`.

**413 catch-all records.** `Betanodavirus sp.` (taxid 1906674, 224 records),
`Nodaviridae sp.` (189) and the rest carry no assigned species, and LowVan does
not build modules for unclassified viruses. Measured before dropping, so the
cost is on the record: the 223 proteins under `Betanodavirus sp.` are **97.6%
mean identity** to a named betanodavirus, so this removes no diversity — it is
a scope rule, not a quality one. That is the opposite of the Ambidensovirus
`sp.` bin the skill documents, which held 27–38% identity between classes.

**The `sp.` test needed a second condition here.** The skill's rule drops an
organism name matching `<taxon> sp.`, and in this family that regex matches two
unrelated things:

```
Betanodavirus sp. SFKS08-2012TYOy2    the VIRUS has no species      -> OUT
Mullus sp. betanodavirus 312.1.2005   the HOST FISH has no species  -> IN
```

Nineteen records across six host genera — *Mullus*, *Gobius*, *Epinephelus*,
*Pagellus*, *Balistapus*, *Pseudorhombus* — are named viruses from fish whose
species was not determined. Testing the whole string would have deleted them.
`assign_lineage.py` applies the test only to the part of the name before the
virus word.

**2 records in the wrong genus.** BV-BRC files *Wuhan nodavirus* (taxid 356127)
and *Mosquito nodavirus MNV-1* (690897) under genus Betanodavirus. Both are
insect noda-like viruses. A rep contig for either would route an insect virus
into the fish module, which has no profile that can call it — the skill's rule
that a reference for an uncovered taxon is worse than no reference, because it
turns a visible gap into a silent one.

**1 record in the wrong family.** `fig|222557.10.peg.1`, deposited under the
MrNV taxon id with the product string `B2 protein`, is 174 aa on an 872 nt
contig (AY655724) ending at nt 610. It matches neither MrNV B2 nor MrNV capsid
at e<1e-3, and its arginine-rich N-terminus is a capsid RNA-binding arm. It is
**extra small virus capsid CP17** — XSV NC_043492.1 is 872 nt with CP17 at
174 aa, CDS 89..613. XSV is a satellite virus in *Sarthroviridae*, genus
*Macronovirus* (PMIDs 14993649, 30507372), so it is outside this family
altogether. It had passed the first length window and was found only because
the length audit flagged it at 174 aa against a 131–133 aa literature range.

**10 records folded IN.** *Atlantic cod nodavirus* has no genus in BV-BRC and
is 88.6% identical to the Betanodavirus canonical polymerase over its full
length, against 50–61% for every other species in the untabled tail. It is the
BFNNV genotype, and BV-BRC already files the near-identical *Atlantic cod
betanodavirus* records under genus Betanodavirus — an inconsistent assignment
on one species, not a new lineage.

### What is NOT covered, and will not be

554 records of metagenomic noda-like viruses (~300 species, most of them
singletons) are **out of scope and tabled**. Measured at the 25-reference
budget they reach **37%**, far below the 60% that would make them marginal.
They are genuinely divergent: the best protein identity any of them reaches to
a module's canonical polymerase is 50–61%, which is between-genus distance.
Routing confirms the separation — 396 of them hit no module reference above
150 bits at all.

77 of them *do* route to a module. That is not misrouting for most: the group
includes species BV-BRC simply left without a genus. The rest are divergent
relatives that a module will call partially, and a reviewer should read the
routing figures as covering the four modules' own genomes, which is what the
788/788 measurement above reports.

---

## 3. The proteins, and where each boundary comes from

### Alphanodavirus — the capsid cuts itself, and the site was measured

The literature gives FHV's maturation cleavage as **Asn363 ‖ Ala364**, yielding
beta (363 aa) and gamma (44 aa) — PMID 1404613 ("provirions … mature by
autocatalytic cleavage of the 407-residue coat precursor protein between
asparagine residue 363 and alanine residue 364"; N363D/T/A abolish cleavage and
infectivity, and all nine revertants restored Asn) and PMID 9765417.

**The residue number is FHV's, not the genus's.** Alpha is 395–407 aa across
the genus, so the cut was placed by *chemistry* on each genome's own precursor:
scan for Asn followed by Ala or Ser with 38–52 residues remaining. Measured
22 Nov 2026 over the distinct alphas in this dump, that finds **exactly one
site in every species and two in none**:

```
flockense (FHV)     407   N363|A   gamma 44      heteronychi (BBV)  407  N363|A  44
American nodavirus  407   N363|A   gamma 44      HB-2007/CHN        407  N363|A  44
boolarraense        403   N359|A   gamma 44      nodamuraense (NoV) 399  N355|A  44
Newington virus     395   N351|A   gamma 44      pariacotoense      401  N361|S  40
```

Pariacoto is the one documented exception and it came out of the measurement
unprompted: PMID 10799587 states "the amino acid sequence at the cleavage site,
which is Asn/Ala in all other alphanodaviruses, is Asn/Ser in PaV". The
derivation reproduces FHV's published numbers without being given them.

Verified end to end on NC_004144: alpha 23–1246, beta 23–1111 (**363 aa**),
gamma 1112–1246 (**44 aa**) — flush, no gap, no overlap, and beta+gamma tile
alpha exactly. Across the whole panel, beta starts at alpha's start 4/4, gamma
ends at alpha's stop 4/4, gamma starts one base after beta ends 4/4.

### B1 — real, and its start codon is unpublished

B1 is the C-terminal fragment of protein A translated from an internal AUG on
subgenomic RNA3, *not* a cleavage product (PMID 15163762: "Protein B1 is in the
same reading frame as protein A and thus represents its C-terminal fragment").
It is therefore a CDS.

BV-BRC holds **3 deposits in Alphanodavirus and 2 across 1,456 Betanodavirus
genomes**, and **no RefSeq annotates a B1 CDS anywhere in the family**. The
collections were built by extraction from the genomes under the one property
the literature does fix — in frame with A, sharing its stop — so an accepted
member must be an **exact C-terminal substring of protein A from its own
genome**. 13/13 passed in Alphanodavirus and 124/124 in Betanodavirus.

**Two defects were found here and both are recorded in the JSON.**

*The extraction let the start roam.* At a Met slack of 40 the extractor took
the nearest Met rather than the aligned one and produced members of 102, 112,
131 and 138 aa from one gene. Requiring the aligned position itself to be a Met
returns members at exactly **102 and 131 aa** — the two documented lengths —
and rejects 8 genomes whose B1 start is not placeable from the literature. A
thin collection of correct starts beats a full one of invented ones.

*`upstream_ext` had to be 0.* B1 lies entirely inside the protein A reading
frame, so there is no stop codon upstream for `scan_to_met_start` to halt on.
With `upstream_ext: 1` the scan walked left through protein A to an arbitrary
Met: on FHV NC_004146 it called B1 at **2620..3036 = 138 aa** where GenBank
annotates **2728..3036 = 102 aa**. This is the Norovirus mature-peptide trap in
a CDS. With it set to 0, FHV now calls A 40–3036, B1 2728–3036 and B2
2738–3058 — all three exactly as GenBank annotates them.

### Betanodavirus — the capsid is NOT cleaved

No beta and no gamma are declared, deliberately. PMID 9672600 measured it: the
RNA2 ORF "carries the catalytic D-75 residue but lacks the site for
autocatalytic proteolysis", and the same paper excluded alternative initiation,
a −1 frameshift and precursor degradation as explanations for the 43/45 kDa
doublet. The T=3 crystal structure is built from the intact 338-residue chain
(PMID 26491970). The probe for beta/gamma in `LITERATURE_GAPS.md` returns
ABSENT, which is the required result: it proves the alphanodavirus rule did not
leak across.

### Macrobrachium_rosenbergii_nodavirus — no function is claimed for B2

RNA silencing suppression is demonstrated for alphanodavirus B2 (PMID 12016316)
and betanodavirus B2 (PMID 16352533) and **has not been shown for MrNV, PvNV or
CMNV**. This module's B2 is annotated `Non-structural protein B2`; the other
two say `RNA silencing suppressor protein`. No subgenomic RNA is claimed
either, because none has been demonstrated for these three.

The one reported MrNV capsid cleavage, at Arg26, is an *E. coli* cysteine
protease artefact of recombinant expression, blocked by E-64 (PMID 34445426).
It is not modelled.

### Le_Blanc_nodavirus — a −1 frameshift and a non-AUG start

**The capsid starts at a codon that is not AUG.** PMID 24503084 Edman-sequenced
the Orsay capsid N-terminus as NKNNT, corresponding to the **AAC codon at RNA2
nt 89–91**, with no AUG 5′ of it in the segment. Verified on HM030971: nt 89–91
are `AAC`. `upstream_ext` is 0, because a terminus with no start codon must not
be scanned for one; `check_met_start.py` reports this feature in 100% of
members and that report is the documented biology. The initiation *mechanism*
is unresolved in the literature and is not claimed.

**CP-delta is a −1 programmed ribosomal frameshift**, so it has no PSSM and is
called from a nucleotide reference set. Slippery heptamer `G_UCA_AAC` at RNA2
nt 1255–1261, immediately 5′ of the capsid stop, stem-loop ~8 nt downstream, no
pseudoknot; the transframe peptide was MS/MS-sequenced and the protein is a
virion head fiber (PMIDs 24503084, 32817218, 28241071). Verified on HM030971:
nt 1255–1261 are `GTCAAAC` and the capsid stop is the `TAG` at 1262–1264.
RefSeq's `join(89..1261,1261..2385)` repeats nt 1261 because that is the
slipped base — a frameshift, not a splice.

> **Two slippery heptamers, not one — measured, and it mattered.**
> Building the reference set on the published Orsay heptamer alone gave 24
> references from Le Blanc and Orsay only, and 23 genomes reported *no
> insertion sits in the slip site* — every Santeuil and Melnik genome in the
> module. Reading the 7 nt immediately 5′ of the capsid stop in each genome
> separates them exactly by species: **Orsay and Le Blanc carry `GTCAAAC` in
> 26 of 26; Santeuil and Melnik carry `GGAAAAC` in 24 of 24.** PMID 24503084
> notes only that Le Blanc shares Orsay's slippery sequence and Santeuil's
> differs, without giving Santeuil's. With both, the set is **47 references
> covering all four species**, every one validated by reproducing its own
> deposited fusion protein byte-for-byte at a confirmed slip site.

**Free delta is a separate initiation and it is real.** Its own AUG at RNA2
nt 1345 (confirmed `ATG` on HM030971), 332–398 aa. PMID 29743360 settled it
with an ATG→CTG mutant that still makes normal CP-delta but no free delta and
has a non-lytic egress defect. It is **absent from RefSeq NC_028098.1** and
present only in HM030971.2, which is why only 3 deposits exist. 46 of 46
extracted members are exact C-terminal substrings of the fusion from their own
genome, beginning at fusion residue 419–434 — Orsay's published 420 among them.

Verified on HM030971: CP **89–1264** and delta **1345–2385**, both exactly as
deposited.

**Two RefSeq annotation traps are defused in the rules.** Santeuil NC_015070.1
labels the *fusion* `protein delta` and gives the *capsid* `/note="protein A"`,
so `protein delta` covers both a 765 aa fusion and a ~346 aa free protein, and
`protein A` appears on an RNA2 capsid. Length routing separates the first;
`protein A` binds nothing in this module and is left to the homology rescue.

---

## 4. What the measurements say

### Held-out panels — complete genomes only

A panel genome here has to clear three conditions, and none implies the others:
it carries both segments; the two contigs are **different** segments; and it is
not a reference. The second is not obvious and cost a correction:
`merge_segments.py` merges any two records sharing a genome_name, so Santeuil
`977913.15` (MN131015 + MN131016, **both 2896 nt and both the polymerase**)
became a 5792 nt "genome" with no capsid, and reported CP, CP-delta and delta
as uncalled. That would have withheld `copy_num` from three features on the
strength of an artefact. A length window cannot catch it — Le Blanc's segment
windows are 2896–3628 and 2525–2937, so two 2896 nt polymerases satisfy both —
so each contig is assigned by tblastn of the module's own collections against
it. The same test caught 13 Betanodavirus pairs and 2 crustacean pairs.

| module | panel | good | poor | features at 100% |
|---|---|---|---|---|
| Alphanodavirus | 4 | 4 | 0 | A, alpha, beta, gamma (B1/B2 at 75%) |
| Betanodavirus | 80 | 80 | 0 | A, B1, B2, CP |
| Macrobrachium_rosenbergii_nodavirus | 4 | 4 | 0 | A, B2, CP |
| Le_Blanc_nodavirus | 19 | 19 | 0 | A, CP, CP-delta, delta |

**107 complete held-out genomes, 107 Good, 0 Poor, 0 duplicate calls.**

The Alphanodavirus and crustacean panels are 4 genomes because only 399 of
1,965 isolate names in this family have both segments deposited. That is a
limitation of the deposits and it is why the cutoffs for those two modules were
left alone (below).

### Nesting is all declared architecture

`qc_duplicate_calls.py` reports **no EXACT duplicates in any module** and no
undeclared nesting. Every containment is real biology: B1 and B2 inside A
(in-frame fragment and +1 overlapping gene), beta and gamma inside alpha
(cleavage products), CP and delta inside CP-delta (frameshift fusion).

`qc_cross_feature.py` reports MISLABEL rows for exactly these pairs — B1 against
A in Betanodavirus, CP and DELTA against CP_DELTA in Le_Blanc. They are 100%
identical over the shorter member's length **by construction**, which is what
the architecture means. The tool has no `--protect`, so the rows stand and are
answered here rather than suppressed. `audit_bins.py`, which does take
`--protect`, reads Unresolved: 0 on all four modules with those pairs declared.

### `copy_num`, set from measurement

Fourteen features are marked essential, and two deliberately are not:

```
Alphanodavirus                       A, ALPHA, BETA, GAMMA
Betanodavirus                        A, B1, B2, CP
Macrobrachium_rosenbergii_nodavirus  A, B2, CP
Le_Blanc_nodavirus                   A, CP, DELTA
```

**Alphanodavirus B1 and B2 are withheld** at 75%. The single non-caller,
`756850.3` (Alphanodavirus HB-2007/CHN), is a *complete* genome, and the rule
is that one complete genome failing to call a feature is enough to withhold it
— otherwise the flag reports a profile gap as a genome defect.

**CP_DELTA never takes `copy_num`**: a special feature is emitted by a
different program.

### Bit cutoffs, raised from measured scores

The generated cutoff is a floor, not a threshold. Each feature's cutoff was set
to **half the lowest score among its genuine panel calls, capped at 500**;
features with fewer than five panel calls were left alone, because one low
score is not a distribution.

| module | feature | lowest panel score | cutoff |
|---|---|---|---|
| Betanodavirus | A | 1854.7 | 500 (capped) |
| | B1 | 206.1 | 103 |
| | B2 | 136.4 | 68 |
| | CP | 606.7 | 303 |
| Le_Blanc_nodavirus | A | 1859.3 | 500 (capped) |
| | CP | 782.3 | 391 |
| | DELTA | 768.8 | 384 |

Alphanodavirus (4 panel calls) and the crustacean module (4) keep their
generated cutoffs, and that is recorded rather than quietly skipped.
**Re-measured after the raise: every feature still calls 100%, every genome
still Good.**

---

## 5. Limits, stated plainly

- **Alphanodavirus B1 is callable only where its start codon can be placed.**
  Four members from three species, after 8 of 13 candidate genomes were
  rejected for having no Met at the aligned start. B1's start is published for
  FHV/BBV (GenBank) and Nodamura (PMID 15163762, protein A residue 908) and for
  nobody else. The profile set is correspondingly thin, and B1 is not marked
  essential.
- **The Pariacoto minor capsid protein is not modelled.** It is a real virion
  protein initiated at the second AUG of the capsid ORF (PMID 13679619) and it
  is alpha minus 24 N-terminal residues — 100% identical to alpha over its whole
  length. No PSSM can separate the two; this is the nested alternative-start
  trap the skill documents for VP1/VP2.
- **Santeuil and Melnik RNA2 deposits are 5′-truncated.** Their contigs begin
  inside the capsid gene, so the proteins start at Lys rather than the verified
  Asn. They are excluded from the CP collection rather than used to train the
  profile, which costs Santeuil representation in that feature.
- **Two of 54 coverage exemplars route and call nothing.** `1168598.3`
  (*Betanodavirus* S1.Tun03.11) and `1110530.25` (*Sparus aurata*
  betanodavirus) hit the references at 951 and 761 bits, far above the 150-bit
  floor, and then return `stop-free 0` for every profile — the annotator finds
  no stop-free ORF to match. Both are partial gene deposits whose reading
  frames are broken relative to the profiles. Not repaired.
- **The BBV RNA2 second ORF ("URF", 72 aa) is not modelled.** PMID 6548308
  predicted it in 1984 and no product has been reported since; one deposit
  exists family-wide.
- **Three modules' rarefaction curves rest on small pools.** Only Betanodavirus
  reaches the 40-protein floor from its held-out panel; the others report the
  raw collapse from the whole-taxon annotation instead, and the pages say so.

## 6. Citations

Every PMID in the four `features.json` files sits in
**`PMID_claude_generated`**, never `PMID`. All were resolved through NCBI
eutils and their titles read against the claim; `check_literature_gaps.py
--verify-pmids` prints them. Expect this module to ship with `PMID` empty on
every feature until a human curator moves them.
