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
    P, X     upstream_ext 1    their own ORFs

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
