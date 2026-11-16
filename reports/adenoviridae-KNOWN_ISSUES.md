# Adenoviridae modules -- known issues (logged 15 Nov 2026, shipped as is per curator)

## 1. Mature-peptide cleaved ends are bound to the tBLASTn HSP
A cleaved terminus has no start/stop to scan to, so it lands wherever the profile
alignment stops. The annotator has no step that moves it onto the AVP motif
([MILV]xGG^x / [MILV]xGx^G). Measured with audit/matpep_ends.py (product must be an
exact substring of its own genome's precursor and each cut end must sit on a motif).

| Module | Product | Ends on the site |
|---|---|---|
| Siadenovirus | VI, VII | 100%, 100% |
| Mastadenovirus_MultiFiber | VI, VII, VIII_N, VIII_C, TP, IIIA | 95-100% |
| Mastadenovirus (full panel) | IIIA | 86% (3225/3750) |
| Mastadenovirus (subset) | VI 90, VII 98, VIII_N 98, VIII_C 79, TP 93 | |
| Aviadenovirus | VI 96, VIII_N 91, **IIIA 57 (72 full panel)** | |
| Barthadenovirus | VI 82, **IIIA 51 (59 full)**, **VIII_C 10 (7/73)** | |

Barth VIII_C typically starts 5-8 aa past the real site (motif at ~144, product at 149
in BAdV-4-like genomes). Off-site IIIa C-termini are lineages where the profile stops
short of the site (BAdV-3, sea lion, CAdV-2, MAdV-3) or the lineage lacks a C-terminal motif.
Fix if revisited: annotator snaps a cleaved end to the nearest motif within +-N aa
(code change -> annotator_regression.py), or per-lineage profiles trimmed at the site.

## 2. HAdV-D IIIa start disagreement (316 genomes, full panel)
PIIIA CDS starts ~29 aa upstream (MSQQQA...) of the IIIA_MAT start (MQRIMA..., the Met
homologous to HAdV-5 MQRAAM). Precursor and product disagree on the protein's start.

## 3. VII_MAT missing in divergent one-fiber lineages
After the routing fix (rep swap, ROUTEFIX_15nov.md), 8 lemur/bat/tree-shrew-type genomes
lost VII: Masta VII_MAT.1 scores 73.6 bits but no stop-free deliverable HSP. Lineages not
represented in Masta VII profiles.

## 4. Routing: 3 one-fiber Masta genomes still route to MultiFiber
47680.6 (tree shrew), 3075435.12 (SAdV CP001), 2050578.14 (Mastadenovirus sp.); the
25-rep budget is used. 7 HAdV-F genomes with broken fibers are assigned Masta by the
fiber census but route (correctly) to MultiFiber.

## 5. Unspliced calls of spliced genes start off-Met
pol, pTP, IVa2 (and others with upstream exons) are called on the body exon; non-M start
is expected there. UXP 3-exon product: see UXP_SET_ASIDE.md.

## Also set aside / not modelled
Mature products not installed (<60% of precursors carry the site): see RESUME.md
(Avi VII/VIII_C/TP; Barth VII/VIII_N/TP; Sia VIII_N/VIII_C/TP/IIIA). Not modelled anywhere:
pVIc, mu, pVIIn, propeptides, iTP, 52K products.

## 6. Essential genes (copy_num), set 16 Nov 2026 by the curator's rules
apply_copy_num.py's genome-count list was NOT used as is (it made every genome Poor:
uncharacterized features share one annotation string, which the quality tool keys on,
and primate-only E3/E4 genes cleared 85% only because Masta is ~85% human genomes by deposited name).
Rules: never on uncharacterized, special, or lineage-specific features; excluded where a
cutoff, copy-number or coverage struggle is on record. Essential now:
Sia 17, Barth 15, Avi 17, MultiFiber 13, Masta 11 (lists in each features.json).
Subset rescore (916 genomes): good Sia 85.5%, Barth 97.3%, Avi 89.1%, MF 48.3%, Masta 74.4%.
Poor genomes are mostly explained by the sequence: C-less (C->T converted) deposits (Sia
10, Masta 14, MF 2; e.g. HAdV-E 130308.18-21, 0.0% C), N gaps, partial genomes, routing
to the other Masta module, and genes interrupted by stop codons (checked on HAdV-F
130309.457: DBP 775 bits and E1B 55K 954 bits, both with internal stops).

## 8. Features per genome vs. deposits (16 Nov 2026): Barth and Avi under-call
Near-complete genomes, medians per genome (deposited proteins = BV-BRC; module = CDS called, 916 subset):

| Module | Genome | Deposited | Module calls | Difference |
|---|---|---|---|---|
| Siadenovirus | 26.0 kb | 23 | 22 | -1 (natural: smallest genomes) |
| Barthadenovirus | 33.0 kb | 29 | 24 | -5 |
| Aviadenovirus | 43.8 kb | 37 | 29 | -8 |
| Mastadenovirus_MultiFiber | 34.1 kb | 29.5 | 32 | +2.5 |
| Mastadenovirus | 35.2 kb | 32.5 | 34 | +1.5 |

Barthadenovirus misses (518 deposited features unassigned): 22K (item 7), plus recurring
lineage names never declared -- "hypothetical protein" 148, Pac22k 12, MucoreX 12, Slatt 12,
HistinVII 12, "core protein 2" 12. Aviadenovirus misses (4,621 unassigned): ORF0 280,
ORF1B 261, ORF20A 257, ORF19 157, ORF1C 142, ORF29 127, "hypothetical protein" 772.
Fix: declare them as uncharacterized lineage-specific orthogroups (literature check first),
plus 22K in Barth; rebuild both modules. Not done -- awaiting curator decision.

## 9. Spliced features need a >=95% nucleotide reference (curator: leave as is, 16 Nov 2026)
33K, E1A 12S/13S and E4ORF6/7 are called by get_splice_variant_features.pl, which calls a
gene only when the genome matches a curated reference transcript at >=95% nt identity over
>=95% of its length. Lineages with no near-identical reference get no call although the gene
is there: HAdV-C AdV-lumc002 (129951.1438) is 80.7% identical over the 895-nt 33K region to
the nearest HAdV-2 reference and has no 33K call; 30 of 48 species C genomes in the subset lack
33K. Module-wide, 33K is called in 92.6% of near-complete Mastadenovirus genomes.
Literature search for more 33K anchors (16 Nov 2026; k33lit/FINDINGS.md) found none: no paper
maps exact 33K junctions for an uncovered lineage. The only near miss is FAdV-4, PMID 36851744,
which gives figure-only exons; a reference there would rescue 130 Avi genomes. The limit is the
≥95% rule. 721 uncalled near-complete genomes carry the locus (MF 33.2% called, Avi 64.8%),
mostly in divergent types of covered species (F41 vs F40 ~360). Closed per curator: literature
only, otherwise move on. Also open: the shipped K33.fasta (11 Nov) predates the literature-
junction rule and was carried from the old module.

## 7. 22K rebuilt in all five modules (RESOLVED 16 Nov 2026)
Cause: the old 33K collection was still a seed in homology assignment, although 33K is a splice
special with no PSSM. 22K deposits share 33K exon 1, so they came out AMBIGUOUS or were binned to
33K. Barthadenovirus was left with 2 K22 sequences and the feature was dropped (UNBUILDABLE.tsv),
and the K22 collections that shipped held spliced 33K proteins (Masta 18, Sia 2, MF 2).
Fix (k22fix/): each candidate was tblastn'd onto its own genome (one HSP over >=90% = unspliced);
22K = unspliced, starting inside that genome's 100K call in another frame. Collections: Sia 12,
Barth 8, Avi 90, MF 36, Masta 331. PSSMs: 3 / 2 / 8 / 5 / 20. Barth K22 is the short unspliced L4
product (59-66 aa, 33K exon 1 read past the donor), the same structure as HAdV 22K (PMID 22811519).
Masta K22 bit_cutoff 50: profile K22.7 hit protein VI at 30.46 in 1123958.3-.5 and made a second,
same-frame K22 call; the lowest real 22K scored 72.87. Other modules stay at 30 (Barth 25).
Test (plain annotator, 218 genomes, k22fix/test/): one K22 call at most per genome in every module.
Avi starts moved 200-400 nt downstream in 22/41, now on the deposited 22K starts (old profile: 7
members, some over-long). Masta 1907210.3299/.3356 lose K22: the 22K frame has an internal stop
(coverage before it 0.38/0.43); the old calls were a partial and a duplicate of the 100K tail.
Remaining (pre-existing, not from this fix): K22 and K100 share a frame and stop in Avi 172859.66,
130663.266/.267 and Masta 1907210.3357/.3305. The real 22K aligns there, so the K100 call has read
through into the 22K frame (a frameshift in the assembly is likely). In 130663.266 the K22 call
also extends upstream to 375 aa, past max_len 244: max_len is not applied after extension.
