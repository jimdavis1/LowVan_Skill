# Aviadenovirus -- build notes (Adenoviridae, 16 Nov 2026)

Aviadenovirus. Fiber orthogroups FIBER_A1 to FIBER_A33, one per fiber class; the most common patterns are A1 alone (183 genomes) and A2-A1 (128). LIPASE/LIPASE2 and ORF14/14A/14B are split by position.

## Adenoviridae: what all five modules share

**Partition, by sequence only (10 Nov 2026).** The genus is the best hexon match
from another genome. Fiber classes are 40% mmseqs classes of intact fiber ORFs,
named by literature orthology (Podgorski 2023, PMID 36610683). A Mastadenovirus
genome carrying any multi-fiber class goes to Mastadenovirus_MultiFiber. Never
trust genome names: 80 genus-blank genomes were placed by hexon, and 3
"unidentified adenovirus" deposited as Mastadenovirus are penguin siadenovirus.
Ichtadenovirus (2 genomes) is out of scope.

| Module | Genomes | Features | Rep contigs | Full-panel routing to own module |
|---|---:|---:|---:|---|
| Mastadenovirus (one fiber) | 3,821 | 54 | 25 | 3,776/3,821 before the rep swap; 36 of the 45 misroutes fixed after it |
| Mastadenovirus_MultiFiber | 591 | 47 | 25 | 100% |
| Aviadenovirus | 452 | 55 | 25 | 100% |
| Barthadenovirus | 74 | 48 | 22 | 100% |
| Siadenovirus | 69 | 25 | 12 | 100% |

**Validation.** Full panel of all 5,007 genomes (14-15 Nov 2026), then the
scan-to-stop and cutoff fixes verified on a 916-genome audit subset. Every
non-spliced CDS call now ends on a stop codon, except calls that stop at an N
gap or an IUPAC base in the assembly (0-73 per module). No CDS feature is called
twice at one locus except the known multi-copy cases below.

**Curator decisions (Jim Davis).**
- CDS features scan to the stop codon (`downstream_ext 1`). Mature peptides set
  both flags from their own termini.
- Fibers: annotation "Fiber protein", symbol FIBER. Lineages are separate
  orthogroup features (Avi FIBER_A*, MultiFiber FIBER_ANC/CL*), not renamed (16 Nov).
- Paralogue families are split by position into orthogroup features, never
  called by location. CR1 layout: one feature per E3 slot (CR1A, CR1B), plus
  CR1G_B, CR1B_E2 and E3_ADPL (signed off 16 Nov).
- No literature function -> "Uncharacterized lineage-specific protein", symbol
  unnumbered (RH, CR1-alpha, ...). Gene symbol E3-10.9K for CR1-delta.
- Splice variants: E1A 12S/13S, E4ORF6/7, 33K (K33); Siadenovirus DBP.
  Junctions come from the literature, extended to other lineages by aligning
  the splice donor and acceptor sites, and every product is translated and
  checked in frame.
- The UXP three-exon product is set aside: get_splice_variant_features.pl
  supports one junction only. UXP is called as the U-exon ORF.

**Mature peptides.** Cut at the HAdV-2/5 protease (AVP) sites from the literature
(PMIDs 28978703, 25421887, 6336325, 34774568, 2691633, 2800334, 24227847).
A site is installed only where at least 60% of the module's precursors carry it.

| Module | Products |
|---|---|
| Mastadenovirus, MultiFiber | VI, VII, VIII-N, VIII-C, TP, IIIa |
| Aviadenovirus | VI, VIII-N, IIIa |
| Barthadenovirus | VI, VIII-C, IIIa |
| Siadenovirus | VI, VII |

## Known issues (logged, shipped as is)

1. **Cleaved ends stop where the tBLASTn HSP stops.** The annotator does not
   move them onto the protease motif. Calls that land on the site:

   | Module | Product | On the site |
   |---|---|---|
   | Siadenovirus | VI, VII | 100% |
   | MultiFiber | all products | 95-100% |
   | Mastadenovirus | IIIa | 86% |
   | Aviadenovirus | IIIa | 72% |
   | Barthadenovirus | IIIa | 59% |
   | Barthadenovirus | VIII-C | 10% |

   Barthadenovirus VIII-C starts 5-8 aa past a site that is present.
2. **HAdV-D IIIa start.** In 316 genomes the PIIIA precursor starts about 29 aa
   upstream of the IIIa product, which starts at the Met homologous to the
   HAdV-5 MQR start.
3. **Masta VII missing in 8 divergent one-fiber genomes** (lemur, bat,
   tree-shrew types). Those lineages are not represented in the VII profiles.
4. **Routing leftovers.** 3 one-fiber genomes still route to MultiFiber:
   47680.6, 3075435.12 and 2050578.14. 7 HAdV-F genomes with broken fibers are
   assigned Mastadenovirus by the fiber census but route, correctly, to
   MultiFiber.
5. **Spliced genes called on the body exon start off-Met.** This is expected
   for pol, pTP and IVa2.
6. **Multi-copy calls that are real.** Barthadenovirus FIBER and LH3; FrAdV-1
   ORF7; Aviadenovirus ORF17 in one genome.
7. **Not modelled.**
   - The mature products whose site is in fewer than 60% of precursors:
     Avi VII/VIII-C/TP; Barth VII/VIII-N/TP; Sia VIII-N/VIII-C/TP/IIIa.
   - Anywhere: pVIc, mu, pVIIn, the propeptides, iTP, and the 52K products.

## Essential genes and pages (added 16 Nov 2026)

**Essential (copy_num 1), set by the curator's rules:** DBP, HEXON, IVA2, K100, K22, K52, ORF22, PENTON, PIIIA, POL, PROTEASE, PTP, PVI, PVIII, UXP, VIII_N, VI_MAT.
Never essential: uncharacterized proteins, special (spliced) features, lineage-specific
genes, and any feature with a cutoff, copy-number or coverage struggle on record.
apply_copy_num.py's genome-count list was not used as is: the quality tool keys essential
features by annotation string, so features sharing "Uncharacterized lineage-specific
protein" collide, and genome-count rates made primate-only genes essential in Mastadenovirus.

**Quality on the 916-genome audit subset:** 115 good / 14 poor
(89.1% good). The subset over-samples unusual genomes. Poor genomes are
mostly explained by the sequence: C-less deposits, N gaps, partial genomes and genes
interrupted by stop codons.

**Pages** (in reports/): aviadenovirus-coverage-audit.html, aviadenovirus-pssm-registry.html,
aviadenovirus-string-collapse.html, aviadenovirus-vocabulary-saturation.html; facts in
reports/generators/aviadenovirus_audit.json.

**Under-calling vs. deposits (16 Nov 2026).** Near-complete genomes: 37 deposited proteins
per genome (median) against 29 module calls. Recurring ORFs never declared: ORF0 280,
ORF1B 261, ORF20A 257, ORF19 157, ORF1C 142, ORF29 127 deposited features (and 772
"hypothetical protein"). Candidates for uncharacterized lineage-specific orthogroups.

**22K rebuilt (16 Nov 2026).** The old 33K collection was still a homology seed although 33K is a splice special with no PSSM, so 22K deposits (which share 33K exon 1) went AMBIGUOUS or to 33K, and the K22 collections picked up spliced 33K proteins. Rebuilt structurally: each candidate tblastn'd onto its own genome (one HSP over >=90% = unspliced), and 22K = unspliced, starting inside that genome's 100K call in another frame. Tested with the plain annotator on a lineage-diverse set: no genome has more than one K22 call. K22: 90 sequences, 8 PSSMs, bit_cutoff 30 (was 7 sequences). 46 test genomes: 22 starts move 200-400 nt downstream onto the deposited 22K starts. K22 and K100 share a frame and stop in 172859.66 and 130663.266/.267 (pre-existing): the 100K call reads through into the 22K frame.
