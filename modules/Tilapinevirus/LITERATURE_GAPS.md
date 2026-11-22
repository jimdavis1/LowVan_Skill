# Literature proteins vs BV-BRC calls — Tilapinevirus

Collections are built from deposited strings, so a protein that is
real biology but unannotated in these genomes is invisible to the
build. This records what was probed for and what was found.

| protein | deposits | verdict | binds to | source | note |
|---|---|---|---|---|---|
| PB1 polymerase subunit | 23 | PRESENT | S1_PB1 | PMID 27048802 | Segment 1, 519 aa, catalytic SDD at 288-290. The only TiLV protein with recognisable homology to anything outside the family. |
| PB2 polymerase subunit | 0 | ABSENT | — | PMID 38066000 | Segment 2, 457 aa. Assigned by fold in the cryo-EM polymerase structure. Deposited only as 'hypothetical protein'. |
| PA polymerase subunit | 0 | ABSENT | — | PMID 38066000 | Segment 3, 419 aa. Assigned by fold. Deposited only as 'hypothetical protein'. |
| nucleoprotein | 0 | ABSENT | — | PMID 35107373 | Segment 4, 356 aa. Identified biochemically then structurally bound to RNA. Deposited only as 'hypothetical protein'. |
| S9-F3 nuclear export protein | 0 | ABSENT | — | PMID 41513822 | THE ELEVENTH PROTEIN. 105 aa from frame 3 of segment 9, overlapping the primary S9 ORF, confirmed by mass spectrometry (4 unique peptides, 49% coverage) and carrying a CRM1-dependent nuclear export signal. Expected ABSENT: BV-BRC holds no deposit of it. Sourced by extraction from the genomes instead. |
| spliced mRNA product | 0 | ABSENT | — | PMID 41513822 | DELIBERATELY ABSENT. Orthomyxoviruses make nuclear export proteins from spliced segment 7/8 mRNAs, but PMID 41513822 looked for splicing in TiLV and reports none: S9-F3 arises by leaky ribosomal scanning of an alternative frame. ABSENT here is the correct result and confirms no influenza splicing rule has leaked in. |
| PA-X frameshift product | 0 | ABSENT | — | PMID 41513822 | DELIBERATELY ABSENT. PMID 41513822 states it could not identify an obvious PA-X-like frameshift signal in segment 3, so no frameshift product is declared. |
| matrix protein | 0 | ABSENT | — | PMID 41513822 | DELIBERATELY ABSENT. No amnoonvirid matrix protein has been described; segments 5-10 have no assigned function beyond predicted motifs. |

## The ABSENT verdicts, and what each one means

This module is the case the gap check was built to surface: **the source
vocabulary and the module's contents disagree, and the vocabulary is wrong.**

| probe | verdict | decision |
|---|---|---|
| PB1 | PRESENT | 23 deposits name it; the only TiLV protein the source names |
| PB2 | ABSENT | **but the module calls it** -- 29 collection members, every one deposited as `hypothetical protein`. Absent from the VOCABULARY, not from the module. Recovered by homology to the Til-4-2011 reference segment set. |
| PA | ABSENT | same -- 37 members, all `hypothetical protein` |
| nucleoprotein | ABSENT | same -- 25 members, all `hypothetical protein` |
| S9-F3 | ABSENT | genuinely absent from BV-BRC: **zero deposits**. Sourced by extraction from frame 3 of each genome's own segment 9. 32 extractions, all exactly 105 aa. |
| spliced mRNA product | ABSENT | **expected and correct.** PMID 41513822 searched for splicing and found none; S9-F3 arises by leaky ribosomal scanning. ABSENT here confirms no orthomyxovirus splicing rule has leaked into this module. |
| PA-X frameshift product | ABSENT | **expected and correct.** PMID 41513822 could not identify a PA-X-like frameshift signal in segment 3. |
| matrix protein | ABSENT | **expected and correct.** No amnoonvirid matrix protein has been described. |

## What was probed, and the edge of this check

Probed: the four proteins with assigned functions (PB1, PB2, PA, NP), the
eleventh protein S9-F3, and the three mechanisms an orthomyxo-like virus might
plausibly use that TiLV does **not** (splicing, PA-X frameshifting, a matrix
protein).  Segments 5-10 were not probed for function because none has one in
the literature; they are declared as uncharacterized lineage-specific proteins
under the literature's own S5..S10 designations.

The check cannot find a protein that is both unannotated in BV-BRC and absent
from this list.  The known residual risk is a twelfth protein in an
alternative frame of a segment other than 9: PMID 41513822 searched a six-frame
translation of the whole genome at >=9 codons and reported eleven proteins, so
that risk is small, but it is not zero.
