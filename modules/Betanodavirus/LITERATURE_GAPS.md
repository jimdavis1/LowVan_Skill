# Literature proteins vs BV-BRC calls — Betanodavirus

Collections are built from deposited strings, so a protein that is
real biology but unannotated in these genomes is invisible to the
build. This records what was probed for and what was found.

| protein | deposits | verdict | binds to | source | note |
|---|---|---|---|---|---|
| protein A (RdRp) | 526 | PRESENT | A, UNASSIGNED | PMID 10580064 | RNA1 replicase, 981-983 aa. Carries mitochondrial targeting signals (GGNNV aa 215-255, PMID 15163743; AHNV aa 1-40 and 225-246, PMID 17602779). No cleavage products described. |
| B1 | 2 | BELOW_FLOOR | B1 | PMID 19136133 | 111 aa (336 nt gene) in RGNNV, in frame with protein A and sharing its stop, from sgRNA3. An anti-necrotic death factor. NO START CODON IS PUBLISHED and no RefSeq annotates a B1 CDS, so this collection is extracted from genomes under a structural gate. |
| B2 | 46 | PRESENT | B2 | PMID 16186236 | 75 aa, +1 frame relative to protein A, from sgRNA3, which begins at RNA1 nt 2730 by premature termination rather than splicing. RNAi antagonist and necrotic death factor. |
| coat protein | 1096 | PRESENT | CP, UNASSIGNED | PMID 9672600 | 338-340 aa, a single UNCLEAVED chain. The RNA2 ORF carries the catalytic D-75 residue but lacks the autocatalytic cleavage site, so there is no alpha/beta/gamma maturation in this genus. |
| capsid beta/gamma maturation products | 0 | ABSENT | — | PMID 26491970 | DELIBERATELY ABSENT. Betanodavirus CP is not proteolytically cleaved (PMID 9672600), and the T=3 crystal structure is built from the intact 338-residue chain (PMID 26491970). ABSENT here is the correct and expected result: it confirms the alphanodavirus cleavage rule has not leaked into this module. |
