# Literature proteins vs BV-BRC calls — Mastadenovirus

Collections are built from deposited strings, so a protein that is
real biology but unannotated in these genomes is invisible to the
build. This records what was probed for and what was found.

| protein | deposits | verdict | binds to | source | note |
|---|---|---|---|---|---|
| E1A | 2067 | PRESENT | E1A, DBP, E1A_243R, UNC6, UNC4, UNASSIGNED, PIX, UNC51 | PMID_claude_generated 1833268 | The first gene expressed; 13S and 12S mRNAs give the 289R and 243R products by alternative 5' splice site within one transcription unit. |
| E1A 243R | 157 | PRESENT | E1A_243R | PMID_claude_generated 1833268 | The 12S product. 243R is 289R minus an internal 46-residue block (conserved region 3), so the two are ~100% identical outside it. Declared as a SPLICE feature, not a second PSSM, because no profile can separate them. |
| E1B 19K | 1494 | PRESENT | E1B19K | PMID_claude_generated 15302210 | Bcl-2 homologue; blocks E1A-induced apoptosis. |
| E1B 55K | 1059 | PRESENT | E1B55K, UNASSIGNED | PMID_claude_generated 17079297, 20861261 | p53-SUMO1 E3 ligase; with E4orf6 forms the cullin-based ubiquitin ligase that degrades p53. |
| protein IX | 2071 | PRESENT | PIX, E1B55K | PMID_claude_generated 2583109, 33789897 | Minor capsid protein, the capsid cement sitting in the hexon facets. Present in Mastadenovirus, Aviadenovirus and Atadenovirus; ABSENT in Siadenovirus. |
| core protein V | 2231 | PRESENT | PV, PVII | PMID_claude_generated 33252742 | Mastadenovirus-specific core protein; bridges the DNA-pVII core to the capsid. Not present in the other genera. |
| E3 12.5K | 1147 | PRESENT | E3_12_5K, UNASSIGNED | PMID_claude_generated 11022786 |  |
| E3 gp19K | 1090 | PRESENT | E3_GP19K, E3_CR1A, UNASSIGNED | PMID_claude_generated 8249282, 8350398 | ER-resident type I membrane glycoprotein that binds MHC class I and retains it in the ER. CARRIES A CLEAVED N-TERMINAL SIGNAL PEPTIDE -- the only signal peptide in the mastadenovirus proteome -- plus a C-terminal ER-retention signal. |
| E3 CR1-alpha | 1442 | PRESENT | E3_CR1A, E3_RIDA, UNASSIGNED | PMID_claude_generated 8627757, 11022786 |  |
| E3 CR1-beta | 1697 | PRESENT | E3_CR1B, UNASSIGNED, DISCARDED | PMID_claude_generated 11022786, 11752165 |  |
| E3 CR1-gamma | 1304 | PRESENT | E3_CR1G, E3_CR1B | PMID_claude_generated 11022786 |  |
| E3 CR1-delta | 440 | PRESENT | E3_CR1D, E3_CR1G | PMID_claude_generated 11022786 |  |
| E3 6.7K | 19 | PRESENT | L1_13_6K, E3_GP19K | PMID_claude_generated 15507617 | Required WITH the RID complex for internalisation and degradation of TRAIL receptor 2. Found by this probe pass, not by the annotation strings. |
| E3 RID-alpha | 1898 | PRESENT | E3_RIDA, K100, E3_RIDB, UNC35, DISCARDED | PMID_claude_generated 14557654, 16227281 | E3 10.4K. With RID-beta removes Fas, EGFR and TNFR1 from the cell surface. |
| E3 RID-beta | 2091 | PRESENT | E3_RIDB, E3_14_7K, PIX, E4ORF2, UNC21, K22, PVIII, UNASSIGNED, DISCARDED | PMID_claude_generated 16227281 | E3 14.5K, the other half of the RID complex. |
| E3 14.7K | 1978 | PRESENT | E3_14_7K, E3_RIDB, E1A | PMID_claude_generated 17024246 |  |
| adenovirus death protein | 17 | PRESENT | E3_CR1B | PMID_claude_generated 15302210 | E3 11.6K. Required for efficient lysis and release at the end of infection. |
| E4orf1 | 968 | PRESENT | E4ORF1, E4ORF6, UNASSIGNED | PMID_claude_generated 38451084 | PI3K-Akt activator; targeted by p53. |
| E4orf2 | 862 | PRESENT | E4ORF2 | PMID_claude_generated 11134268 |  |
| E4orf3 | 728 | PRESENT | E4ORF3 | PMID_claude_generated 16103189, 16254336 | Disrupts PML nuclear bodies and drives the MRN complex into cytoplasmic aggresomes. |
| E4orf4 | 746 | PRESENT | E4ORF4 | PMID_claude_generated 31792953 |  |
| E4orf5 | 20 | PRESENT | E4ORF5, E4ORF6 | PMID_claude_generated 11134268 | Found as accessory cluster UNC1 in triage pass 1 (59 sequences, median 151 aa, labelled 'Orf5' and 'E4 17 kDa protein'); promoted to a declared feature in pass 2. |
| E4orf6 | 1421 | PRESENT | E4ORF67, E4ORF6, DBP, E4ORF4 | PMID_claude_generated 17079297, 20566845 | E4 34K. With E1B 55K forms the p53 ubiquitin ligase and stimulates late mRNA nuclear export. |
| E4orf6/7 | 1663 | PRESENT | E4ORF67 | PMID_claude_generated 1836562, 8035503 | A SPLICED product: exon 1 of orf6 joined to orf7. Two exons, one intron. |
| long fiber | 310 | PRESENT | FIBER1, FIBER | PMID_claude_generated 7975882 | The enteric species HAdV-F40 and F41 carry TWO fiber genes, a long and a short. This is the claim that put FIBER1/FIBER2 in the Mastadenovirus module; triage pass 1 had left them in the accessory pool as UNC2 (median 555 aa) and UNC5 (median 368 aa). |
| short fiber | 271 | PRESENT | FIBER2, FIBER1, FIBER, UNC30 | PMID_claude_generated 7975882 |  |
| U exon protein | 1007 | PRESENT | UXP, DISCARDED | PMID_claude_generated 17881437, 20739525 | A late l-strand product overlapping E3 antisense; unrecognised until 2007, which is exactly the shape of protein this probe pass exists to catch. |
| i-leader protein | 235 | PRESENT | I_LEADER, L1_13_6K | PMID_claude_generated 16439556 |  |
| L1 13.6K | 474 | PRESENT | L1_13_6K, E4ORF3, UNASSIGNED, E4ORF4 | PMID_claude_generated 3238923 |  |
| VA RNA I/II | 0 | ABSENT | — | PMID_claude_generated 37434327 | NOT MODELLED AND NOT A GAP: the VA RNAs are non-coding Pol III transcripts. They have no protein product and cannot carry a PSSM. Recorded here so the absence is a decision rather than an oversight. |
