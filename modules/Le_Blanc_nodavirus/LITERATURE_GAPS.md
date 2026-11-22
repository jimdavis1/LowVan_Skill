# Literature proteins vs BV-BRC calls — Le_Blanc_nodavirus

Collections are built from deposited strings, so a protein that is
real biology but unannotated in these genomes is invisible to the
build. This records what was probed for and what was found.

| protein | deposits | verdict | binds to | source | note |
|---|---|---|---|---|---|
| protein A (RdRp) | 50 | PRESENT | A | PMID 24503084 | The ONLY CDS on RNA1 in this group. Orsay 937, Le Blanc 942, Santeuil 982 aa. |
| capsid protein | 70 | PRESENT | CP, CP_DELTA | PMID 24503084 | 388-391 aa, initiated at a NON-AUG codon: the Edman-sequenced mature N-terminus NKNNT corresponds to the AAC codon at RNA2 nt 89-91, with no AUG 5' of it in the segment. |
| CP-delta fusion | 51 | PRESENT | CP_DELTA, DELTA | PMID 24503084 | Made by a -1 programmed ribosomal frameshift at the slippery heptamer G_UCA_AAC, RNA2 nt 1255-1261, with a stem-loop about 8 nt downstream and no pseudoknot. The transframe peptide was sequenced by MS/MS and the protein is a virion head fiber. |
| free delta | 51 | PRESENT | CP_DELTA, DELTA | PMID 29743360 | Initiated at its own AUG (Orsay RNA2 nt 1345) and equal to the C-terminal part of the fusion. Existence as a free protein settled by an ATG-to-CTG mutant that still makes CP-delta but has a non-lytic egress defect. Absent from RefSeq NC_028098.1 and present only in HM030971.2. |
| B2 | 0 | ABSENT | — | PMID 21283608 | DELIBERATELY ABSENT. PMID 21283608 searched RNA1's 3' end for a +1 ORF with B2's properties and found none; PMID 23468484 detected no silencing-suppressor activity for either Orsay protein proposed as one. ABSENT here is the expected and required result and is what separates this module from every other nodavirus module. |
| B1 | 0 | ABSENT | — | PMID 21283608 | DELIBERATELY ABSENT, for the same reason as B2: RNA1 carries a single CDS and no subgenomic RNA3 has been described in this group. |
