# Literature proteins vs BV-BRC calls — Barthadenovirus

Collections are built from deposited strings, so a protein that is
real biology but unannotated in these genomes is invisible to the
build. This records what was probed for and what was found.

| protein | deposits | verdict | binds to | source | note |
|---|---|---|---|---|---|
| p32K | 27 | PRESENT | P32K | PMID_claude_generated 28943338, 33789897 | Genus-specific MINOR COAT protein; alpha-helical clusters beneath the capsid facet. OAdV-7 NC_004037 complement(238..1098), 286 aa. |
| LH1 | 29 | PRESENT | LH1, LH3 | PMID_claude_generated 33789897 | NC_004037 1143..1505, 120 aa. No function established -- keeping the number. |
| LH2 | 23 | PRESENT | LH2 | PMID_claude_generated 28943338, 33789897 | STRUCTURAL component, previously unrecognised as one: an extended polypeptide driving a large conformational change in vertex protein IIIa relative to mastadenoviruses. NC_004037 1490..1873, 127 aa. |
| LH3 (deposited as 'E1B 55K') | 22 | PRESENT | LH3, LH2, UNASSIGNED | PMID_claude_generated 28943338, 33789897 | CAPSID CEMENTING PROTEIN. Trimeric beta-helix, bacteriophage-tailspike fold; binds the capsid via the same triskelion motif mastadenovirus pIX uses, and mastadenovirus E1B-55K was EXAPTED from it. BOTH reference genomes annotate this locus as 'E1B 55K' (NC_004037 1916..3064, 382 aa; NC_009989 2169..3290, 373 aa) and atadenoviruses have NO E1B region at all. 13 such deposits hit LH3 at 42-77% identity full length. |
| RH0 | 3 | PRESENT | RH0 | PMID_claude_generated 12504581 | MISSING from the feature set until now. Not in either reference genome; 3 deposits in the dump, one as 'F-box related protein RH0' which fits the RH family. No paper covers it -- keeping the number. |
| RH1 | 11 | PRESENT | RH1, RH4 | PMID_claude_generated 12504581 | F-BOX PROTEIN. MISSING from the feature set until now -- 7 deposits went to UNASSIGNED. NC_004037 complement(29130..29459), 109 aa. |
| RH2 | 22 | PRESENT | RH2 | PMID_claude_generated 12504581 | F-BOX PROTEIN. NC_004037 complement(28486..29016), 176 aa. 7 of its 22 deposits had been consolidated away into RH4. |
| RH3 | 4 | PRESENT | RH3 | PMID_claude_generated 12504581 | Not in OAdV-7 and not an F-box member. 4 deposits (160 aa) in the dump, ALL of which consolidation had moved into RH4, leaving RH3 empty. Not an absent gene -- a paralogue wrongly merged. |
| RH4 | 24 | PRESENT | RH4 | PMID_claude_generated 12504581 | F-BOX PROTEIN. NC_004037 complement(27845..>28318), 157 aa. |
| RH5 | 21 | PRESENT | RH5 | PMID_claude_generated 12504581 | NC_004037 complement(27244..27840), 198 aa. Carries NO F-box note in the reference annotation, unlike RH1/2/4/6. No function established -- keeping the number. |
| RH6 | 7 | PRESENT | RH6 | PMID_claude_generated 12504581 | F-BOX PROTEIN. MISSING from the feature set until now. NC_004037 complement(26688..27161), 157 aa. |
| E4.1 | 34 | PRESENT | E4_1 | PMID_claude_generated 12504581 | NC_004037 complement(25225..25653), 142 aa. No function established. |
| E4.2 | 31 | PRESENT | E4_2 | PMID_claude_generated 12504581 | Shares a conserved motif with E4.3 and with mastadenovirus E4 34K, per the reference annotation. NC_004037 complement(24563..25225), 220 aa. |
| E4.3 | 32 | PRESENT | E4_3 | PMID_claude_generated 12504581 | Shares a conserved motif with E4.2 and mastadenovirus E4 34K. NC_004037 complement(23853..24566), 237 aa. |
| U exon protein | 37 | PRESENT | UXP | PMID_claude_generated 17881437 | NC_004037 complement(<22106..22266), 54 aa; 'continues upstream, probably spliced'. |
| protein IX | 0 | ABSENT | — | PMID_claude_generated 28943338, 33789897 | ABSENT, and this CORRECTS an earlier note in this file which claimed pIX was present because 'the atadenovirus capsid structure resolved it together with p32K'. It did not. The structure resolves the GENUS-SPECIFIC proteins LH2, p32K and LH3, and its point is that LH3 occupies pIX's structural role. Neither reference genome annotates a pIX. |
