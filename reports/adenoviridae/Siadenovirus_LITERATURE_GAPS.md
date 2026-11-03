# Literature proteins vs BV-BRC calls — Siadenovirus

Collections are built from deposited strings, so a protein that is
real biology but unannotated in these genomes is invisible to the
build. This records what was probed for and what was found.

| protein | deposits | verdict | binds to | source | note |
|---|---|---|---|---|---|
| sialidase (ORF1) | 48 | PRESENT | SIALIDASE | PMID_claude_generated 36016389 | Genus-defining gene -- the 'sia' in Siadenovirus. High similarity to BACTERIAL sialidases; confirmed structural component of the TAdV-3 virion. Whether the domain is catalytically active is explicitly open. TAdV-3 NC_001958 399..1952, 517 aa. |
| hyd (hydrophobic protein) | 20 | PRESENT | HYD | PMID_claude_generated 41960683 | Genus-specific. Function UNKNOWN -- named for a property. mRNA uses the tripartite leader and shares a polyA signal with ORF1. |
| E3 | 46 | PRESENT | E3 | PMID_claude_generated 36016389 | A single E3 of its own, not the mastadenovirus cassette. NC_001958 21213..22115, 300 aa, but the reference's OWN note says 'identification based on very limited similarity'. No functional study found. Not asserting a function. |
| ORF7 | 15 | PRESENT | UNASSIGNED, ORF7, PIIIA | PMID_claude_generated 36016389 | One of the five genus-specific genes. NC_001958 24509..25165, 218 aa; the discovered cluster's median is 225 aa. Function UNKNOWN -- keeping the number. |
| ORF8 | 4 | PRESENT | ORF8 | PMID_claude_generated 36016389 | One of the five genus-specific genes. NC_001958 complement(25201..25698), 165 aa; cluster median 174 aa. Function UNKNOWN -- keeping the number. |
| fiber / sialic acid receptor | 46 | PRESENT | FIBER | PMID_claude_generated 32459612 | TAdV-3 binds sialic acid on N-linked membrane glycoproteins of B cells; required for infection. Fibre head is a twisted beta-sandwich binding sialyllactose. |
| 100K / host shut-off | 44 | PRESENT | K100 | PMID_claude_generated 41960683 | SHUTOFF WAS DECLARED AS A SEPARATE FEATURE AND THAT WAS WRONG. The 8 'shutoff protein' deposits are 703 aa and all consolidated into K100 -- the late 100K protein IS the host shut-off protein. Synonym, not a gene. NC_001958 18184..20226, 680 aa. |
| U exon protein | 37 | PRESENT | UXP | PMID_claude_generated 17881437 | Late l-strand product. |
| protein IX | 0 | ABSENT | — | PMID_claude_generated 33789897 | ABSENT, and decided. Siadenovirus is the one genus with no pIX; NC_001958 annotates none and the cross-genus capsid-cement comparison says so directly. Not a gap. |
