# Literature proteins vs BV-BRC calls — Okavirus

Collections are built from deposited strings, so a protein that is
real biology but unannotated in these genomes is invisible to the
build. This records what was probed for and what was found.

| protein | deposits | verdict | binds to | source | note |
|---|---|---|---|---|---|
| pp1a replicase polyprotein | 108 | PRESENT | PP1AB, PP1A | PMID 12376758 | 4060 aa GAV / 4070-4076 aa YHV; carries 3CLpro, PLP1, PLPx and HD1-HD4. |
| pp1ab replicase polyprotein | 106 | PRESENT | PP1AB | PMID 10811931 | 6673-6692 aa, from a -1 ribosomal frameshift at 23-25% efficiency. |
| p20 nucleocapsid | 20 | PRESENT | N | PMID 15280504 | 144 aa GAV / 146 aa YHV, ORF2, upstream of the glycoprotein gene -- unique among nidoviruses. |
| pp3 glycoprotein precursor | 53 | PRESENT | PP3, NOT_MODELLED | PMID 12655087 | 1640-1666 aa, from subgenomic mRNA 3, cleaved into p25 + gp116 + gp64. |
| gp116 envelope glycoprotein | 4 | RECOGNISED | NOT_MODELLED | PMID 12655087 | 845-899 aa. Deposits exist but are 267-338 aa FRAGMENTS, so the real product is derived by cutting pp3 at its measured [Ala-X-Ala] site rather than binned. |
| gp64 envelope glycoprotein | 3 | RECOGNISED | NOT_MODELLED | PMID 12655087 | 538-539 aa. Same: the deposits are fragments; the product is derived. |
| p25 triple membrane-spanning protein | 0 | ABSENT | — | PMID 12655087 | 228 aa, the N-terminal pp3 product. ABSENT from the database: no deposit anywhere. Derived from pp3. NOTE the source disagreement -- the ICTV profile calls it p25 and the Walker & Sittidilokratna review calls it p22 at 228 aa; 228 aa is about 25 kDa, so p25 is used. |
| px / ORFX | 0 | ABSENT | — | PMID 20017924 | 86 aa, alternative frame overlapping the 5' end of ORF3. ABSENT from the database: zero deposits. Extracted from the genomes; all 15 extractions exactly 86 aa. |
| ORF4 | 2 | BELOW_FLOOR | ORF4 | PMID 12376758 | 83 aa, GAV only; truncated after 20 codons in YHV. |
| nsp mature peptides of pp1ab | 4 | RECOGNISED | NOT_MODELLED | PMID 12502857 | DELIBERATELY NOT DECLARED. The ICTV profile states the number and composition of the nsps released from pp1ab is not yet known. Deposits annotated 3C-like protease, helicase and RNA polymerase exist but are domain fragments with no published junctions; only the 3CLpro SPECIFICITY VxHE/(L/V) is characterised, and a specificity is not a coordinate. |
| membrane (M) protein | 0 | ABSENT | — | PMID 12376758 | DELIBERATELY ABSENT. Roniviruses have no discrete M gene; p25 is the protein that resembles a coronavirus M. |
