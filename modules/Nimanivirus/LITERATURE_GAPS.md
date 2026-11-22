# Literature proteins vs BV-BRC calls — Nimanivirus

Collections are built from deposited strings, so a protein that is
real biology but unannotated in these genomes is invisible to the
build. This records what was probed for and what was found.

| protein | deposits | verdict | binds to | source | note |
|---|---|---|---|---|---|
| ORF1a polyprotein | 33 | PRESENT | PP1A, PP1AB | PMID 35011135 | 4933 aa, invariant. |
| ORF1ab polyprotein | 1 | BELOW_FLOOR | PP1AB | PMID 35011135 | 7797 aa, -1 frameshift at 5'-GGGUUUU-3'. |
| ORF2 nucleocapsid | 32 | PRESENT | N | PMID 35011135 | 125 aa, invariant. |
| ORF3 glycoprotein precursor | 44 | PRESENT | ORF3 | PMID 33023199 | 1503 aa, invariant. |
| p25 triple membrane-spanning protein | 0 | ABSENT | — | PMID 35011135 | DELIBERATELY ABSENT, and this is the reason this module exists separately. The ICTV profile records that MrGV lacks the triple-membrane-spanning p25 at the N-terminus of its precursor glycoprotein. ABSENT here confirms no okavirus cleavage rule has leaked into this module. |
| px / ORFX | 0 | ABSENT | — | PMID 35011135 | DELIBERATELY ABSENT. The ICTV profile states ORFX is absent from MrGV. |
| gp116 / gp64 cleavage products | 0 | ABSENT | — | PMID 33023199 | NOT DECLARED, and recorded as a gap. ICTV asserts MrGV ORF3 encodes homologues of gp116 and gp64, but no cleavage site is published for MrGV and the okavirus sites cannot be transferred -- measured here, MrGV ORF3 against YHV pp3 gives a best HSP of 20.4 bits at 23.7% identity over 38 residues, which is noise. |
| ORF4 | 0 | ABSENT | — | PMID 35011135 | Not reported for MrGV; the ICTV profile does not state whether ORF4 is present. No deposit exists. |
