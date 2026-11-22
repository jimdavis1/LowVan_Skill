# Literature proteins vs BV-BRC calls — Macrobrachium_rosenbergii_nodavirus

Collections are built from deposited strings, so a protein that is
real biology but unannotated in these genomes is invisible to the
build. This records what was probed for and what was found.

| protein | deposits | verdict | binds to | source | note |
|---|---|---|---|---|---|
| protein A (RdRp) | 86 | PRESENT | A, UNASSIGNED | PMID 32659308 | RNA1 replicase, 1033-1045 aa. A methyltransferase and the RdRp separated by a predicted transmembrane region, from a structural model only. |
| B2 | 12 | PRESENT | B2 | PMID 23318596 | 131-133 aa, overlapping the 3' end of protein A in a shifted frame. NO subgenomic RNA has been demonstrated for MrNV, PvNV or CMNV, and RNAi-suppressor activity has NOT been shown for any of the three, so the annotation claims neither. |
| capsid protein | 53 | PRESENT | CP | PMID 30346944 | 371 aa in MrNV, 368 in PvNV, 437 in CMNV. A single uncleaved chain with dimeric P-domain spikes, unlike the trimeric spikes of alpha- and betanodaviruses, and with NO gamma peptide. |
| B1 | 0 | ABSENT | — | PMID 30346944 | DELIBERATELY ABSENT. No B1 is documented for MrNV, PvNV or CMNV and none is annotated in NC_005094.1, NC_014978.1 or MT270124.1. Its absence is the single gene that separates this module from Betanodavirus, so ABSENT here is the expected and required result. |
| capsid maturation cleavage product | 0 | ABSENT | — | PMID 34445426 | DELIBERATELY ABSENT. The one reported MrNV capsid cleavage, at Arg26, was produced by an E. coli cysteine protease in recombinant protein and is blocked by E-64 -- a heterologous-expression artefact, not viral maturation. The virion capsid is a single 43 kDa polypeptide (PMID 15660790). |
