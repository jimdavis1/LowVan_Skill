# Literature proteins vs BV-BRC calls — Alphanodavirus

Collections are built from deposited strings, so a protein that is
real biology but unannotated in these genomes is invisible to the
build. This records what was probed for and what was found.

| protein | deposits | verdict | binds to | source | note |
|---|---|---|---|---|---|
| protein A (RdRp) | 16 | PRESENT | A, UNASSIGNED | PMID 7815535 | The replicase, translated from genomic RNA1. Its N-terminal 46 residues are sufficient for mitochondrial localisation and membrane insertion (PMID 12208963); replication occurs on outer mitochondrial membranes (PMID 11689648). No cleavage product of protein A is described in any nodavirus. |
| B1 | 3 | PRESENT | B1 | PMID 12573589 | In-frame C-terminal fragment of protein A, translated from an internal AUG on subgenomic RNA3. NOT a cleavage product. PMID 15163762: 'Protein B1 is in the same reading frame as protein A and thus represents its C-terminal fragment'; detected as protein for Nodamura by SDS-PAGE and MALDI-TOF. |
| B2 | 12 | PRESENT | B2 | PMID 12016316 | RNA silencing suppressor in the +1 frame relative to protein A, on sgRNA3. Nodamura makes two forms, B2-134 and B2-137, from the first and second AUG of RNA3. |
| capsid precursor alpha | 20 | PRESENT | ALPHA, GAMMA, UNASSIGNED, BETA | PMID 6548308 | The RNA2 product, 395-407 aa, which assembles and then cleaves itself. |
| mature capsid beta | 4 | PRESENT | BETA | PMID 1404613 | alpha residues 1 to the Asn of the autocatalytic site; 363 aa in FHV. Derived here by chemistry on each genome's own alpha. |
| gamma peptide | 4 | PRESENT | GAMMA | PMID 9765417 | alpha from the residue after the cut to its C-terminus; 44 aa in FHV, 40 aa in Pariacoto at an Asn/Ser site. Forms a pentameric helical bundle at the fivefold axes and is required for membrane disruption during entry. |
| Pariacoto minor capsid (second-AUG form) | 0 | ABSENT | — | PMID 13679619 | A real virion protein in Pariacoto only, initiated at the SECOND AUG of the capsid ORF and 24 residues shorter than alpha. It is alpha minus an N-terminal run, so it is 100% identical to alpha over its entire length and no profile can separate the two. Deliberately not declared; recorded in NOT_CALLABLE.md. |
| BBV RNA2 second ORF (URF) | 1 | BELOW_FLOOR | UNASSIGNED | PMID 6548308 | A 72 aa open reading frame predicted at base 1110 of BBV RNA2 in 1984, with no product reported in the 40 years since. One deposit family-wide. Excluded by name, not modelled. |
