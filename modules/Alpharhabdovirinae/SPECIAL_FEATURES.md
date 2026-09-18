# MERHA_L_SPLICED — moved to the Merhavirus module, and built

This module no longer declares a `special` feature. `MERHA_L_SPLICED` was
deferred here for a long time as "real, buildable, not yet built"; it is now
built, in a module of its own.

**Why it could not work here.** `get_splice_variant_features.pl` reads its
hand-curated references from `Splice-Variants/<Module>/<FEAT>.fasta`, keyed on
the module. A reference set for CTRV's L would have had to sit beside a module
covering ~38 genera. Genus *Merhavirus* is now its own module, the same split
as `Orthopneumovirus_muris` from `Orthopneumovirus`.

**What moved.** `Alpharhabdovirinae.8.dna` (NC_025384.1, Culex
tritaeniorhynchus rhabdovirus) is removed from this module's reference set, so
the genus routes to the module that can call its L. This module keeps 9 rep
contigs.

**What did not move.** The profiles. 16 clusters across G, G_MAT, G_SP, L, M,
N and P contain Merhavirus sequences and were copied — not moved — into
Merhavirus and regenerated under that module's name. Those sequences were part
of what these profiles were built from, and removing them would change this
module's calls on everything else.

See `modules/Merhavirus/SPLICE_NOTES.md` for the intron, the reference set and
the fragment-retirement behaviour.
