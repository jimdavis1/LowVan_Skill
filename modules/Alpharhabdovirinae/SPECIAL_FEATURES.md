# MERHA_L_SPLICED has no reference set

`Alpharhabdovirinae` declares one feature with `"special": "splice"`:
`MERHA_L_SPLICED`. A spliced feature is called by
`get_splice_variant_features.pl` from a curated **nucleotide** reference set
under `Splice-Variants/<Module>/<FEAT>.fasta`, not by a PSSM.

**That reference set was never built.** `Splice-Variants/` in the runtime
checkout holds only the four influenza genera. So the feature is declared and
never fires — which is the legitimate "declared as a marker that the taxon is
known and unbuilt" case from `references/module-partitioning.md`, but it is
only legitimate while it is written down. This file is that record.

To build it, see `references/special-features.md` §"Build the reference set":
one reference per genome rather than one per clade, no ambiguous bases, and a
self round-trip check that every donor genome recovers its own protein
byte-identically.
