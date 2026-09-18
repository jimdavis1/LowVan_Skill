# Transcript-Editing — hand-curated references for `special: transcript_edit`

`annotate/get_transcript_edited_features.pl` calls ribosomal-frameshift and
RNA-editing products from the references here. Same contract as
`Splice-Variants/`: no reference set, no call, no error.

Layout mirrors the runtime: `Transcript-Editing/<Module>/<FEATURE>.fasta`.

18 modules depend on this — the coronaviruses (ORF1ab, NSP12), the
paramyxoviruses (V, W, P, I), Orthoebolavirus (GPC, Gn, ssGP) and
Togaviridae (TF). These are hand-curated and were previously only in the
runtime directory, which has no remote.
