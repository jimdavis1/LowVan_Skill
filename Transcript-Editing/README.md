# Transcript-Editing — hand-curated references for `special: transcript_edit`

`annotate/get_transcript_edited_features.pl` calls ribosomal-frameshift and
RNA-editing products from the references here. Same contract as
`Splice-Variants/`: no reference set, no call, no error.

Layout mirrors the runtime: `Transcript-Editing/<Module>/<FEATURE>.fasta`.

## Case marks the inserted bases

A reference here is **not** collinear with the genome: it is the genome with
one or two bases duplicated at the ribosomal slip site, or inserted at an
RNA-editing site, so that two reading frames translate as one. Which bases
those are is the whole content of the curation, and in an unmarked file they
are unrecoverable without going back to the genome.

So: **genome-derived bases lowercase, the added bases UPPERCASE.** Every
sequence in an insertion-type set carries exactly one or two uppercase
characters, and `grep -o '[ACGT]'` finds them.

```
Torovirinae/ORF1AB.fasta    ...gatttaaacCtgt...   the UUUAAAC slip, C duplicated
Orthoavulavirus/W.fasta     ...tgctaaaaaGGgg...   two G's added at the editing site
```

Nothing in the pipeline reads case. `get_transcript_edited_features.pl` builds
the edited sequence from blastn's own `qseq`/`sseq`, and blastn returns those
uppercase whatever the input file holds — measured, not assumed: a marked
Torovirinae reference blasted against its own genome returns a `qseq` with zero
lowercase characters. The twenty hand-curated sets have been mixed-case in
production all along. That the mark changes nothing is exactly why a program
has to write it: no test will ever fail because it is missing.
`build_transcript_edit_refs.py` and `special_cardio.py` apply it, and
`skill/scripts/mark_transcript_edit_inserts.py` repairs a set written before
they did, re-deriving each insertion from the genome (or from the deposited
pre-shift protein) and refusing any sequence it cannot confirm.

A set whose edit REMOVES a base has nothing to mark and stays lowercase
throughout: the Closteroviridae RdRp sets are of that kind, because a +1
frameshift skips a base rather than adding one.

The mark is also a check on the builder. In the two Astroviridae sets, built
before `--heptamer` existed, the uppercase base sits 2–19 nt away from the
`AAAAAAC` slippery site in most references — a valid insertion point that
reproduces the deposited protein but is not where the ribosome slips. That is
obvious at a glance in a marked file and invisible in an unmarked one.

## Modules

18 modules depend on this — the coronaviruses (ORF1ab, NSP12), the
paramyxoviruses (V, W, P, I), Orthoebolavirus (GPC, Gn, ssGP) and
Togaviridae (TF). These are hand-curated and were previously only in the
runtime directory, which has no remote.
