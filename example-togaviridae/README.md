# Worked example: Togaviridae

The second worked example, and the one to read if your taxon has a protein that
**a PSSM cannot call**. Rhabdoviridae covers the ordinary path; everything novel
here is about proteins that are not collinear with the genome.

```
build_collections.py          source annotation strings -> feature collections
build_togaviridae_json.py     collections + curated alignments -> the module JSON
build_tf_refs.py              a -1 frameshift product -> a transcript-edited reference set
extract_combos.py             recover a precursor from the interval between its flanking calls
extract_nsp3.py               rebuild a feature whose profiles stop at a readthrough codon
run_bespoke.sh                build the three features that do not come from the normal pipeline
batch_query.pl                the BV-BRC export (one batched request, not one per genome)
```

`build_collections.py` and `build_togaviridae_json.py` are **taxon-specific by
design** — copy them, replace the rule tables at the top, leave the machinery
below alone. The same split as Rhabdoviridae: in `build_collections.py`
everything above `def compile_rules` is Togaviridae, everything below is general.

## What this example adds over Rhabdoviridae

### A protein no profile can call

Alphavirus **TF** is made by a −1 ribosomal frameshift at a conserved UUUUUUA
heptamer inside 6K. tblastn splits it into two HSPs in frames 1 and 3, so no
PSSM will ever call it. It is declared `special: transcript_edit` with no
profile and reconstructed by `get_transcript_edited_features.pl` from a curated
nucleotide reference set that `build_tf_refs.py` builds.

Three things in that script cost a debugging pass each and are commented in
place:

- The slippery heptamer sits at **phase 2** of the reading frame, not phase 0 —
  `X_XXY_YYZ` puts the lone X at the third position of a codon. Testing for
  phase 0 rejects every genome and looks like a finding.
- Anchor the search to a **window** (`WINDOW = (100, 170)`), not to the called
  feature. The site is at a near-constant offset; searching the called interval
  misses it whenever an N-run truncated that call.
- **Keep the stop codon** in the reference, because TF really terminates at one.
  A product that ends at a cleavage instead takes no stop codon and no trailing
  `*`. The convention is per-product, not per-feature-type.

### A readthrough codon inside a protein

Most alphaviruses carry an in-frame opal six codons before the end of nsP3.
`internal_stop: 1` tells the annotator the match may legitimately span it, but
the flag alone is not enough: profiles trained on sequences that all stop at the
opal have never seen what follows. `extract_nsp3.py` rebuilds that feature from
the end of nsP2 to the start of nsP4 so the readthrough is inside the training
data. Note the readthrough position is written `X`, because mafft rewrites `*`
to `-` and psiblast refuses an MSA containing `*`.

### Precursors the source cannot supply

pE2 (E3+E2) and P123 are real proteins with eight and two usable sequences in
the whole export. `extract_combos.py` recovers them from the genome interval
between their flanking calls instead. Note `need_met` differs between them —
P123 starts at the polyprotein's own initiating Met, pE2 starts at a cleavage.

## Read the comments

As with Rhabdoviridae, much of the value is in the comments recording why a rule
exists and what measurement prompted it: why the nonstructural polyprotein is
annotated `Nonstructural polyprotein` and not `... P1234` (naming it P1234 would
be wrong for 86% of the genus), why 227 genomes BV-BRC files under Togaviridae
are excluded as not togaviruses at all, and why the P123 feature legitimately
ends 18 nt short of nsP3's end.

The general treatment of all of this is in
`skill/references/special-features.md`.
