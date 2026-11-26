# Special features: readthrough, frameshift and splicing

Most viral proteins are collinear with the genome: one interval, one reading
frame, one PSSM. This file is about the ones that are not, and about the one
JSON flag that handles a protein which *is* collinear but contains a stop.

Three mechanisms, three different answers:

| the genome does this | the protein is | how to handle it |
|---|---|---|
| in-frame stop read through ~10% of the time | still collinear, one frame | `internal_stop: 1` on the feature, normal PSSM |
| one codon is edited, changing what it means | still collinear, one frame | `transl_except`, normal PSSM |
| ribosome shifts frame mid-gene | **not** collinear | `special: transcript_edit`, no PSSM |
| polymerase inserts non-templated bases | **not** collinear | `special: transcript_edit`, no PSSM |
| introns removed from the transcript | **not** collinear | `special: splice`, no PSSM |

Getting the first one wrong costs you a truncated protein on most of the taxon.
Getting the second wrong costs you the protein entirely, because no profile can
ever call it.

The second row is the one that looks like the third and is not. A substitution
edit changes a codon's *meaning* and adds no bases, so the protein stays
collinear and a PSSM calls it perfectly well. `special: transcript_edit` cannot
do it: its gap-filling loop takes the genome's own base at every ungapped
column, so an edited codon survives into the translation unchanged, and a
substitution that creates or removes a stop then fails the caller's
internal-stop check outright.

---

## Part 1 — Readthrough: `internal_stop`

### The symptom

A mature peptide comes out consistently short, by a small and *constant* number
of residues, across most of the taxon. Alphavirus nsP3 came out 7 residues short
on 8 genomes in 9.

### The cause

An in-frame stop that the ribosome reads through part of the time. In
alphaviruses it is an opal (UGA) sitting **six codons before the end of nsP3** —
inside the protein, not between two proteins. About 10% of ribosomes read
through it and make the full polyprotein; the rest terminate there.

Two consequences, and it is easy to conflate them:

- **tblastn breaks at a stop by design**, returning two HSPs instead of one, so
  the profile match stops at the opal and the call is cropped.
- **Both products are real.** The terminated form and the read-through form are
  different proteins with different jobs — in alphaviruses, P123+nsP4 makes
  minus-strand RNA and the fully cleaved set makes plus-strand.

### The fix

Declare `internal_stop: 1` on every feature whose match may legitimately contain
one, and leave it off everywhere else:

```json
"NSP3":     { "internal_stop": 1, ... },
"NSP1234":  { "internal_stop": 1, ... }
```

**Omit the flag rather than writing `internal_stop: 0`.** Then
`grep internal_stop` returns exactly the features that permit one, which is the
question an auditor actually asks. Writing the zeros makes that grep return
everything.

**The collection needs `-x`, and `build_features.py` now adds it for you when
the feature declares `internal_stop`.** Where the read-through codon is
deposited as an `X` — which is how Coltivirus segment 9 arrives, in 17 of 17
members — `fasta-cluster-pssm-2.pl` drops every one of them, because by default
it discards any sequence matching `/[BJXZ]/`. The whole collection then empties,
and what you see is not a message about the feature but mmseqs dying on a
0-byte temp file:

```
Cannot seek at the end of file 2c00cafae6.fasta
VP9RDT  FAILED -- no output directory carries VP9RDT profiles
```

which reads like an mmseqs fault. Tobamovirus got this right by hand, setting
`"flags": ["-x"]` on REP126 and REP183; nothing made the next module do the
same. It is now derived from the declaration, so a feature that permits a
read-through keeps the sequences that demonstrate one.

`-mis N` (default 1) caps how many stops are tolerated. One opal is expected; a
run of stops is a broken genome, and the annotator crops rather than inventing a
readthrough across it.

### Build the profiles so they span the stop

The flag alone is not enough. A profile trained on sequences that all stop at
the opal has never seen what comes after it and will not extend past it. Rebuild
the feature's alignments from **the end of the preceding protein to the start of
the following one**, so the readthrough is inside the training data. For
alphavirus nsP3 that meant nsP2-end to nsP4-start.

Two traps when you do this:

- **mafft rewrites `*` to `-`, and psiblast refuses an MSA containing `*`.**
  Write the readthrough position as `X` in the build FASTA.
- Some pipelines drop sequences matching `/(B|J|X|Z)/i`. `fasta-cluster-pssm-2.pl`
  needs `-x` to keep them (`FCP_Main_Utils.pm:244`).

### Verify behaviourally, not structurally

Byte-identical profiles prove nothing about whether readthrough works. Call a
genome and look at the protein:

```
nsP3  556 aa   ...RRRRRSRRTEY*LTGVGG
nsP4 starts    YIFSTDTG...
```

The stop should sit the documented distance from the C-terminus (six residues
here), on every genome that has one, and be **absent** on the genomes that lack
the stop entirely. A single exemplar tells you nothing; check both classes.

### A precursor that terminates at the stop is a different feature

If the terminated form is itself a real protein, declare it separately — but its
C-terminus is **the stop**, not the full-length product's C-terminus. Alphavirus
P123 ends 18 nt short of nsP3's end and that is correct, not a truncation: P123
*is* the opal-terminated product, while nsP3 is the readthrough mature peptide
that runs six codons further. Two proteins, two C-termini.

### When the match carries more stops than it may read through

A readthrough match can pick up a second stop, and the usual cause is not a
broken match. If the profile is **longer than this genome's protein**, the
alignment overruns the real terminator, so the HSP ends up holding the
readthrough stop *and* the terminator, and `n_stops` exceeds `-mis`.

The annotator used to call that broken and crop the gene at the HSP edge, which
put the C-terminus at the **first** stop — the readthrough site — and so
returned the short product under the long product's name. Odontoglossum
ringspot virus (BV-BRC 12238.214) has its amber at residue 1112 and its
terminator at 1596, against a 1616-residue Tobamovirus REP183 profile: it was
called at 1112 aa and flagged "Feature is too short". Raising `-mis` to 2 did
not help either — it read straight through the real terminator and gave 1610 aa,
stopping at the next stop in the 3' UTR.

`internal_stop : 1` means *read through one stop and end at the next*, so
`annotate_by_viral_pssm.pl` now does exactly that: it keeps `-mis` stops and
terminates the gene at the one after them, including the stop codon by the same
convention `scan_to_stop_codon` uses. ORSV now comes out at **1596 aa**, its
real terminator. Genomes whose match carries at most `-mis` stops are
unaffected — the path is not reached — and CHIKV S27 and a control tobamovirus
both annotate byte-identically before and after.

**It is gated on `min_len`, and it has to be.** A genome carrying frameshifting
indels has stops scattered from early on, and terminating at the second of
those returns a fragment rather than a protein. The first version of this fix
had no gate and produced a 23-residue "RNA-dependent RNA polymerase" on
Tobacco mild green mosaic virus (12241.68) and a 167-residue one on Bottle
gourd mottle virus (2034161.3) — it turned two already-flagged genomes into
worse-flagged genomes and cost more than it gained. The rule is now: take the
read-through termination only if the result reaches the feature's declared
`min_len`; below that the match really is broken and the old behaviour is the
honest answer. Get this wrong in the permissive direction and the module
invents short features on exactly the genomes least able to withstand it.

Two things this does **not** do. It does not rescue a genuinely broken match —
that is what the `min_len` gate refuses. And it does not stop the terminated
precursor and the full-length product from both being called on the same
genome — that is a separate `copy_num` question, and the two are different
features with different C-termini anyway.

---

## Part 1b — Substitution edits: `transl_except`

### The symptom

A protein the literature describes is never emitted correctly, and the reason
is a single residue. Either the translation carries a `*` in the middle of it,
or an `X`, and which one you get depends on how the deposit happened to be
sequenced rather than on any difference in the virus.

### The worked case

Hepatitis delta virus makes two proteins from one ORF. The small antigen ends
at an amber codon; ADAR1 edits that adenosine on the antigenome, UAG becomes
UGG, and translation runs on for 19 more codons into the large antigen — a
different protein, carrying the farnesylation box that the small one lacks and
that virion assembly needs. Both are real, both are made by every genome, and
neither is optional.

What a genome in the database holds at that position is not a property of the
virus but of the sequencing:

| the deposit holds | translates as | what it means |
|---|---|---|
| `TAG` | `*` | the unedited genome |
| `TRG` | `X` | a mixed population — direct evidence of editing |
| `TGG` | `W` | the edited genome |

Measured across 1,043 HDV genomes: 508, 392 and 143 respectively. The protein's
residue is `W` in all three.

### The fix

```json
"L-HDAg": {
  "internal_stop": 1,
  "transl_except": { "replace_residue": ["X", "*"],
                     "insert_residue": "W",
                     "offset": [19, 20],
                     "from": "C-term" }
}
```

Read it as: within 19 to 20 residues of the protein's C-terminus, a residue
that translated as `X` or `*` is really `W`. `from` takes `C-term` or `N-term`.

### Four things that are easy to get wrong

**The offset is measured on the emitted protein, never on the PSSM query.** A
column recorded against a profile is valid only for the build that produced it
and silently wrong after the next reclustering, and it would put PSSM
filenames into the curator's JSON, which knows only feature names. Measuring
from the protein's own terminus makes the declaration independent of how the
collection happened to cluster.

**Use a range, not a number, and make it tight.** HDV's site sits 19 residues
from the C-terminus in 1,028 deposits and 20 in 42 — genotype III, whose
extension is one residue longer and which has a different editing RNA
structure behind it (PMID 12097551). Both are real. But widen the window
further and ambiguity arrives: at offsets 19–20, *no* deposit of 1,078 holds
more than one `X` or `*`; at 18–21, eight do; at 15–20, twenty-four do. The
window is what makes the rule deterministic, so set it from the data and no
wider.

**It is independent of `internal_stop` and must stay so.** Only the unedited
genome has a stop to read through; a `TRG` genome has none and reaches its
terminator by the ordinary path. Applying the substitution inside the
readthrough branch would fire on half the dump and leave an `X` in the other
half's product, at correct coordinates, with nothing reporting a problem.

**Exactly one residue in the window must qualify.** If two do — an N-masked
genome — the honest answer is to substitute nothing and log it, because
choosing between them is guessing. The feature is still called; it carries
what the genome actually encodes.

### The other half of the problem

The *short* form ends at the editing site. On an unedited genome a stop marks
that; on the other half of the dump nothing does. So the short form needs
`downstream_ext: 0`: scanning to the next stop would run it on to the long
form's terminator and emit the long protein under the short one's name, two
features on identical coordinates. An end that no codon marks must not be
scanned for. The two features also need declaring in `protected_pairs.json`,
since one is a strict prefix of the other by biology.

---

## Part 2 — Frameshift and editing: `special: transcript_edit`

### First, know which biology you have

`transcript_edit` is the kit's **mechanism**, not a claim about the virus. It
serves two genuinely different things:

| | paramyxovirus V/W | alphavirus TF, coronavirus ORF1ab |
|---|---|---|
| what happens | polymerase inserts a **non-templated** base | ribosome re-reads a **templated** base |
| mRNA differs from template? | **yes** — a real edit | **no** — identical |
| would you find it by sequencing mRNA? | yes | no |

They are informationally identical — a +1 templated duplication and a −1
ribosomal slip produce the same reading — which is why one mechanism serves both.
But do not write "the transcript is edited" into an annotation for a frameshift
product. There is no edited transcript to find.

### Diagnosing non-collinearity

tblastn the protein against its own genome. A frameshift product returns **two
HSPs in different frames**, overlapping by a few residues because each extends
through the junction on mismatches:

```
q(aa)    s(nt)        ident   frame
1-50     9772-9921    96.0%   1      <- 0-frame portion
46-76    9906-9998    93.5%   3      <- shifted portion
```

**The frame pair tells you the direction.** On the plus strand BLAST numbers
frames 1/2/3, so −1 reads as 1→3 by wraparound and +1 reads as 1→2. You can
diagnose direction without looking at a base.

A collinear protein returns one clean HSP. If yours does, you do not need any of
this.

### Build the reference set

There is no PSSM. `get_transcript_edited_features.pl` blastn's a curated
**nucleotide** reference carrying the extra base against the genome, fills the
resulting gap in the *subject* from the query, and translates. Everything except
that one base comes from the genome, so each genome gets its own protein rather
than a copy of the reference.

References live at `<LOWVAN_DATA_DIR>/Transcript-Editing/<Module>/<FEAT>.fasta`.

Constructing one, for a −1 frameshift at a slippery heptamer:

1. **Find the slip site, at the right phase.** A slippery heptamer is
   `X_XXY_YYZ` — the lone X is the **third** position of a codon, so the motif
   begins at **phase 2** of the reading frame, not phase 0. Testing for phase 0
   will reject every genome and look like a finding. It is not.
2. **Anchor the search to a window, not to the called feature.** The site sits at
   a near-constant offset (alphavirus TF: codon 39–46, with 190 of 279 at codon
   46 exactly). Searching only inside the called interval misses it whenever an
   N-run truncated that call, and picks up spurious upstream matches.
3. Read the 0-frame up to the slip, **duplicate the slipped base**, continue in
   the shifted frame to the stop.

   **Cite the evidence for that junction, and state how strong it is.**
   Duplicating the slipped base gives the tandem-slippage product: both P-site and
   A-site codons are decoded, so `U_UUA_AAC` reads ...Leu-Asn.... That was shown by
   mass spectrometry of the **SARS-CoV** product (Baranov 2005, PMID 15680415).
   For any other virus with the same heptamer it is an inference from a shared
   mechanism, and for a different heptamer — Serpentovirinae's `A_AAA_AAC` —
   it rests on the general tandem-slippage rule alone. Say which of the three
   it is in the build notes and on the pages. No tobanivirus pp1ab product has
   ever been sequenced, so all four Tobaniviridae modules ship the inferred
   junction and say so. Where submitters' own joins disagree, that is a
   different convention, not a defect in their records; a GenBank or RefSeq
   join is never the evidence.
4. **Decide on the stop codon** — see below.
5. Reject any construct containing an ambiguous base. A reference with an N in
   it is worse than no reference.

### The trailing `*` convention

`gjoseqlib::translate_seq` emits `*` only for an unambiguous stop codon (an
ambiguous one gives `x`) and does **not** strip a terminal stop. The caller
stores the translation verbatim. So the reference decides:

- The product really terminates at a stop → **keep the stop codon** in the
  reference; the protein carries a trailing `*`. (V, W, ORF1ab, ssGP, TF)
- The product ends at a protease or signal-peptidase cleavage → **no stop codon**;
  no `*`. (NSP12, Gn)

Check the existing sets before inventing a convention: they already encode this
correctly, and the split tracks biology, not sloppiness.

Note the star reaches the GTO and the feature table but not the protein FASTA —
`annotate_by_viral_pssm.pl` strips it for `@prot_seqs` only.

### CDS or mat_peptide?

By where the **N-terminus** comes from, not by how exotic the protein is. TF has
no initiating Met — its N-terminus is made by signal peptidase at the E2/6K
junction, the same cut that makes 6K's — so it is a `mat_peptide` and its
reference starts `GCC`, exactly as coronavirus NSP12's start `TCA`/`AGT`. V and W
initiate at their own Met and are `CDS`.

A product can be a hybrid: TF has a cleaved N-terminus like NSP12 and a real stop
codon like V/W. Take each end on its own terms.

### How many references?

**One per genome, not one per clade.** The correction gate is 95% identity with
`-word_size 28`, which demands a close relative. The shipped sets are dense for
this reason — Betacoronavirus NSP12 has 3006, Orthoavulavirus V has 923.
Togaviridae TF has 234 from 191 strains.

Genomes too divergent to correct degrade gracefully rather than failing: at
80–95% identity the annotator emits a `partial_cds` annotated
`Uncorrected <product> encoding region`, with no protein translation.

### Validate with a self round-trip

Every genome that donated a reference must recover its own protein
**byte-identically** through the annotator's own gate-and-fill logic. Anything
less means the construction and the caller disagree. Then run the real pipeline
on a species-spanning panel and check the protein, not just the presence.

---

## Part 3 — What will not tell you the truth

Three tools are blind to special features, in three different ways. None of them
is broken; all of them will quietly mislead you.

- **`validate_calls.py`** runs `annotate_by_viral_pssm.pl` only. It now says so,
  but its feature count is still lower than the pipeline's.
- **`evaluate_coverage.py`** annotates via the flat feature table, which has no
  representation for a special feature. It now warns. Coverage for those
  features is not low in that report — it is absent.
- **`evaluate_module.py`** compares **CDS** features only, and matches one-to-one
  by call index. A `mat_peptide` cannot match a reference CDS however correct it
  is, and where several reference CDS overlap one call they collapse to a single
  entry — so its "matched N (x%)" is not a coverage rate. Read the **missed**
  column instead.

Only `New-annotate-viral-taxon.pl` runs the whole chain.

## Part 4 — Two measurement traps

**Do not blast many queries against one combined database** to simulate the
annotator. `-max_target_seqs` caps reported subjects per query, so genome ×
reference pairs go unreported and the result undercounts — by 236 vs 293 on
Togaviridae TF. The annotator blasts one genome at a time.

**The annotator's gates describe the nucleotide match, not the protein.** It
screens identity, coverage, gap count and gap runs, then (as of this revision)
also screens the translation for ambiguous residues and internal stops, demoting
a failure to `partial_cds`. Before that screen existed, 31 of 293 TF calls
shipped with an `X` in them — because the gap-fill draws its bases from the
*subject*, so an N-masked target produces a bad protein no amount of reference
curation can prevent. If you are counting coverage, count clean calls separately.
