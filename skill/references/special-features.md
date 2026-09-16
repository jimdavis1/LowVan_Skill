# Special features: readthrough, frameshift and splicing

Most viral proteins are collinear with the genome: one interval, one reading
frame, one PSSM. This file is about the ones that are not, and about the one
JSON flag that handles a protein which *is* collinear but contains a stop.

Three mechanisms, three different answers:

| the genome does this | the protein is | how to handle it |
|---|---|---|
| in-frame stop read through ~10% of the time | still collinear, one frame | `internal_stop: 1` on the feature, normal PSSM |
| ribosome shifts frame mid-gene | **not** collinear | `special: transcript_edit`, no PSSM |
| polymerase inserts non-templated bases | **not** collinear | `special: transcript_edit`, no PSSM |
| introns removed from the transcript | **not** collinear | `special: splice`, no PSSM |

Getting the first one wrong costs you a truncated protein on most of the taxon.
Getting the second wrong costs you the protein entirely, because no profile can
ever call it.

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

Expect these two to collide on the occasional genome. Where the full-length
match happens to include the ORF's own terminator as well as the internal stop,
`n_stops` exceeds `-mis`, readthrough is refused, and the long feature is cropped
back onto the short one's coordinates — a duplicate. Raising `-mis` does not fix
it; it reads through the real terminator instead. Treat as an outlier unless it
is common in your taxon.

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
