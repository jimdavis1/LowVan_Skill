
## Captured 30 September 2026, before a machine restart

Three pipeline edits were live in a working copy under `/private/tmp` and in no
patch here. They are now captured. Two of them are the same defect.

- **`viral_genome_quality.pl.patch`** — `$genome_in->{features}->[0] or die`
  killed the scorer on any genome with zero features. A genome that called
  nothing is not unscoreable; it is the worst case and the one the "no features
  called" bucket exists for. **Five of 97 Tobamovirus panel genomes hit it**,
  were removed from the scored set, and the summary then reported
  *"no features called: 0"* while five genomes had called none. Now scored as
  Poor with every essential correctly flagged absent.
- **`annotate_by_viral_pssm-GTO.pl.patch`** — the other half of the same bug.
  `annotate_by_viral_pssm.pl` emits a `no_features_called` sentinel row so the
  taxon assignment is not lost; that row has no coordinates, so `$feature`
  stayed undef and pushing it undef made GenomeTypeObject die with
  "No feature location", losing the whole genome. Now skipped.
- **`get_splice_variant_features.pl.patch`** — retires the fragments a spliced
  feature was assembled from. Adding is the right default and is why this did
  not exist: in influenza every spliced product has a *different* annotation
  from its unspliced parent (M1/M2, NS1/NS2, PA/PA-X) and both are real.
  Merhavirus is the other case — CTRV's L gene carries a 76-nt intron, so the
  ordinary L profile hits two exons as two fragments, 1333 aa and 801 aa, and
  neither is a protein; only the 2123 aa join is. Leaving all three gave one
  genome three features with one annotation, which `viral_genome_quality.pl`
  counts against `copy_num` and reports as "too many HSPs". A fragment is
  recognised by sharing an exon boundary with the join *and* carrying the same
  function.

Each was verified to reverse cleanly against the working tree it came from.
