# Splice-Variants — hand-curated references for `special: splice`

`annotate/get_splice_variant_features.pl` calls a spliced feature by BLASTn
against the references here and then resolving the junction from coordinates
carried in the FASTA header. Without the reference set the feature is
declared and never fires, silently.

Layout mirrors the runtime: `Splice-Variants/<Module>/<FEATURE>.fasta`, which
is what the program looks for.

Header format, coordinates **local to the reference sequence**:

```
>genome_id SD:SD_start-SD_end;last_nt_kept SA:SA_start-SA_end;first_nt_kept
```

`SD_start-SD_end` and `SA_start-SA_end` are windows matched **exactly**
between reference and subject — that is the conservation check. The two
semicolon values are the last nucleotide kept on the 5' side and the first
kept on the 3' side, so the product is `seq[..sd_last] + seq[sa_first..]`.

**The reference must be the gene region, starting at its own start codon, not
the whole genome.** The program builds the 5' piece from the start of the
*aligned* region and rejects any that translates with a stop, so a
whole-genome reference makes that piece everything upstream of the gene and
the splice is skipped with only a warning on stderr.

These are hand-curated and were previously only in the runtime directory,
which has no remote. They are the input the programs cannot work without.
