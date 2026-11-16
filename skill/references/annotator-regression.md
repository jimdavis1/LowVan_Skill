# Re-testing an annotator change

`annotate_by_viral_pssm.pl` is shared by every LowVan module. A change made to
fix one family is a change to all of them, and the family you were fixing is
not the one most likely to break.

**Curator's rule (12 Nov 2026):** any change to the annotator is run against a
fixed genome set, old against new, and every difference is explained before
the change is committed. Keep the set and these instructions with the kit, so
the next change is tested the same way.

## What happened the first time

The adenovirus build needed the annotator to prefer a profile whose hit can be
delivered without an internal stop. The first version counted *any* stop in an
HSP as a defect, including the CDS's own terminal stop, which is how a full
coding sequence is called when `downstream_ext` is off. On the family it was
meant to fix it made E4ORF5 137 → 79 aa, K100 695 → 616 aa, E1A 193 → 140 aa,
and many features lost their last residue. The second version tests the crop's
own rule (an HSP is deliverable if it has no stop, or its pre-stop length still
clears the coverage cutoff) and changed nothing except the calls it was meant
to fix. Only running both versions over the same genomes showed the
difference.

## The set

`skill/assets/annotator_regression/genomes.tsv`, with the contigs in
`contigs/` beside it (about 1 MB):

- **one complete genome per installed non-adenovirus module**, every segment
  included for segmented taxa. These were each module's first representative
  contig on 14 Nov 2026, copied so the set does not drift with a data release.
- **11 adenovirus genomes** that the 12 Nov ranking change either altered
  (E4ORF1, PV, IVA2, PVII, E4ORF5, I_LEADER, E1A, PIX) or had to leave alone.

Add a genome when a change is made for it; never remove one to make a run
pass.

## Running it

```bash
export LOWVAN_DATA_DIR=/path/to/Viral_Annotation      # the SAME data for both
python3 skill/scripts/annotator_regression.py \
    --old $LOWVAN_DATA_DIR/annotate_by_viral_pssm.pl.pre_change \
    --new annotate/annotate_by_viral_pssm.pl \
    --out regression_run/ --jobs 3 --threads 2
```

About 30 minutes on a laptop. Outputs are reused, so an interrupted run
resumes; delete `<out>` to start clean.

## Reading the report

`<out>/REPORT.md` lists, per genome, the feature count and every feature whose
calls differ, with old and new (start, end, strand, length in aa, profile).

- **0 differences outside the genomes the change targeted** is the result you
  want.
- **A difference is a finding, not a failure.** Explain each one: what the old
  call was, what the new one is, and why the new one is right. Put that in the
  commit message.
- **A genome with no output** from one version is a crash: read `<out>/*.err`.
- A change to *data* (profiles, cutoffs, `Viral_PSSM.json`) is not tested
  here. Both versions read the same `LOWVAN_DATA_DIR` on purpose, so the only
  variable is the code.
