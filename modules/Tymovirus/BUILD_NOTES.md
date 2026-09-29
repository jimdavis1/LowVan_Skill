# Tymovirus — binning decisions

Dump: 250 unique proteins, 342 features. Binned 329 of 342 = **96.2%**.
The cleanest taxon in this group: three overlapping ORFs.

    ORF1  replicase polyprotein   ~1830 aa   MT-Pro-Hel-RdRp
    ORF2  overlapping / movement   ~628 aa   overlaps ORF1, 5'-proximal
    ORF3  coat protein             ~189 aa   3'-proximal, subgenomic

**MP is bound before REP**, because in tymoviruses the movement protein is
defined by overlapping the replicase ORF — any string naming both belongs to
MP (`overlapping protein/movement protein`, `overlapping out-of-phase
protein`).

**TYMVgp1/gp2/gp3 are positional labels** — the U-number trap — and are not
trusted as names. They are admitted only inside the three explicit length
windows, which is what actually decides them: gp1 is 628 aa (MP), gp2 1844
(REP), gp3 189 (CP). A future record numbering them differently is rejected
by the window rather than misfiled.

Source spellings absorbed include two misspellings, `RNA-dependend RNA
polymerase` and `RNA-depedent RNA polymerase`, and the
`unnamed protein product; 206K ORF (AA 1-1844)` form.

| feature | collection | profiles | window (aa) |
|---|---|---|---|
| REP | 77 | 5 + 12 leftover | 1548–2412 |
| MP | 66 | 5 + 11 leftover | 541–731 |
| CP | 69 | 5 + 10 leftover | 176–214 |

Routing: 25 references cover 94/102 genomes = **92.2%**; 8 clusters holding 8
genomes uncovered, to be named in §07 of the coverage audit.

One record, `Piper methysticum tymovirus 1` at 9,491 nt, is 40% longer than a
normal tymovirus genome. It is classified as Tymovirus in BV-BRC and is kept;
it widens the segment window but is not a rep contig.

## The evaluation is a sample, not the whole taxon

Annotation costs ~86 seconds per genome, so the evaluation that decides
`copy_num` and reports per-genome quality was run on a **random sample of up
to 150 genomes per module** (seed 29; Tymovirus has only 102 and so is
complete). Every genome already scored before the cap was kept in the sample
rather than discarded.

150 genomes is ample for the decision it feeds: the test for `copy_num` is
whether a feature is called on at least 85% of routed genomes, and at n=150 a
true rate of 85% is estimated to within about +/-6 points. Capillovirus was
already unambiguous at n=101, all three features calling at 100.0%.

**Whole-taxon coverage is measured separately** by `evaluate_coverage.py`,
which clusters the module's entire BV-BRC holding and annotates one exemplar
per cluster. That, not this sample, is the figure the coverage audit reports
for reach across the taxon.
