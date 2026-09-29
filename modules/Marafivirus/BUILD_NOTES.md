# Marafivirus — binning decisions

Dump: 219 unique proteins, 252 features. Binned 235 of 252 = **93.3%**.

Small taxon, thin feature set by design. ORF1 is a ~2,100 aa replicase
polyprotein with the coat protein fused at its C-terminus; CP is also
expressed separately.

**Why not merged with Tymovirus.** Marafiviruses lack the tymovirus
overlapping/movement protein and carry CP inside ORF1, so the gene layout
differs — a module is defined by genome organisation, so they are separate.

**CP is one collection, not two.** The source offers `21 kDa capsid protein`,
`minor coat protein`, `major capsid protein` and `minor capsid protein`. The
p21 major and p31 minor forms are the same reading frame, p31 being the
N-terminally extended product, so they are one homology group; clustering can
separate them into profiles if it wants to. Splitting on the words
"minor"/"major" would have made two features of 4 and 5 sequences, neither of
which can build a profile.

`RNA-dependent RNA polymerase` at 601 aa is an ORF1 domain footprint, not a
separate CDS, and is excluded as in Capillovirus.

| feature | collection | profiles | window (aa) |
|---|---|---|---|
| POLY | 142 | 5 + 7 leftover | 1798–2398 |
| CP | 14 | 3 (`-m 2`) | 183–309 |
| MP | 11 | 2 (`-m 2`) | 250–450 |

Routing: 22 references cover 161/161 genomes = **100.0%**.

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
