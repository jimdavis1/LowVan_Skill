# Capillovirus — binning decisions

Dump: 779 unique proteins, 924 features, 515 genome records, one genus.
Binned 896 of 924 = **97.0%**.

## Three annotation strings are one protein

In the 280-360 aa band the source offers `movement protein`,
`putative movement protein`, `serine protease` and `36K protein`. An
all-vs-all blastp of one representative each:

    serine protease       vs putative movement protein   96.9%  E=0.0   642 bits
    36K protein           vs putative movement protein   98.4%  E=0.0   648 bits
    36K protein           vs serine protease             96.9%  E=0.0   639 bits
    movement protein      vs 36K protein                 42.5%/212      188 bits

The first three are the same protein. `serine protease` is simply wrong:
capilloviruses encode a papain-like **cysteine** protease as an ORF1 domain,
and nothing at 320 aa is a standalone protease. Binning on the word
"protease" would have created a 25-sequence phantom feature and taken those
sequences out of MP.

## MP is bimodal, and the default window would have deleted half of it

MP has a 320 aa class and a 463 aa class. They are 42% identical over 212
residues — divergent, but the same ORF2 homology group, and both real. The
pooled median is 320, so the default 0.7–1.3x window is 224–416 and all 119
sequences of the 463 class are discarded as overlong.

Set explicitly: **MP 260–480**. Same failure mode as Citrivirus CP, arriving
through the length filter instead of the string rules.

## CP is fused into ORF1, and that is visible in the data

ORF1 is a ~241 kDa polyprotein carrying the coat protein at its C-terminus;
CP is also expressed from a subgenomic RNA. 11 of 56 records whose CP was
annotated inside the polyprotein failed the Met rule — 8 of them beginning at
Ile at exactly 205 aa. They are dropped for profile building, correctly: the
CP profile must start at the subgenomic Met.

## Not modelled

`RNA-dependent RNA polymerase core region` — 26 features, all exactly 270 aa,
a Pfam domain footprint on ORF1, not a separate CDS. Likewise bare
`helicase` / `methyltransferase`. Declaring them would create features no
genome separately encodes.

## Result

| feature | collection | profiles | window (aa) |
|---|---|---|---|
| POLY | 320 | 2 + 3 leftover | 1537–2675 |
| MP | 278 | 2 + 2 leftover | 248–492 |
| CP | 45 | 1 (`-m 2`) | 203–328 |

Routing: 13 references cover 334/334 genomes = **100.0%**.

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
