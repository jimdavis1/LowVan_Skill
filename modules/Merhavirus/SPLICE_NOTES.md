# Merhavirus L_SPLICED — built and working

Culex tritaeniorhynchus rhabdovirus (CTRV) carries a **76-nt GU-AG intron** in
its L coding region (Kuwata et al. 2011, *J Virol*, PMID 21507977). CTRV
replicates in the nucleus, which is what makes splicing possible at all —
unique among Mononegavirales outside Bornaviridae.

## Why this needed its own module

`get_splice_variant_features.pl` is keyed on the module: it reads
`Splice-Variants/<Module>/<FEAT>.fasta`. While the feature lived in
`Alpharhabdovirinae`, a reference set for it would have had to sit beside a
module covering ~38 genera. Splitting the genus out is the same move as
`Orthopneumovirus_muris` out of `Orthopneumovirus`.

Seven features are inherited from the parent unchanged, built from the
parent's own clusters that contain Merhavirus sequences — 16 of them across
G, G_MAT, G_SP, L, M, N and P — regenerated under this module's name. **The
parent keeps its copies**: Merhavirus sequences were part of what those
profiles were built from, and removing them would change the parent's calls
on everything else. `Alpharhabdovirinae.8.dna` (NC_025384.1, CTRV) *is*
removed from the parent's reference set, or the genus routes to a module that
cannot call its L.

## The reference set

12 of 13 CTRV genomes. The intron was located **empirically in each genome**,
not transcribed from the paper: search for a 76-nt `GTAAGT`…`AG` intron whose
excision extends the L ORF past 2000 aa. That independently reproduced the
published coordinates — NC_025384 gives **8648-8723**, exactly Kuwata's
figure — and the junction windows are identical in all 12
(`CTTTGGTAAGT` / `TTTAGGTTTAC`).

Two things that are easy to get wrong:

**The coordinates in the header are local to the reference sequence**, and
the reference must be **the L gene region from its start codon**, not the
whole genome. The program builds the 5' piece with
`substr($sseq, 0, $splice_left_end+1)` — from the start of the *aligned*
region — and then refuses any 5' piece that translates with a stop. With an
11 kb genome as the reference that piece is everything upstream of L, so it
always has stops and the splice is silently skipped. With the gene region it
is clean, and all 12 references then share the same local coordinates:
`SD:3962-3972;3966 SA:4038-4048;4043`.

**PQ305977 is deposited in genome-sense**, the reverse complement of the
others. It is stored reverse-complemented here so its coordinates are
consistent with the rest.

**OL700055 is excluded, and stays excluded.** It is a metagenome-assembled
genome. It carries the canonical `GTAAGT` donor at 8638 but a single base
destroys the acceptor: where all twelve others read `TTTAG|GTTTAC`, this reads
`TTTTG|GTTTAC` — A to T, so there is no AG. The nearest AG that would rescue
the ORF sits 106 nt downstream and gives 2113 aa instead of 2123. In a MAG
that is far more likely an assembly error than biology, and adding a
reference with a non-canonical acceptor would encode that error into the set
where it could mis-splice other genomes. It is the one CTRV genome that still
gets two L calls.

## The fragment-retirement question, answered

The concern recorded before this was built: does the splice call *replace* the
truncated L, or duplicate it? It duplicated. The ordinary L profile hits both
exons as two fragments — 1333 aa and 801 aa — and neither is a protein, so a
CTRV genome ended up with three features all annotated
`RNA-dependent RNA polymerase`, which `viral_genome_quality.pl` counts against
`copy_num` and reports as "too many HSPs".

Adding was nonetheless the right default and is why no removal existed: in
influenza every spliced product has a **different** annotation from its
unspliced parent — M1/M2, NS1/NS2, PA/PA-X — because both are real proteins.
`Merhavirus L_SPLICED` is the **only** splice feature in the whole installed
set that shares an annotation with a sibling, so the retirement added to
`get_splice_variant_features.pl` provably cannot fire anywhere else.

A fragment is recognised by sharing an exon boundary with the join **and**
carrying the same function. Containment alone is not enough: the truncated L
runs past its exon into the intron to reach its stop.

## Result

    CTRV NC_025384    N 472 · P 264 · M 173 · G 503 · G_mat 484
                      L 2123 aa   4682..8647 + 8724..11129

    Merida virus      L 2136 aa   single exon, no splice, nothing retired

**12 of 13** CTRV genomes give exactly one L call with two exons. The
thirteenth is OL700055, described above: a MAG whose acceptor site is
disrupted by one base, which still gets the two fragment calls. Non-CTRV
Merhavirus is untouched — Merida virus keeps its single-exon 2136 aa L and
nothing is retired.
