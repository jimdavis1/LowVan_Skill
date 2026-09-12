# Installing and running the module

The point of this step is that a module which has never been run is not
finished. Every defect described in the QC references was found by running the
annotator, not by inspecting the alignments.

## Get the code

```bash
git clone https://github.com/CEPI-dxkb/Viral_Annotation.git
export LOWVAN_DATA_DIR=$PWD/Viral_Annotation
```

`LOWVAN_DATA_DIR` overrides the hard-coded default paths, so you do not need to
edit the scripts.

## Dependencies

The main annotator (`annotate_by_viral_pssm.pl`) needs only:

```
perl  JSON::XS  File::Slurp  Getopt::Long  Cwd  gjoseqlib
blastn  tblastn        (2.16.0+ — the JSON output format is version-specific)
```

Building PSSMs additionally needs `psiblast makeblastdb mmseqs mafft`, and
`BlastInterface.pm`.

```bash
for m in JSON::XS File::Slurp Getopt::Long Cwd gjoseqlib BlastInterface; do
  printf "%-24s " $m; perl -M$m -e 'print "ok\n"' 2>/dev/null || echo MISSING
done
for b in blastn tblastn psiblast makeblastdb mmseqs mafft; do
  printf "%-14s " $b; command -v $b >/dev/null && echo ok || echo MISSING
done
```

`gjoseqlib.pm` and `BlastInterface.pm` come from
<https://github.com/TheSEED/seed_gjo/> — put them on `PERL5LIB`. Install missing
CPAN modules with `cpanm --local-lib=~/perl5 File::Slurp` and add
`~/perl5/lib/perl5` to `PERL5LIB`.

The GTO pipeline (`annotate_by_viral_pssm-GTO.pl`, transcript editing, splice
variants, quality) additionally needs `GenomeTypeObject.pm` and `P3DataAPI` from
BV-BRC. You do not need them to build and test a module from FASTA.

## Install

```bash
python3 scripts/install_module.py --workdir . --repo $LOWVAN_DATA_DIR --check
python3 scripts/install_module.py --workdir . --repo $LOWVAN_DATA_DIR
```

Copies PSSMs to `Viral-PSSMs/<Module>.pssms/<FEAT>/`, alignments to
`PSSM-Alignments/<Module>/<FEAT>/`, rep contigs to `Viral-Rep-Contigs/`, and
merges the JSON block into `Viral_PSSM.json` (backing it up first). Idempotent —
re-run after any rebuild.

The master is written in the canonical JSON::XS format. If the checkout's master
was not already canonical, the first install reformats every line and says so.
That is one-time: afterwards a diff shows only real changes, and Perl tools stop
reordering the file. Because the master is shared, mention that reformat in the
commit message rather than letting a reviewer discover a 7,000-line diff.

`--check` refuses to install when the JSON and the PSSMs disagree. Fix that
before proceeding; both directions of mismatch are silent at runtime.

## Representative contigs

Routing is a **BLASTn** of the incoming contig against `Viral-Rep-Contigs/`.
This is the step people under-provision, because nucleotide similarity falls off
far faster than protein similarity.

```bash
python3 scripts/check_rep_contigs.py --repdir Rep-Contigs --acc-file accs.txt
```

Rice yellow stunt virus returned **zero BLASTn hits** against five
Betarhabdovirinae rep contigs, one per genus. It was rejected before any PSSM
ran, and the fix was another rep contig, not a better profile. Budget one per
genus as a starting point and add more wherever a test genome fails to route.

Rep contig format is FASTA with the header `>accn|<ACCESSION>   <description>`.
Fetch from NCBI and never guess an accession — look each up and record what
confirmed it.

## Run it

```bash
export LOWVAN_DATA_DIR=/path/to/Viral_Annotation
perl $LOWVAN_DATA_DIR/annotate_by_viral_pssm.pl -i contigs.fasta -p out -threads 4
```

Produces `out.faa`, `out.ffn`, `out.feature.tbl`. Useful debugging flags:

```
-no -tbl     feature table to stdout, no files
-no -aa      proteins to stdout
-tmp         keep the temp directory
-mcb 150     minimum contig bitscore to accept the genome at all
```

Running with `-no -tbl` prints, for every feature, the score of every PSSM
tried and which one was chosen. That output is how you diagnose a wrong call —
it shows you the competing profile that won.

## Validate one genome

```bash
python3 scripts/validate_calls.py --repo $LOWVAN_DATA_DIR --fasta genome.fna
```

No reference annotation needed, which is the point: this is the check you can
run on a genome nobody has annotated. It reads the proteins the annotator
actually emitted and asks whether they are sound — M start, no internal stop,
length inside the JSON's declared range, coordinates on the contig, no two
features on the same span.

`DERIVED_GAP` is the subtle one. A mature peptide is derived from a parent CDS,
and the derived profiles often cover fewer clusters than the parent does. On
Maraba virus the G CDS was called at 512 aa and **every one of the ten G_MAT
PSSMs scored 0** — not a threshold miss, no tblastn hit at all. The genome got a
glycoprotein and no mature peptide, silently. Nothing else in the toolchain
notices that, because the CDS call is perfectly correct.

## Score it

```bash
python3 scripts/evaluate_module.py --repo $LOWVAN_DATA_DIR --acc-file accs.txt
```

Fetches each accession's sequence and its GenBank CDS coordinates, runs the
annotator, and compares. GenBank CDS come from submitters, so they are a
yardstick and not truth — the recall percentage matters less than two specific
outcomes:

**Duplicates are always a defect.** Two features called on identical
coordinates means a cluster is binned under the wrong feature and its profile
out-scores the correct one. Go to `qc_cross_feature.py`.

**Misses are either unbuilt features or a cutoff set too high.** Check which
before touching a cutoff.

Worked result from the Rhabdoviridae module across four taxa and eight held-out
genomes: 36 of 41 reference CDS matched, **0 duplicates**, 0 unmatched calls.
The five misses were all accessory ORFs annotated `hypothetical protein` in
plant rhabdoviruses — features the module deliberately does not declare.

## The loop that found the real bug

Worth internalising, because it is the reason this step exists:

1. ran the annotator on a held-out rabies genome
2. it emitted `Nucleocapsid protein` on the **glycoprotein** coordinates, and
   never called N at all
3. `-no -tbl` showed `Alpharhabdovirinae.N.10.pssm` scoring 1049 bits at the G
   locus, beating the correct N profile at 890
4. that PSSM's query was 524 aa — G-sized, not N-sized
5. all nine members of N cluster 10 were 98–99% identical to glycoproteins, and
   all nine were annotated `nucleoprotein` in BV-BRC
6. retired the cluster; re-ran; all five rabies genes called correctly

No amount of looking at the N alignment would have revealed this. It looked
fine, because it *was* a fine alignment — of the wrong protein.
