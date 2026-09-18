# LowVan module kit

**Name the proteins in a viral genome, consistently, across a whole taxonomic
family.**

Public sequence databases hold tens of thousands of viral genomes whose protein
names were written by thousands of different submitters, so the same protein
appears under dozens of spellings and many carry no useful name at all. A LowVan
*module* replaces that with one controlled set: a collection of PSSM profiles
plus a JSON declaring which proteins a taxon has, how long they are, and how
confident a match must be.

This repository is everything needed to build such a module for a new viral
taxon, run it, and find out whether it actually works — plus the agent skill that
walks through doing so.

Building one is mostly judgement: which annotation strings mean the same protein,
which clusters are real, which features the data cannot support. **The skill
exists because those judgements are easy to get wrong in ways nothing downstream
detects.**


## Where the work lives

**This repository is the source of truth.** `build/` and `annotate/` are the
pipeline and the annotator; `skill/` is the method; `modules/` holds the built
modules; `reports/` holds the published pages.

`CEPI-dxkb/Viral_Annotation` is a *downstream consumer*, not an input. It has
not received the annotator work done here — at the time of writing it has no
`internal_stop` support at all, so a module relying on a readthrough is
silently cropped against it. Do not build against a fresh clone of it and do
not treat `patches/` as the current state: patches are diffs against pristine
upstream and drift behind `build/` and `annotate/` as work accumulates.

To set up a runtime directory for the annotator, copy `annotate/*` and
`build/*` over it first. See `modules/README.md`.

## Does it work?

The Rhabdoviridae module built with this kit, measured against **every**
Rhabdoviridae genome BV-BRC holds between 1 and 50 kb — 25,237 records collapsed
to 1,356 distinct genomes, one per 95%-identity cluster:

| | |
|---|---|
| genomes assigned a module and annotated | **88%** |
| genomes carrying all five core proteins | **69%** |
| genomes annotated acceptably | **70%** |
| distinct annotation strings produced | **20**, against 255 for the same genomes in BV-BRC |

`reports/` holds four standalone HTML pages with the full evidence. GitHub will
not render them in place — download the file, or use a raw-HTML viewer.

| Report | What it answers |
|---|---|
| `rhabdoviridae-coverage-audit.html` | the whole-taxon validation: what the module covers, how it was built, and the rules that governed it. **Start here.** |
| `bvbrc-string-collapse.html` | how many source annotation strings each controlled annotation absorbs |
| `rhabdoviridae-pssm-registry.html` | every declared feature, its profiles, and what fraction of its own collection each recovers |
| `lowvan-module-handbook.html` | the build method in long form |

The audit is also worth reading for what a whole-taxon run *finds*: four bugs in
the annotator, a module worth dropping, five reference genomes worth removing,
and a mis-declared cleavage site that had been misplacing a protein's C-terminus
in 40% of genomes. None of it was visible from the eight-genome reference panel,
which the module passed at 95% throughout.

## What is in here

~1.3 MB. The PSSM and alignment data a *finished* module ships is not included;
see **Getting the reference data**.

```
skill/                 the lowvan-module skill: 15 steps, 10 reference documents,
                       29 scripts
build/                 clustering and profile construction            (patched)
annotate/              the annotator and genome-quality scoring       (patched)
evaluate/              reference-panel scoring, whole-taxon coverage, download
example-rhabdoviridae/ the two taxon-specific programs, as a worked example
example-togaviridae/   a second example: a -1 frameshift product no PSSM can
                       call, a readthrough codon inside a protein, and two
                       precursors recovered from genome intervals
example-matonaviridae/ a third example: what to do when the SOURCE ANNOTATION is
                       the problem -- 126 wrong-frame records that built two
                       convincing junk profiles, a quarter of the vocabulary
                       naming no protein, and a binning proved from sequence
patches/               7 diffs against upstream Viral_Annotation
vendor/bv-brc/         the BV-BRC/SEED modules the whole toolchain needs, so a
                       clone runs with nothing else installed  (SEED licence)
check-environment.sh   verifies binaries, CPAN modules, and that it all compiles
reports/               the artifacts these produced
CHANGELOG.md           what differs from upstream and why, with the measurements
```

## What you adapt, and what runs unchanged

Two programs are taxon-specific: `build_collections.py` (which annotation strings
mean which protein) and the module JSON generator (what the taxon has, how long,
how confident a match must be). Both are in `example-rhabdoviridae/` as working
code with the parts to replace marked, and its README lists them table by table.
`example-togaviridae/` is the second case: read it when a protein is not
collinear with the genome, or comes out short by a constant number of residues.
`example-matonaviridae/` is the third: read it when the export you were handed is
mislabelled, fragmentary or frame-shifted, and when you need to prove a
mature-peptide binning without trusting the annotation you are replacing.

Everything else runs unchanged — clustering, profile construction, curation QC,
signal-peptide derivation, install, reference-panel scoring and whole-taxon
evaluation.

That division is real rather than a shortcoming. Deciding that `outer coat
protein` is the glycoprotein, that `M1` is the phosphoprotein, or that `P6` must
never be matched by regex at all is the work, and it is different for every
taxon. The kit carries the machinery and the method; the judgement is yours, with
roughly a thousand lines of comments in the example recording how each call went
here and what measurement settled it.

**The organising rule: anything algorithmic is a script, anything requiring
judgement is written instructions.** Every analysis in the coverage report is a
program you can re-run on another taxon. What stayed prose is what does not
transfer — which genes flank a gap, how many reference contigs are worth
maintaining, whether a taxon has enough data to model at all.

## Requirements

Run this first; it checks everything below and tells you what is missing:

```bash
./check-environment.sh
```

**External binaries** — conda is easiest:

```bash
conda install -c bioconda blast mmseqs2 mafft     # plus python 3.9+
```

**Perl modules from CPAN** — these are the only ones not bundled:

```bash
cpanm File::Slurp Data::UUID JSON::XS Getopt::Long::Descriptive \
      Class::Accessor IPC::Run Math::Round
```

`Data::UUID` (or `UUID`) is what `GenomeTypeObject` mints feature ids with. It
loads inside an `eval`, so without it a run dies partway through with
`No UUID generator found` rather than at startup — `make_gto.py` and
`check-environment.sh` both check for it up front.

**Everything else is bundled.** `vendor/bv-brc/` carries the BV-BRC and SEED
modules — `GenomeTypeObject`, `IDclient`, `P3DataAPI`, `BlastInterface`,
`gjoseqlib` and their closures. **You do not need a BV-BRC installation or a
`seed_gjo` checkout.** All eight Perl entry points compile against
`vendor/bv-brc/lib` and `build/` alone, which `check-environment.sh` verifies.

```bash
export PERL5LIB="$PWD/vendor/bv-brc/lib:$PWD/build:$PERL5LIB"
export PATH="$PWD/annotate:$PWD/vendor/bv-brc/bin:$PATH"
```

## Getting the reference data

The profiles and alignments for already-built modules live in the upstream
repository, <https://github.com/CEPI-dxkb/Viral_Annotation>. Clone it and point
`LOWVAN_DATA_DIR` at it.

You do **not** need it to build a module for a new taxon. You need it to run the
annotator against taxa that already have modules.

## Quick start

```bash
export LOWVAN_DATA_DIR=/path/to/Viral_Annotation
export PERL5LIB="$PWD/vendor/bv-brc/lib:$PWD/build:$LOWVAN_DATA_DIR:$PERL5LIB"
export PATH="$PWD/annotate:$PWD/vendor/bv-brc/bin:$PATH"

# annotate one genome
perl annotate/annotate_by_viral_pssm.pl -i genome.fna -p out

# score a module against reference genomes
python3 evaluate/evaluate_module.py --repo $LOWVAN_DATA_DIR --acc NC_001542

# score it against everything the taxon has in BV-BRC
perl evaluate/New-annotate-viral-taxon.pl -i <Taxon> -download-only -d Contigs \
     -metaout <Taxon>.metadata -min 1000 -max 50000
python3 evaluate/evaluate_coverage.py --contigs Contigs --repo $LOWVAN_DATA_DIR \
        --json $LOWVAN_DATA_DIR/Viral_PSSM.json
```

### The GTO path, for quality scoring

Coverage above tells you which proteins were called. It cannot tell you whether a
genome is *well* annotated, because the flat feature table cannot represent
special features. Quality scoring runs on GTOs:

```bash
# rejoin segmented genomes first, or every segment counts as a broken genome
python3 skill/scripts/merge_segments.py --json module.json --fasta-dir Contigs \
        --metadata <Taxon>.metadata --out Contigs-merged

# fasta -> GTO. ids are ISSUED by the ID server; never hand-write this JSON
python3 skill/scripts/make_gto.py --fasta-dir Contigs-merged \
        --metadata <Taxon>.metadata --out gto --jobs 24

# annotate and quality-score
python3 skill/scripts/run_gto_eval.py --gto-dir gto --repo $LOWVAN_DATA_DIR \
        --out gto_out --report quality.tsv \
        --perl /Applications/BV-BRC.app/runtime/bin/perl \
        --perl5lib /Applications/BV-BRC.app/deployment/lib
```

Three environment problems are near-universal, and all three are flags rather
than edits:

| symptom | fix |
|---|---|
| `Can't locate GenomeTypeObject.pm` / `IDclient.pm` | `--perl5lib <devkit>/deployment/lib` |
| `No UUID generator found` | `GenomeTypeObject` needs `Data::UUID` or `UUID`. The BV-BRC runtime perl has one; a conda perl usually does not — `--perl <devkit>/runtime/bin/perl` |
| `rast-create-genome: command not found` | the macOS bundle ships `plbin/rast-create-genome.pl` without its `bin` wrapper; `make_gto.py` writes it from the bundle's own pattern |

Only the GTO wrapper and the quality script run under the dev-kit perl. The inner
`annotate_by_viral_pssm.pl` resolves through `PATH`, so it keeps the interpreter
its own dependencies are installed under.

### Analysis scripts

These live in `skill/scripts/` — one copy, so they cannot drift from the skill
that documents them:

| Script | Purpose |
|---|---|
| `merge_segments.py` | rejoin multi-segment genomes, gated on the declared count |
| `make_gto.py` | contig FASTAs to GTOs via `rast-create-genome` |
| `run_gto_eval.py` | annotate and quality-score GTOs; good-vs-poor breakdown |
| `minimal_refs.py` | fewest reference contigs that close a routing gap |
| `analyze_feature_gap.py` | why is one protein called less often than its neighbours? |
| `check_cleaved_ends.py` | are `cleave=` termini crisp, or is the alignment drifting? |
| `annotation_rarefaction.py` | vocabulary growth, controlled versus free text |
| `build_leftover_pssms.py` | second-pass profiles from the sequences no cluster took |
| `check_dump.py` | verify the BV-BRC protein dump is complete and line-aligned |

## Licensing, and what is deliberately not here

This repository is MIT licensed. The annotator and clustering code in
`annotate/` and `build/` is redistributed from
[CEPI-dxkb/Viral_Annotation](https://github.com/CEPI-dxkb/Viral_Annotation),
which is also MIT and the same copyright holder, so it travels cleanly.

**SignalP is not included and cannot be.** SignalP 6.0 is academic-licensed and
distributed only through a form you submit yourself. `skill/scripts/apply_signalp.py`
is a wrapper that runs SignalP and consumes its output; it contains no SignalP
code, binaries, or model weights. Get it from
<https://services.healthtech.dtu.dk/services/SignalP-6.0/> if you want it.

**You do not need it.** Step 6 of the skill carries a decision table for building
without SignalP. Where a mature N-terminus is documented in the literature or
recorded as a `mat_peptide` in RefSeq, `apply_signalp.py --motif` cuts on that
motif with no predictor involved — and that is the better answer even when
SignalP is available, because a motif re-locates itself after a rebuild while a
column number does not. Where no such site exists, the signal-peptide and mature
products are dropped and the module ships without them; the parent protein is
still called correctly. Nothing else in the workflow depends on SignalP.

**The BV-BRC modules are included**, in `vendor/bv-brc/`: the dependency closure
of `GenomeTypeObject.pm`, `IDclient.pm` and `rast-create-genome`, 11 modules and
one script, so the GTO path runs from a clone with nothing else installed. They
are SEED Toolkit files under the **SEED Toolkit Public License, which permits
redistribution**; copyright is University of Chicago and the Fellowship for
Interpretation of Genomes, 2003-2013. That licence rather than MIT governs those
files, so their notices stay with them. `vendor/README.md` and
`vendor/bv-brc/MANIFEST.md` carry the attribution and exact provenance.

Everything else in this repository is either original to it or MIT from
upstream. The controlled vocabulary ships as `skill/assets/annotation-vocabulary.tsv`
rather than as the manuscript spreadsheet, so there is no publisher-copyright
question; `query_PATRIC_bob.pl` and `P3DataAPI` are both fully open.

## Using the skill

Copy `skill/` to `~/.claude/skills/lowvan-module/` and invoke it, or read
`skill/SKILL.md` directly — it is written to be followed by hand. Fourteen steps
run from taking delivery of a BV-BRC export through to verifying that every
cleaved protein end lands where it should.

Three documents carry most of the reasoning:

- `skill/references/annotation-triage.md` — turning source annotation strings
  into feature collections, and what must never be grouped by regex
- `skill/references/curation.md` — what to fix in an alignment and what to leave
- `skill/references/module-partitioning.md` — deciding where one module ends and
  the next begins, which is the decision everything downstream depends on

## Relationship to upstream

This kit carries **patched** copies of several upstream programs. `patches/`
holds the diffs and `CHANGELOG.md` explains each one with the measurement that
prompted it. The patches are not upstreamed; if you are working from a fresh
clone of [CEPI-dxkb/Viral_Annotation](https://github.com/CEPI-dxkb/Viral_Annotation),
apply them or use the copies here.
