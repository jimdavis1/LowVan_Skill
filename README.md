# LowVan module kit

Everything needed to build a LowVan annotation module for a viral taxon, run it,
and find out whether it actually works — plus the skill that walks an agent
through doing so.

Roughly 900 KB. The PSSM and alignment data that a finished module ships is not
here; see **Getting the reference data** below.

```
skill/                 the lowvan-module skill: workflow, references, scripts
build/                 clustering and profile construction  (patched)
annotate/              the annotator and genome-quality scoring  (patched)
evaluate/              reference-panel scoring, whole-taxon coverage, download
example-rhabdoviridae/ the two taxon-specific programs, as a worked example
patches/               diffs against upstream Viral_Annotation
reports/               the four artifacts this produced, as standalone HTML
CHANGELOG.md           what differs from upstream and why, with the measurements
```

## What you write, and what you get

Two programs are taxon-specific and you will adapt them:
`build_collections.py` (which annotation strings mean which protein) and the
module JSON generator (what the taxon has, how long, how confident a match must
be). Both are in `example-rhabdoviridae/` as working code with the parts to
replace marked, and its README lists them table by table.

Everything else runs unchanged: clustering, profile construction, curation QC,
signal-peptide derivation, install, and both evaluation paths.

That division is real rather than a shortcoming. Deciding that `outer coat
protein` is the glycoprotein, that `M1` is the phosphoprotein, or that `P6`
must never be matched by regex at all is the work — and it is different for
every taxon. The kit carries the machinery and the method; the judgement is
yours to make, with roughly a thousand lines of comments in the example
recording how those calls went here and what measurement settled each one.

## What this is for

A LowVan module is a set of PSSM profiles plus a JSON that tells the annotator
which proteins a taxon has, how long they are, and how confident a match must
be. Building one is mostly judgement: which annotation strings mean the same
protein, which clusters are real, which features the data cannot support. The
skill exists because those judgements are easy to get wrong in ways nothing
downstream detects.

## Requirements

```
perl        JSON::XS File::Slurp Getopt::Long Cwd
            gjoseqlib.pm, BlastInterface.pm   github.com/TheSEED/seed_gjo
binaries    blastn tblastn psiblast makeblastdb mmseqs mafft
python      3.9+
```

For anything that talks to BV-BRC — the protein dump that feeds step 1
(`query_PATRIC_bob.pl`, which uses `P3DataAPI`), the genome download in
`New-annotate-viral-taxon.pl`, the GTO path, and `viral_genome_quality.pl` — you
need the BV-BRC dev kit:

```
GenomeTypeObject.pm, IDclient.pm, P3DataAPI.pm, rast-create-genome
  from BV-BRC-CLI-<ver>.tgz  github.com/BV-BRC/BV-BRC-CLI/releases
  or, on macOS, the BV-BRC.app desktop bundle, which carries the same tree
  under /Applications/BV-BRC.app/{deployment,runtime}
  put dev_container/modules/*/lib on PERL5LIB
  put the CLI scripts and this kit's annotate/ on PATH

plus  Class::Accessor  Getopt::Long::Descriptive  Params::Validate  URI  LWP
```

`Params::Validate` is XS. Without a compiler, install it from conda-forge
(`perl-params-validate`) rather than CPAN.

## Getting the reference data

The profiles and alignments for already-built modules live in the upstream
repository, <https://github.com/CEPI-dxkb/Viral_Annotation>. Clone it and point
`LOWVAN_DATA_DIR` at it.

You do **not** need it to build a module for a new taxon. You need it to run the
annotator against taxa that already have modules.

## Quick start

```bash
export LOWVAN_DATA_DIR=/path/to/Viral_Annotation
export PERL5LIB=$LOWVAN_DATA_DIR:/path/to/seed_gjo:$PERL5LIB
export PATH=$PWD/annotate:$PATH

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

## What it produced

`reports/` holds four self-contained HTML pages. Open them in a browser; they
need no server and no network.

| Report | What it answers |
|---|---|
| `rhabdoviridae-coverage-audit.html` | the whole-taxon validation: what the module covers across every Rhabdoviridae genome in BV-BRC, how it was built, and the rules that governed it. **Start here.** |
| `bvbrc-string-collapse.html` | how many distinct BV-BRC annotation strings each controlled annotation absorbs |
| `rhabdoviridae-pssm-registry.html` | every declared feature, its profiles, and what fraction of its own collection each recovers |
| `lowvan-module-handbook.html` | the build method in long form |

The headline numbers from the coverage audit, for a module of 955 profiles over
1,356 distinct genomes drawn from 25,237 BV-BRC records:

```
annotated, in scope                  88%
all five core proteins               69%     (was 62% before the last rebuild)
annotated acceptably                 70%     (excluding sequence-quality flags)
annotation vocabulary                20 strings vs 255 for the same genomes
```

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

**The BV-BRC dev kit is not included either.** `GenomeTypeObject.pm`,
`IDclient.pm` and `rast-create-genome` come from
[BV-BRC-CLI](https://github.com/BV-BRC/BV-BRC-CLI/releases). They are required
for the GTO path only.

Everything else in this repository is either original to it or MIT from
upstream. The controlled vocabulary ships as `skill/assets/annotation-vocabulary.tsv`
rather than as the manuscript spreadsheet, so there is no publisher-copyright
question; `query_PATRIC_bob.pl` and `P3DataAPI` are both fully open.

## Using the skill

Copy `skill/` to `~/.claude/skills/lowvan-module/` and invoke it, or read
`skill/SKILL.md` directly — it is written to be followed by hand. The thirteen
steps run from taking delivery of a BV-BRC export through to chasing any feature
that loses 20% or more of its genomes.

Two documents carry most of the reasoning:

- `skill/references/curation.md` — what to fix in an alignment and what to leave
- `skill/references/annotation-triage.md` — turning source annotation strings
  into feature collections, and what must never be grouped by regex

## Relationship to upstream

This kit carries **patched** copies of several upstream programs. `patches/`
holds the diffs and `CHANGELOG.md` explains each one with the measurement that
prompted it. The patches are not upstreamed; if you are working from a fresh
clone of `Viral_Annotation`, apply them or use the copies here.
