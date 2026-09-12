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

For the GTO path — `annotate_by_viral_pssm-GTO.pl`, `viral_genome_quality.pl`
and the whole-taxon evaluation — you also need the BV-BRC dev kit:

```
GenomeTypeObject.pm, IDclient.pm, rast-create-genome
  from BV-BRC-CLI-<ver>.tgz  github.com/BV-BRC/BV-BRC-CLI/releases
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
