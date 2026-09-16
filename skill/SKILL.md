---
name: lowvan-module
description: Build a complete LowVan viral annotation module for a taxon - triage BV-BRC annotation strings into feature collections, cluster and align them, curate the alignments, build PSSMs, write the Viral_PSSM.json block, then install into a Viral_Annotation checkout and run the annotator to test it. Use when asked to build, extend, rebuild, or QC a LowVan/PSSM annotation module for any viral family, subfamily, or genus.
---

# Building a LowVan module

A LowVan module makes one viral taxon annotatable. You start from a BV-BRC
protein dump and finish with four things the annotator loads at runtime:

| Artifact | Path in the repo | What it does |
|---|---|---|
| JSON block | `Viral_PSSM.json` → `"<Module>"` | per-feature cutoffs, lengths, annotation strings |
| PSSMs | `Viral-PSSMs/<Module>.pssms/<FEAT>/` | the profiles that call features |
| Rep contigs | `Viral-Rep-Contigs/<Module>.<n>.dna` | BLASTn subjects that route a genome to the module |
| Alignments | `PSSM-Alignments/<Module>/<FEAT>/` | provenance; not read at runtime |

**It is not a discovery tool.** It calls proteins that a curator decided exist.
Every judgement is yours; the scripts only make judgements cheap to check.

## The one rule that matters most

**A module is defined by genome organisation, not phylogeny.** One module per
group of taxa that share a gene layout. Genera whose genomes carry the same five
core genes collapse into one module even across subfamilies; a genus with an
extra segment or a different accessory set gets its own. Getting this wrong is
expensive later, so settle it before building collections
(`references/module-partitioning.md`).

## Workflow

Work through these in order. Do not skip 7 or 9 — they catch defects the earlier
steps cannot see.

### 0. Prerequisites

```bash
git clone https://github.com/CEPI-dxkb/Viral_Annotation.git
export LOWVAN_DATA_DIR=$PWD/Viral_Annotation
```

Needs `perl` with `JSON::XS File::Slurp Getopt::Long Cwd`, plus `gjoseqlib.pm`
and `BlastInterface.pm` from <https://github.com/TheSEED/seed_gjo/> on
`PERL5LIB`; and `blastn tblastn psiblast makeblastdb mmseqs mafft` on `PATH`.
Full check and install notes: `references/install-and-test.md`.

Optional, and only if the taxon has a cleaved precursor: **SignalP 6.0**, for
step 6. It is academic-licensed and cannot be installed unattended — the user
downloads it themselves through a form. Raise it early if the taxon looks like it
will need one, because the download is not instant.

**Never block the build on it.** Step 6 has a decision table for the case where
the user does not have SignalP and cannot easily get it: a documented cleavage
site beats a prediction anyway, and where neither exists the cleaved products are
dropped and the module ships without them. Everything else in the workflow runs
unchanged.

### 1. Take delivery of the inputs

Five tab-separated files from BV-BRC, named for the taxon. The three `uniq.*`
files are **line-index aligned** — that alignment is the join key, so never sort
them independently.

```
<T>.id_md5        feature_id \t md5                     every feature
<T>.id_name_gs    genome_id \t genome_name \t family \t genus
<T>.uniq.md5      md5                                  one per unique sequence
<T>.uniq.id_ann   representative_feature_id \t annotation
<T>.uniq.seq      md5 \t sequence
```

```bash
python3 scripts/check_dump.py --workdir .
```

**Run that before anything else.** It checks that every protein md5 named in
`id_md5` actually has a sequence, and that the three `uniq.*` files are
line-aligned. The step that builds `uniq.md5` selects protein-coding features
with a positive filter on the feature id, and omitting one shape — `peg` — cost
this project 45% of its protein diversity with no error anywhere. The commands
that produce the dump, and the safer form of that filter, are in
`references/inputs.md`.

### 2. Triage the annotation strings

Histogram every annotation string by genus, then write rules mapping strings to
feature keys. This is the step that determines everything downstream, and the
step where mistakes are least visible.

Read `references/annotation-triage.md` before writing rules — it covers the
**U-number trap** (positional labels like `U1`/`ORF3`/`VP2` are not homology
groups), genus-dependent synonyms, and why every rule needs a genus column.

Output: `collections/<Module>/<FEAT>.fasta`, plus `synonyms.tsv`,
`UNASSIGNED_TRACKING.tsv` and `NOT_MODELLED.tsv` so nothing vanishes silently.

### 3. Build the alignments and PSSMs

```bash
cat collections/<Module>/<FEAT>.fasta \
  | fasta-cluster-pssm-2.pl -a "<Annotation string>" -p <Module>.<FEAT> \
      -m 5 -f 0.33 -n 0.75 -c 0 -fd 0.20 -e 3 -efo 0.15 -mi 0.8 -mc 0.8 \
      -p_nterm 65 -n_nterm 3 -nterm_eval_length 10 -nterm-mmseq-id 0.7
```

Those are the documented defaults. `-m 5` is a floor to aim for, not a law: drop
to `-m 2` rather than lose a taxon, and record every departure in a
`BUILD_PARAMS` file next to the feature. When coverage is poor, `-mi` (clustering
identity) is usually the binding constraint, not `-m`.

The upstream pipeline has seven bugs worth knowing about, with patches in
`assets/patches/`. Read `references/pipeline.md` before the first run.

#### Then build profiles for the leftovers

`Leftover_Seqs.aa` holds everything that reached no cluster. It is tempting to
read that as noise. It is the opposite: the leftovers are the **divergent
members**, which is precisely the population a profile set most needs to cover,
and dropping them silently caps what the module can ever call.

```bash
python3 scripts/build_leftover_pssms.py --workdir . --module <M> --key <FEAT>
python3 scripts/build_leftover_pssms.py ... --write
```

It re-clusters the leftovers on their own at the **same** `-mi` — the identity
floor is never lowered — so they get their own chance instead of competing with
the dense clusters that won the first pass. New profiles are numbered `lo1`,
`lo2` … so they cannot collide with first-pass cluster ids.

**The member floor drops to 2 here, and only here.** That is not the same
concession as lowering `-mi`. The leftover pool is *defined* by having failed
the feature's `-m`, so reapplying that `-m` guarantees nothing comes back:
adding 21 sequences to Alpharhabdovirinae P moved exactly one group over the
`-m 5` line in the main pass. Identity still does the quality work — two
sequences that cluster at 0.6 are genuinely related — while the member count
is only asking whether a profile can be built at all.

The evidence this matters is not theoretical. Alpharhabdovirinae P held 2727
sequences, 2231 inside one of 24 profiles and 282 left over, and the module
called a phosphoprotein in only 68% of genomes against 99% for N. Taking the
115 genomes that had a called N and M but no P, and ORF-calling the interval
between them:

```
N->M interval present           115/115   median 1018 nt
contains an ORF                 115/115   median  296 aa
ORF matches the P collection    109/115   median bitscore 513

closer to a LEFTOVER P            85      (no profile covers it)
closer to a CLUSTERED P            0      (a profile exists)
```

Zero. Every uncalled gene was full length, unambiguous, recognisably the right
protein — and none of them resembled anything a profile was built from. The
uncovered population and the uncalled genes are the same set, so the leftovers
are a direct readout of what the module will miss in the wild.

**Check the leftovers before concluding a feature is finished**, and note that
neither self-recall nor a reference panel would have shown this: the profiles
find their own training data perfectly, and the eight-genome panel scores 39/41.

### 4. Curate the alignments

The pipeline produces drafts. Curation is the work, and it is manual — open
`corrected_alis/*.fa` in Jalview and fix them. Non-negotiable checks, from the
LowVan curator's checklist:

- every profile starts at the right Met, or you know why it does not
- no PSSM is a truncated version of a longer one
- the sequences actually encode the protein the key claims
- ends are as similar as possible; a jagged N-terminus means split the cluster

`references/curation.md` has the mechanics: N-terminal splitting, the
truncation QC, identical-sequence collapse, mature-peptide and signal-peptide
derivation, and the Met rule (**flag, never delete**). A precursor often yields
several mature peptides — GPC to signal peptide + Gn + Gc, HA0 to HA1 + HA2 —
and each product's `upstream_ext` / `downstream_ext` follow from its own two
termini; see `references/json-schema.md`.

### 5. Rebuild PSSMs after curating

```bash
python3 scripts/rebuild_pssms.py            # audit
python3 scripts/rebuild_pssms.py --write    # rebuild the stale ones
```

Rebuilds every PSSM from its current alignment and compares bytes, so it finds
alignments edited since their profile was built. Do not trust file timestamps —
cloud-synced directories flatten them. Do not trust the PSSM's embedded query
either: deleting a non-master sequence leaves the query untouched and still
changes every score.

### 6. Derive the cleaved products

Signal peptides and mature chains have no collection of their own — nothing in
BV-BRC is annotated as "mature G". They are a coordinate operation on the parent
CDS you just curated, so they come after the parent is final and must be redone
whenever it changes.

**Re-stage the queries every time the parent is rebuilt.** The ids carry cluster
numbers, and re-clustering renumbers them. A prediction keyed to a cluster that
has since been split would be applied to the wrong sequences, silently.

```bash
python3 scripts/apply_signalp.py --workdir . --parent G \
        --stage-queries Alignments/G_SIGNALP_QUERIES.faa
```

#### First decide whether you need SignalP at all

**SignalP is optional, and a module without cleaved products is still a valid
module.** Work down this list and stop at the first one that applies:

| situation | what to do |
|---|---|
| the mature N-terminus is documented for a reference protein | use `--motif`. No predictor needed, and this is the *better* answer even when SignalP is available — see below |
| RefSeq carries a `mat_peptide` on a reference genome of this taxon | read the cut site off it, turn it into a motif, use `--motif` |
| neither, and SignalP is unavailable | **drop the `_SP` / `_MAT` features** from the JSON and say so in the build notes. The parent CDS is still called correctly; you lose two derived features, not the module |
| neither, and SignalP is available | run it, then apply it with the checks below |

If the user does not have SignalP, do not stall the build waiting for it. Ask
whether a documented cleavage site exists for any protein in the taxon; if one
does, the motif path is strictly better. If none does, drop the features, finish
every other step, and record the omission so it can be revisited.

**Do not skip the decision silently.** A `_MAT` feature declared in the JSON with
no profiles behind it is worse than no feature at all — the annotator will report
it as missing on every genome.

**Getting SignalP, if you go that way.** It is academic-licensed: the tarball is
behind a form the user submits themselves (name, institution, acceptance of the
terms). **Do not attempt to download it for them and do not accept a licence on
anyone's behalf** — walk them through the form and wait.
`scripts/apply_signalp.py --help-install` prints the steps. Once they have it:

```bash
conda create -y -n signalp python=3.10 pip
conda run -n signalp pip install signalp-6-package/
conda run -n signalp pip install "numpy<2"   # torch <2 cannot use numpy 2
conda run -n signalp signalp6 -ff Alignments/G_SIGNALP_QUERIES.faa -od out \
     --organism eukarya --mode fast --format none \
     --model_dir signalp-6-package/models
```

`--model_dir` avoids copying 1.5 GB of weights into the working tree. Without the
numpy pin every prediction dies with `RuntimeError: Numpy is not available`.

**Applying it is where the judgement lives.** A prediction is made from one
consensus sequence; the cut is then applied to every member of the cluster, so it
has to be checked against them:

```bash
python3 scripts/apply_signalp.py --workdir . --signalp out/prediction_results.txt \
        --parent G --self-cutoff 30 --corroborate '[KED]FP[ILM]YTIP'
python3 scripts/apply_signalp.py ... --write
```

The tool refuses a cluster when SignalP's probability is below 0.70, when the
first mature residue is not the modal residue in at least 90% of the cluster,
when any member lacks residues on both sides of the cut, or when the rebuilt
consensus no longer matches what SignalP was given. On Rhabdoviridae that
refused 11 of 53 — the correct outcome, not a shortfall.

`--self-cutoff` is the one that catches a real shipping error. A 15-25 residue
profile can fail to score its own training sequences above the feature's
`bit_cutoff`, which means the annotator could never call it. Two Alpharhabdovirinae
clusters scored 27.5 and 28.4 bits against a cutoff of 30 and were withheld; the
mature chains were still built. Measure with `-comp_based_stats 0` — see the CBS
trap in `references/curation.md`.

**Prefer a measured site to a prediction, and record which you used.** This is
also the whole no-SignalP path: where a mature N-terminus is known
experimentally, or is recorded as a `mat_peptide` in RefSeq, `--motif` cuts by
locating that motif and needs no predictor at all:

```bash
python3 scripts/apply_signalp.py --workdir . --parent G --write \
        --motif 'KFP[ILM]YTIP' --motif-evidence "RABV G, 19 aa signal + 505 aa
        mature, N-terminus confirmed by direct sequencing (PMID 6897030)"
```

This matters beyond provenance. A cleavage site recorded as a column number is
true of exactly one clustering; expressed as a motif it re-locates itself after
any rebuild. The Rhabdoviridae `G_MAT` derivation named ten clusters, and one
rebuild later five had ceased to exist and the other five had different
membership — the profiles were being carried forward against a clustering that
was gone. `--corroborate` marks clusters whose predicted mature N-terminus
matches a measured one, so PROVENANCE distinguishes a prediction that agrees with
an experiment from one standing alone.

**A precursor cleaved more than once** — GPC into signal peptide + Gn + Gc, HA0
into HA1 + HA2 — needs each internal product bounded by *both* neighbours. This
tool derives only the leading product. Align the rest by hand from published
sites and set both extensions to 0; see `cleave=` in `references/json-schema.md`.

### 6b. Handle anything that is not collinear with the genome

Before writing the JSON, settle whether the taxon has a protein that a PSSM
cannot call. Three symptoms, all cheap to check and all expensive to miss:

- **A mature peptide is short by a constant number of residues across most of
  the taxon.** A readthrough stop is cropping the match. Set `internal_stop: 1`
  and rebuild that feature's alignments so they span the stop.
- **tblastn of a protein against its own genome returns two HSPs in different
  frames.** A ribosomal frameshift or an edited transcript. Declare
  `special: transcript_edit`, build a nucleotide reference set, no PSSM.
- **A reference genome annotates a `join(...)` CDS.** Same thing, already
  documented by someone else.

```bash
# the frameshift test, on one genome
tblastn -query protein.faa -db genome -outfmt "6 qstart qend sstart send pident sframe"
# two HSPs, frames 1 and 3 -> -1 frameshift.  1 and 2 -> +1.  One HSP -> collinear.
```

**Do not skip this because the source annotation is silent.** BV-BRC held four
records for alphavirus TF and two of them were mis-annotated 6K-E1 fusions in
frame 0; the protein is nonetheless real, mass-spec confirmed in three species,
and is the form predominantly incorporated into the virion. The source cannot
be used to check this class of work, and its absence is not evidence.

Everything about both cases — phase of the slip site, the trailing `*`
convention, CDS vs mat_peptide, how many references, and which evaluation tools
are blind to the result — is in `references/special-features.md`.

### 7. Write the JSON

**Order matters here and it is not the obvious one.** Collections decide what
the JSON declares (uncharacterized homology groups are discovered from
`collections/UNC*.fasta`), and the JSON decides what the PSSM build targets. So
the cycle is **build collections -> write the JSON -> build features**, not the
reading order of these steps. Building features first simply finds nothing for
any feature the collections just created.

One block per module: `close_genomes`, `segments`, `features`. Set `min_len` and
`max_len` from the collection's own length distribution, and keep `bit_cutoff`
high enough that a genus-specific accessory cannot fire on a sibling genus.
Schema and every field: `references/json-schema.md`.

**Write it in one canonical format** — JSON::XS `pretty` + `canonical`: 3-space
indent, `" : "`, sorted keys. Perl randomises hash order, so without sorted keys
any Perl tool that rewrites the file reorders every line on every run and the
diff is useless. `scripts/json_canon.py --check` verifies it; `--write` fixes
it; `install_module.py` already writes the master this way.

Then check the annotation strings and the citations:

```bash
python3 scripts/check_annotations.py --json <T>_Viral_PSSM.json
python3 scripts/check_dlits.py       --json <T>_Viral_PSSM.json
```

**Never invent a designation.** The `U1` in
`Uncharacterized lineage-specific U1 protein` goes in only where the source
data or a DLIT already uses it for that taxon; keep the tag in the feature key
and gene symbol instead. Expect many features to share
`Uncharacterized lineage-specific protein` — that is correct, not lazy. A
functional claim likewise needs a citation a human has read; if the only DLIT
is model-proposed, fall back to the uncharacterized form and keep the citation
for checking. **Any DLIT a model proposed must be flagged.** A citation you did not read is
not a DLIT, even when the paper is real and on-topic, so it goes in `PMID` and
in `PMID_claude_generated` alongside it. Drive that from a registry in the
generator so it cannot drift, and clear entries only when a human has actually
read the paper. `check_dlits.py` also catches DOIs in the PMID field and ids
that do not resolve in PubMed.

Reuse a vocabulary string wherever one exists — shrinking the vocabulary is the
point of LowVan, and the whole controlled set is only 153 strings. The
annotation records **what the protein does**; the short community name goes in
`gene_symbol`, never in the string ("RNA-dependent RNA polymerase" + symbol
"L", not "L protein"). The vocabulary carries **no taxon prefix** — that was removed
in the current revision, so "Paramyxoviridae C protein" is now just "C protein".
Add genuinely new strings to `assets/annotation-vocabulary.tsv` with Taxon / Annotation / Symbol /
Feature Type / Segment / Used for Genome Quality / PubMed IDs.
`references/annotation-vocabulary.md` covers the conventions.

### 8. QC before installing

```bash
python3 scripts/qc_cross_feature.py --workdir .
```

Finds clusters binned under the wrong feature — a mislabelled source annotation
puts a protein in the wrong collection, and its PSSM then out-scores the correct
one at the *other* protein's locus. On Rhabdoviridae this found nine rabies
glycoproteins labelled "nucleoprotein" in BV-BRC; the resulting profile called
`Nucleocapsid protein` on the G coordinates and N was never called at all.
Length gates do not catch this class — 524 aa is inside the range a real
nucleocapsid can occupy.

### 9. Install

```bash
python3 scripts/install_module.py --workdir . --repo $LOWVAN_DATA_DIR --check
python3 scripts/install_module.py --workdir . --repo $LOWVAN_DATA_DIR
```

`--check` refuses to install a module whose JSON and PSSMs disagree. A feature
declared with no PSSM is never called; a PSSM with no JSON entry is never
loaded. Both are silent at runtime.

### 10. Run it and score it

```bash
python3 scripts/check_rep_contigs.py --repdir Rep-Contigs --acc-file accs.txt
python3 scripts/validate_calls.py   --repo $LOWVAN_DATA_DIR --fasta genome.fna
python3 scripts/evaluate_module.py  --repo $LOWVAN_DATA_DIR --acc-file accs.txt
```

`validate_calls.py` asks whether one genome's output is *internally* sound —
does every protein start with M, is any truncated by an internal stop, is the
length inside what the JSON declares, did two features land on the same
coordinates. It needs no reference annotation, so it works on a genome nobody
has annotated before. `evaluate_module.py` then asks the different question of
whether the calls agree with GenBank across a held-out set.

Check routing first. Routing is a **nucleotide** BLASTn against rep contigs, so
a divergent group needs more of them than you expect — one per genus was not
enough for the plant rhabdoviruses, and an unrouted genome is rejected before
any PSSM runs.

Then score the calls on held-out genomes. Two outcomes matter more than the
recall percentage:

- **duplicates** — two features on identical coordinates. Always a defect; go
  back to *Write the JSON* and check `cleave=` on the mature peptide, then
  rerun `qc_cross_feature.py` and `filter_unchar.py`.
- **misses** — a reference CDS with no call. Either an unbuilt feature or a
  `bit_cutoff` set too high.

### 11. Evaluate coverage over the whole taxon

Step 10 scores a handful of reference genomes. That is a correctness check, and
it flatters the module: the panel is small, curated, and usually made of the
same genomes the profiles were built from. Before shipping, run the module over
everything BV-BRC holds for the taxon.

```bash
# download every genome in the length window, and stop
perl Other_Scripts/New-annotate-viral-taxon.pl -i <Taxon> \
     -download-only -d Contigs -metaout <Taxon>.metadata -min 1000 -max 50000

# cluster, pick exemplars, annotate, report
python3 scripts/evaluate_coverage.py --contigs Contigs \
        --repo $LOWVAN_DATA_DIR --min-complete 10000 --out coverage.json
```

**Merge segmented genomes first.** A segmented virus is stored one record per
segment, so every genome counts as several and a fully annotated one looks like
a handful of partials. 56 Dichorhavirus records are 28 genomes; scoring them as
56 halves the apparent completeness of a module that did nothing wrong.

```bash
python3 scripts/merge_segments.py --json <T>_Viral_PSSM.json --fasta-dir Contigs \
        --metadata <T>.metadata --out Contigs-merged/
```

Merging is gated on the segment count the module declares (`"segments": N` in
the module block). More parts than declared is a name collision between
isolates; fewer is a genuinely incomplete genome. Only an exact match merges —
everything else passes through untouched and is logged.

**Cluster before annotating.** BV-BRC holds thousands of near-identical genomes
for the well-sequenced species — 25,237 Rhabdoviridae records collapse to 1,983
clusters at 95% nucleotide identity. Annotating all of them measures sequencing
effort, not module coverage. One exemplar per cluster gives each distinct genome
type one vote.

**Restrict to near-complete genomes for the protein-distribution question.**
Most records in a family-wide download are single-gene partial submissions: the
median Rhabdoviridae record is 1,575 bases, and only 9,007 of 25,237 reach
10 kb. Partials are flagged poor for missing essential features, which is
correct and says nothing about the module. Of the 1,983 clusters, 657 had no
near-complete member at all.

Four questions, and what each one tells you:

**1. How many genomes does BLASTn routing miss — and can it be fixed?** It is
not enough to count them. Work out the *minimal set of new reference contigs*
that would close the gap, because that turns a complaint into a work order.

```bash
python3 scripts/minimal_refs.py --unrouted unrouted.ids --clusters clu_cluster.tsv \
        --json <T>_Viral_PSSM.json --metadata <T>.metadata --max 10
```

Greedy set cover: take the exemplar routing the most still-unrouted genomes,
remove what it covers, repeat. The script excludes any candidate whose taxon
the module has no profiles for — **a reference contig for an uncovered taxon is
worse than none**, because it routes the genome to a module that then calls
nothing, converting a visible gap into a silent one.

It prints the curve and stops there. **The cutoff is yours**: the shape of the
curve decides whether the gap is worth closing, and how many references you are
willing to carry forever is a judgement the script cannot make.

On Rhabdoviridae the curve was almost flat — 5 references bought 17% of the
gap, 50 bought 44%, and full coverage needed 243. That is the signature of a
long tail of singletons: 200 of 346 unrouted genomes hit no other unrouted
genome at all. You cannot reference-your-way out of that, and trying would
double the rep-contig set to chase the last few percent.

So this is a judgement call about where the curve stops paying: take the
references whose marginal coverage is still large, stop when it flattens. The
first one alone covered 31 genomes; by the tenth the increment was down to two
or three. Add the head of the curve, leave the tail, and record the decision.

Also check the opposite lever before adding anything. Blast the unrouted
genomes against the **existing** rep contigs and look at the best-score
distribution. If they cluster just under the threshold, lowering
`-min_contig_bit` is one cheap change instead of many new files. If they have
no hit at all, lowering it buys nothing. On Rhabdoviridae 313 of 346 had **no
BLASTn hit whatsoever** to the existing references, so dropping the threshold
to 50 would have recovered 10% and loosened routing for everything else — a bad
trade, and only measurement showed it.

**Beware the genomes that hit a reference and still produce nothing.** A genome no rep contig claims
is never annotated at all, however good the profiles are. This is the failure
that hid Rice yellow stunt until a sixth Betarhabdovirinae rep contig was added.
Break the misses down by genus: scattered singletons are divergent outliers,
but a cluster of them in one genus is a missing rep contig, and that is cheap
to fix.

**2. Good versus poor quality** — and this must come from the GTO, never from
the feature table.

`annotate_by_viral_pssm.pl` writes a flat table, and reimplementing the quality
rules against it looks easy: compare each row's length to `min_len`, count
copies. It is wrong twice over.

The table is one row per called region, not per feature. A polymerase broken by
a frameshift returns several rows with separate feature ids — L was fragmented
in **18%** of Rhabdoviridae genomes, G in 5% — so comparing each fragment to the
full-length floor produces a flood of "Feature is too short" that is an artefact
of the comparison, not a property of the genome.

Worse, the flat table cannot represent **special features**. Splice variants and
transcript-edited products are absent from it, so a taxon that has them is
scored as though those features do not exist. Rhabdoviridae happens to ship no
special features in its four modules, which is luck, not design — the shortcut
appeared to work here and would have failed silently on Paramyxoviridae.

So run the real pipeline:

```bash
python3 scripts/make_gto.py --fasta-dir Contigs-merged \
        --metadata <T>.metadata --out gto --jobs 24
python3 scripts/run_gto_eval.py --gto-dir gto --repo $LOWVAN_DATA_DIR \
        --out gto_out --jobs 8 --report quality.tsv \
        --perl /Applications/BV-BRC.app/runtime/bin/perl \
        --perl5lib /Applications/BV-BRC.app/deployment/lib
```

The Quality GTOs carry `genome_quality`, `feature_quality` and
`feature_quality_flags` as the annotator computed them.

**Use `rast-create-genome`; never hand-write the GTO JSON.** It issues the
genome id from the ID server, and `new_feature_id` numbers every feature the
annotator adds off that id, so a GTO carrying a reused or invented id produces
feature ids that collide with a real genome's. `make_gto.py` fails loudly rather
than falling back. Most of each call is the round trip to the ID server, so
`--jobs 24` is reasonable.

The BV-BRC side of this is **bundled** in `vendor/bv-brc/`, so no BV-BRC
installation is required; `make_gto.py` finds the vendored `rast-create-genome`
on its own and refuses to start if the perl it will use is missing `File::Slurp`
or a UUID module. Three things still go wrong, and all three are flags rather
than edits:

| symptom | cause and flag |
|---|---|
| `Can't locate GenomeTypeObject.pm` / `IDclient.pm` | the kit bundles these in `vendor/bv-brc/lib`; `make_gto.py` finds them itself. Only needed if you are calling the perl directly — `--perl5lib vendor/bv-brc/lib` |
| `No UUID generator found` | `GenomeTypeObject` needs `Data::UUID` or `UUID`. The BV-BRC runtime perl has one; a conda perl usually has neither — `--perl <kit>/runtime/bin/perl` |
| `rast-create-genome: command not found` | the macOS bundle ships `plbin/rast-create-genome.pl` without its `bin` wrapper; `make_gto.py` writes it from the bundle's own pattern |

The `--perl` split matters: only the GTO wrapper and the quality script run
under the dev-kit perl. The inner `annotate_by_viral_pssm.pl` resolves through
`PATH`, so it keeps the interpreter its own dependencies are installed under.

**3. Which flags fire on the poor ones.** This is the audit of `min_len`,
`max_len` and `copy_num`. "Feature is too short" dominating means the length
bounds are wrong, not the genomes — see the percentile rule in
`references/json-schema.md`.

**4. The distribution of called proteins.** In complete genomes the core
proteins should appear in near-equal numbers, because every genome has one of
each. A core protein lagging the others is a coverage hole in that feature's
profiles, and the gap size estimates how many genomes it costs.

Question 4 is the one that catches what nothing else does. Self-recall says a
profile finds its own training data; per-genome counts say whether it finds the
protein in genomes nobody curated.

### 12. Publish the artifacts

```bash
python3 scripts/collect_synmap.py .                       # then gen the collapse page
python3 scripts/collect_registry.py --workdir . --out registry.json
python3 scripts/gen_registry.py --registry registry.json --taxon <Family> \
    --notes notes.json --out registry.html
```

- **String collapse** — how many BV-BRC strings each annotation absorbs. This is
  the number the LowVan manuscript reports.
- **Vocabulary rarefaction** — the same claim as a curve, and the more
  convincing form of it. Shuffle the genomes, accumulate distinct annotation
  strings, average over replicates. A controlled vocabulary saturates, because
  the Nth genome reuses names the first N-1 established; free text climbs
  roughly linearly, because every submitter spells the same protein a new way.

  ```bash
  python3 scripts/annotation_rarefaction.py --new coverage_eval/ann \
          --old bvbrc_products.tsv --out rarefaction.json --replicates 100
  ```

  Report the ratio, not just the picture: strings per 100 genomes at the end of
  each curve says how much collapse the vocabulary actually bought.
- **PSSM registry** — every declared feature, its profiles, and what fraction of
  its own collection it recovers.

Write the `--notes` prose. Auto-generated counts without commentary explain
nothing. `references/artifacts.md` has the note schema and what to say.

### 13. Chase any feature losing 20% or more

Step 11 gives you a per-protein call rate. Where one core protein trails the
others by 20 points or more, **do not accept it and do not guess at the
cause.** The gene is either absent from those genomes or invisible to the
profiles, and those demand opposite responses. The check is cheap and it is
decisive.

Take the genomes that called the protein's two genomic neighbours but not the
protein itself, and look in the interval between them.

```bash
python3 scripts/analyze_feature_gap.py --workdir . --module <Module> \
        --key P --left N --right M \
        --ann-dir coverage_eval/ann --genome-dir coverage_eval/ex \
        --save-extracted
```

It selects the genomes, cuts the interval, ORF-calls it in six frames, blasts
the longest ORF against the feature's collection and then against
`Leftover_Seqs.aa` versus the clustered members, and prints a verdict. What it
cannot do is pick the flanking features: **you supply `--left` and `--right`
from the genome organisation**, which is taxon-specific knowledge and does not
transfer.

Four outcomes, four different conclusions:

| what you find in the interval | means |
|---|---|
| interval short, or full of ambiguous bases | genome quality — nothing to fix |
| no ORF | the gene really is absent |
| ORF present, matches the collection, closest to a **clustered** member | profile threshold or coverage cutoff is wrong |
| ORF present, matches the collection, closest to a **leftover** | the profile set does not cover that variant |

The last one is the common case and the easy fix: build leftover profiles
(`build_leftover_pssms.py`, step 3) rather than touching any threshold.

Worked example, Alpharhabdovirinae P at 68% against N at 99% — a 31-point gap:

```
N->M interval present           115/115   median 1018 nt   (P needs ~600-1100)
ambiguous bases                   8/115
contains an ORF                 115/115   median  296 aa   (P bounds 174-359)
ORF matches the P collection    109/115   median bitscore 513, 82 above the 80-bit cutoff

closer to a LEFTOVER P            85
closer to a CLUSTERED P            0
```

The prior expectation was genome quality. It was wrong: every gene was there,
full length and unambiguous, and **not one** resembled a sequence a profile was
built from. The fix was 19 more profiles from sequences already in the
collection, not a threshold change and not more reference contigs.

ORF-calling also recovers proteins the source does not have at all — 21 of
those 115 phosphoproteins appear in no BV-BRC protein record. Put them in
`Extracted/<Module>.<KEY>.faa`, which `build_collections.py` merges into the
collection, so the recovery survives the next rebuild.

**Then go back to step 3 and rebuild the feature.** This is the loop that ends
the work: measure, find the gap, learn which of the four things it is, fix that
one, measure again.

### 14. Verify that every cleaved end is actually crisp

`cleave=` sets `upstream_ext` and `downstream_ext` to 0, which tells the
annotator not to scan outward for a start or a stop. That is right for a
protease cut site, whose position no codon marks — and silently wrong for any
terminus that is not one, where the boundary then drifts with whatever the
tblastn alignment happened to do.

```bash
python3 scripts/check_cleaved_ends.py --json <T>_Viral_PSSM.json \
        --ann-dir coverage_eval/ann --precursor G --products G_SP,G_MAT
```

A clean cut is a spike at one offset. A spread is creep. An offset of +3 at a
precursor's C-terminus is its stop codon and is expected.

Worked example: Alpharhabdovirinae `G_MAT` carried `downstream_ext: 0`, which
asserts its C-terminus is a cut site. It is not — signal peptidase removes the
leader and nothing else, so the mature chain runs to the precursor's own stop.
Over 1,108 genomes, **40% of mature G calls ended somewhere other than the
precursor's stop**, median 28 codons short and out to 1.2 kb. Setting
`downstream_ext: 1` fixed it.

**Check this every time, and check it hardest on polyproteins.** A precursor
with three internal products has six cleaved termini, and every one left to
drift compounds along the chain. Only genuine protease sites keep a 0; give
every other terminus a 1 and let the annotator find the codon.

## Reference files

| File | Read it when |
|---|---|
| `references/inputs.md` | taking delivery of the BV-BRC dump |
| `references/module-partitioning.md` | deciding what the modules are (step 2, before rules) |
| `references/annotation-triage.md` | writing classification rules |
| `references/annotation-vocabulary.md` | choosing annotation strings |
| `references/pipeline.md` | first pipeline run; parameter tuning; upstream bugs |
| `references/curation.md` | curating alignments; splitting; mat_peptides |
| `references/json-schema.md` | writing the JSON block |
| `references/special-features.md` | a protein is short by a constant amount, or is not collinear with the genome |
| `references/install-and-test.md` | installing, running, scoring |
| `references/artifacts.md` | writing the two reports |

## Scripts

All take `--workdir` pointing at the module working directory, which looks like:

```
<workdir>/
  <T>.id_md5  <T>.id_name_gs  <T>.uniq.{md5,id_ann,seq}   inputs
  <T>_Viral_PSSM.json                                     the module JSON
  collections/<Module>/<FEAT>.fasta                       binned sequences
  Alignments/<Module>/<FEAT>/{corrected_alis,pssms}/      pipeline output
  Rep-Contigs/<Module>.<n>.dna                            routing subjects
  Reports/                                                artifacts
```

| Script | Purpose |
|---|---|
| `rebuild_pssms.py` | find alignments edited since their PSSM was built; rebuild |
| `build_sp_pssms.py` | derive per-cluster signal-peptide PSSMs from a parent CDS and its mature peptide |
| `evaluate_coverage.py` | cluster the taxon's whole BV-BRC holding, annotate one exemplar per cluster, and report routing, quality and protein distribution |
| `apply_signalp.py` | stage consensus N-termini for SignalP, then cut verified `_SP` / `_MAT` products from a parent CDS |
| `norm_pssm.py` | normalise `psiblast -out_pssm` output to the pipeline's format |
| `qc_cross_feature.py` | find clusters binned under the wrong feature |
| `check_annotations.py` | check annotation strings and gene symbols against the vocabulary |
| `check_dlits.py` | verify DLIT citations resolve, and that model-proposed ones are flagged |
| `check_dump.py` | verify the BV-BRC dump is complete and line-aligned |
| `json_canon.py` | one canonical JSON format so diffs show only real changes |
| `check_rep_contigs.py` | check that rep contigs route the taxon's genomes |
| `install_module.py` | validate and install into a Viral_Annotation checkout |
| `validate_calls.py` | is one genome's output internally sound? no reference needed |
| `evaluate_module.py` | run the annotator on held-out genomes and score |
| `collect_registry.py` | measure per-feature self-recall |
| `gen_registry.py` | render the PSSM registry artifact |
| `collect_synmap.py` | gather data for the string-collapse artifact |
| `build_leftover_pssms.py` | second-pass profiles from `Leftover_Seqs.aa` at a member floor of 2 |
| `filter_unchar.py` | split uncharacterized proteins into `UNC1..UNCn` by homology |
| `analyze_feature_gap.py` | why is this feature called less often than its neighbours? |
| `check_cleaved_ends.py` | are `cleave=` termini crisp, or is the alignment drifting? |
| `merge_segments.py` | reassemble multi-segment genomes split one record per segment |
| `minimal_refs.py` | greedy set cover for the fewest reference contigs that close a routing gap |
| `annotation_rarefaction.py` | vocabulary growth curve, controlled versus free text |
| `make_gto.py` | contig FASTAs -> GTOs, with genome ids issued by the ID server |
| `run_gto_eval.py` | annotate and quality-score GTOs; the good-vs-poor breakdown |

## What "done" looks like

- every module in the JSON has PSSMs, or is deliberately declared and recorded
  as unbuildable with the reason
- `qc_cross_feature.py` reports no MISLABEL rows
- `rebuild_pssms.py` reports 0 stale
- `install_module.py --check` reports no problems
- `json_canon.py --check` reports every JSON canonical
- `check_rep_contigs.py` routes every test genome
- `evaluate_module.py` reports **0 duplicates**, and the misses are all
  features you know you did not build
- `check_cleaved_ends.py` shows a single spike at every `cleave=` terminus
- no feature trails its neighbours by 20 points with the cause unexplained
- the artifacts published, with prose

If you cannot call all of a taxon's proteins, you do not have a module for that
taxon. Say so and record why, rather than shipping a module that silently
misses half a genome.
