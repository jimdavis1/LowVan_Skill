# Worked example: Matonaviridae

The third worked example, and the one to read if **the source annotation is the
problem**. Rhabdoviridae covers the ordinary path; Togaviridae covers proteins a
PSSM cannot call. Everything novel here is about a taxon whose BV-BRC export is
largely mislabelled, fragmentary or in the wrong reading frame — and about
proving a binning is right without trusting that export.

```
build_collections.py           source annotation strings -> feature collections
                               (+ homology resolution of generic strings,
                                + a corrupt-record filter)
build_matonaviridae_json.py    collections + measured cutoffs -> the module JSON
measure_cutoffs.py             the signal/noise window per feature
check_tiling.py                PROVE the mature-peptide binning from sequence alone
qc_polyprotein.py              cross-feature QC that knows precursor from product
download_contigs.py            the taxon's contigs (the kit's -download-only ships broken)
run_pipeline.sh                the clustering/alignment/PSSM driver
gen_rarefaction.py             render rarefaction.json as a report page
validate_palette.py            the dataviz colour checks, without a JS runtime
annos.tsv                      feature key -> annotation string
```

`build_collections.py` and `build_matonaviridae_json.py` are **taxon-specific by
design** — copy them, replace the tables at the top, leave the machinery below
alone. In `build_collections.py` the `WINDOW`, `NOT_MODELLED` and `rules()`
blocks are Matonaviridae; `resolve_by_homology`, `drop_corrupt` and everything
after are general.

Run with `LOWVAN_KIT` pointing at this clone and `LOWVAN_WORK` at your build
directory.

## The taxon in one paragraph

Three species, one genome organisation, seven features. Two ORFs: a
nonstructural polyprotein cut **once** into p150 and p90, and a structural
polyprotein cut by host signal peptidase into capsid, E2 and E1. In all three
species the mature peptides sum **exactly** to their precursor — rubella
1301+815=2116 and 300+282+481=1063; Rustrela 1082+828=1910 and 333+324+490=1147;
Ruhugu 1237+811=2048 and 317+295+486=1098.

Nothing here is bespoke. No `-x` build, no extracted feature set, no special
feature, and no SignalP: `run_pipeline.sh` at the documented defaults builds all
seven features from the export alone.

## What this example adds

### The feature set comes from papers, not from `mat_peptide` records

Both decisions that shape this module are literature decisions, and both would
have gone the other way if read off the database:

- **p150 is a terminal product, not a precursor.** It carries a
  methyltransferase, a Y domain, a proline hinge, an X (macro) domain *and* the
  papain-like protease, which invites five features. There is exactly one
  cleavage of p200 (PMID 9557742, "at a single site"; PMID 10823845, "two mature
  products … NH₂-p150-p90-COOH"; PMID 10799588 maps the scissile bond by
  mutagenesis), and the NS-protease review states the negative outright: *no
  other cleavages within P200 have been detected.*
- **There is no signal-peptide feature, and SignalP is never needed.** Both
  structural signal sequences are *retained on the upstream product* — the E2
  signal sequence stays on the capsid's C-terminus in the mature virion
  (PMID 2214022, by antipeptide antiserum; PMID 11160697). Nothing is removed,
  so C/E2/E1 already tile p110 and there is no `_SP` product to derive.

Naming follows the same test. p150's protease activity is demonstrated, its
methyltransferase is assigned by homology, so the string is
`Protease p150 protein` and not `Methyltransferase and protease p150 protein`.
p90's polymerase is demonstrated and its helicase is homology-assigned, so it
reuses the existing controlled string `RNA-dependent RNA polymerase` and spends
no new vocabulary. Six of seven strings are reused; one is new.

### `check_tiling.py` — prove the binning without the source annotation

The strongest check in this build costs ten lines and needs no reference
annotation. If the mature peptides really are the products of the precursor,
then in any genome carrying both, the precursor must **equal** the concatenation
of its products:

```
p200 = p150 + p90      114 exact, 1 mismatched   (of 115 genomes)
p110 = C + E2 + E1      71 exact, 0 mismatched   (of  71 genomes)
```

This matters because the source annotation is the thing being replaced, so it
cannot also be the thing that validates the replacement. The termini then
reproduce the published chemistry independently: p150 ends `...PLSRGG` in 161
sequences and p90 begins `GTCAAT` in 87 — the SRGG/GTCA site of PMID 10799588 —
and E2 ends `...PAYG`, which is the `278`PPAY`281` late domain at exactly the
residue a 282-aa E2 puts it (PMID 39207105).

**Run this for any module that declares mature peptides inside a CDS.**

### Generic strings resolved by homology, not by length

A quarter of this export's vocabulary names no protein: `envelope
glycoprotein`, `structural protein`, `polyprotein`, `hemagglutinin`,
`NS4 protein` — 227 sequences. Length cannot resolve them either, because the
200–299 aa band holds **both** mature E2 (282 aa) and the 246-aa WHO E1
genotyping amplicon, which is the single most common record in the taxon.

So `resolve_by_homology` blasts them against the sequences the specific strings
already binned confidently and assigns by best hit, writing every call with its
bitscore to `GENERIC_RESOLVED.tsv`. 226 of 227 resolved. Two results worth
seeing: `hemagglutinin` → E1 (correct — rubella E1 *is* the haemagglutinin), and
`NS4 protein` → p90 at 505 bits, which is the U-number trap appearing in a
family that has no NS4 at all.

Guessing from length would have mis-binned 55 partial E1 fragments as E2.

### `drop_corrupt` — because corrupt records build convincing profiles

BV-BRC holds **126 Matonaviridae records whose translation is in the wrong
frame**. Two shapes, and only the first announces itself:

1. internal stop codons — 116 sequences, up to 8 stops each
2. no stops, full length, wrong frame — seven E1 records beginning
   `GFHLPLYCTGVRHADT` where every genuine rubella E1 begins `EEAFTYLCTAPGCATQ`

In the first pass **both shapes built PSSMs**: `E1.lo1` from five 8-stop
sequences and `E1.3` from the frameshifted seven. Both looked healthy — each
recovered its own members at ~1000 bits, which is what a good profile does. The
tell was that each recovered **nothing else**, and that together they were
precisely the twelve sequences the dominant E1 profile could not see.

Nothing else in the workflow catches this. `qc_cross_feature.py` does not, the
truncation QC does not, self-recall actively hides it, and a small reference
panel would never show it. The filter is two tests: reject an internal stop, and
reject anything below 100 bits against the feature's own modal-length core.
That floor is safe — genuine divergent species clear it easily (Rustrela E1 513
bits, Ruhugu 596) while a wrong-frame translation scores nothing.

### A species that cannot have a profile, and is called anyway

Ruhugu virus has two genomes that are **sequence duplicates**, so after dedup
there is one unique sequence per feature — below any clustering floor, including
the leftover floor of 2. No parameter reaches it and none was tried.

It is called regardless, by the rubella profiles, at 246–2379 bits, because it
is rubella's closest known relative. The cost is specific and worth knowing:
`C.downstream_ext` and `E2.upstream_ext` are `0` because that junction is a
signalase cut, so the boundary is wherever the calling profile's alignment
ended. On rubella and Rustrela that is exact; on Ruhugu the C/E2 junction lands
**54 codons early**, while p200, p150, p110 and E1 are called exactly. Not
fixable by a threshold — it needs a second distinct *R. ruteetense* sequence to
exist.

### Two routing measurements that were wrong before they were right

Worth reading before you measure routing yourself, because both wrong answers
looked entirely reasonable.

**First:** plain `blastn`, whose default task is megablast with a ~95% identity
floor. It reported all 13 divergent records rejected. The annotator searches with
`-word_size 11 -reward 2 -penalty -3 -soft_masking false -evalue 0.5`
(`annotate_by_viral_pssm.pl:209`), which is far more sensitive.

**Second:** after fixing the flags, one `blastn` call with all 681 genomes as
subjects — where `-max_target_seqs` defaults to **500**, so hits past the cap
vanished. That invented a gap of 67 "unrouted" rubella partials, complete with a
plausible cluster of scores just under the 150-bit floor, detailed enough to
plan work around.

Ground truth is the annotator, one genome at a time: **668 of 668 rubivirus
genomes route, at every length band** (217/217 near-complete, 49/49, 19/19,
383/383). There was never a gap. The nine rejections are the deliberately
excluded taxa, which is the outcome the exclusion decision requires.

A routing measurement that is not the annotator is not a routing measurement.

### `qc_polyprotein.py` — containment is not mislabelling

`skill/scripts/qc_cross_feature.py` flags any two features sharing near-identical
sequence. In a polyprotein taxon that is the architecture: every mature peptide
is a substring of its precursor. Matonaviridae produced 13 such rows and
Togaviridae 174, **none of them mislabels**, and a real one would be invisible
among them. This wrapper subtracts the pairs the genome organisation predicts and
reports only what is left (0 here).

Two taxa have now needed the same wrapper, which is an argument for teaching
containment to the kit script directly.

## Read the comments

As with the other examples, much of the value is in the comments recording why a
rule exists and what measurement settled it: why length cannot separate p150
from p110 (their windows overlap across 1000–1200 once Rustrela and Ruhugu are
included), why fourteen species prefixes are excluded from training but not from
the module's scope, and why the p110/capsid second in-frame AUG — real, and used
— gets no feature of its own.
