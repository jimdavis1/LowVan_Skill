# Annotation strings

Every string a module emits becomes a permanent annotation in BV-BRC. LowVan
exists to *shrink* the annotation vocabulary — the manuscript reports 1.5- to
19.4-fold reductions in unique strings, and the whole controlled set is 153
strings — so coining a near-duplicate of an existing string works directly
against the point of the project.

The canonical list is **S1 Table** of the LowVan manuscript, shipped here as
`assets/annotation-vocabulary.tsv` (478 rows, 35 taxa, 153 distinct strings) and
`assets/S1-Table.xlsx`.

```bash
python3 scripts/check_annotations.py --json <T>_Viral_PSSM.json
```

Columns: `Taxon | Annotation | Symbol | Feature Type | Segment | Used for Genome
Quality | PubMed IDs*`.

## The governing rule

From the manuscript:

> In LowVan, the functional description is recorded as the annotation string,
> while the commonly used symbolic name is recorded as the gene symbol. For
> example, the measles virus polymerase is annotated as "RNA-dependent RNA
> polymerase" with the gene symbol, "L."

So the annotation says **what the protein does**; the symbol is a separate
field. `RNA-dependent RNA polymerase` + symbol `L` — never `L protein` as the
annotation. In S1 Table the symbol lives in its own column; in the module JSON
it lives in `gene_symbol`.

## Rule one: reuse before you coin

If S1 Table already has a string for this protein, use it verbatim. The core
mononegavirus proteins are all there and shared across families:

```
Nucleocapsid protein          Matrix protein            Phosphoprotein
RNA-dependent RNA polymerase  Movement protein          Small hydrophobic protein
Fusion glycoprotein           Receptor binding protein
C protein                     V protein                 W protein
Uncharacterized lineage-specific protein
```

Rhabdoviridae reuses nine of these unchanged.

## No taxon prefix

**The current revision of S1 Table removed the family and genus prefixes.**
Earlier versions carried them; they are gone, and adding one now creates a
duplicate of an existing string:

| Was | Is now |
|---|---|
| `Paramyxoviridae C protein` | `C protein` |
| `Pneumoviridae small hydrophobic protein` | `Small hydrophobic protein` |
| `Filoviridae matrix protein VP40` | `Matrix protein VP40` |
| `Sarbecovirus cytotoxic interferon antagonist ORF6` | `Cytotoxic interferon antagonist ORF6 protein` |
| `Fimoviridae uncharacterized protein` | `Uncharacterized lineage-specific protein` |

The taxon is already recorded in the Taxon column and in the module key, so
repeating it in the string only adds strings.

## Conventions for new strings

**Sentence case.** `Envelope glycoprotein precursor`, not Title Case and not
lower case.

**Numbered proteins read `<description> <SYMBOL> protein`.** The symbol goes
before the trailing word, not after the description:

```
Helicase NSP13 protein                    RNA-binding NSP9 protein
Main protease involved in polyprotein cleavage NSP5 protein
Interferon antagonist NS1 protein         Transcription activator protein VP30
RNA silencing suppressor P7 protein       Virulence factor PB1-F2 protein
```

**Unnamed ORFs are `Uncharacterized lineage-specific <symbol> protein`.** Or
`Uncharacterized <symbol> protein` where the protein is not lineage-restricted.
Both are real S1 Table forms. A bare `Uncharacterized protein` is too generic.

**Mature peptides begin `Mature `.**

```
Mature envelope glycoprotein              Mature C-terminal envelope glycoprotein Gc
Mature nonstructural membrane protein     Mature secreted nonstructural 38-kDa protein
```

**Signal peptides name their parent: `Signal peptide of <symbol>`.** S1 Table
has `Signal peptide of F0` and `Signal peptide of GPC`; a bare `Signal peptide`
is no longer used.

**No qualifiers.** The manuscript is explicit that LowVan does not add
`incomplete`, `truncated`, `partial` and the like, precisely to stop the
vocabulary proliferating. A truncated protein gets the same string as the
full-length one.

**Say what it does when you know.** `Transcription elongation factor M2-1` beats
`M2-1 protein`. Where a function is unknown, the manuscript's stated practice is
to use the number alone as both annotation and symbol, with the intent of adding
a function later.

## Never invent a designation

A **designation** is the short positional tag inside a string: the `U1` in
`Uncharacterized lineage-specific U1 protein`, the `Gx`, the `Px`, the `alpha1`.

**Put one in an annotation string only when the source data or a DLIT already
uses it for that taxon.** A reader has no way to tell an invented `U1` from a
real one, and the number implies a correspondence between taxa that usually
does not exist.

This project got it wrong first time round. BV-BRC genuinely uses U-numbers for
four rhabdovirus genera — Hapavirus U1/U2/U3, Tibrovirus U1/U2/U3, Sripuvirus
U1/U2, Sunrhavirus U3 — traceable to Walker 2015. That nomenclature was then
extended to twelve genera that use nothing of the kind: Arurhavirus calls its
ORFs `pAG1`/`pAG2`, Stangrhavirus calls them `VP1`/`VP2`/`VP3`, Almendravirus
calls one `SH protein`. One invented `U2` was assigned to a genus that has a
real `U3` and no `U2` at all. Thirty-five annotation strings carried a tag
nothing supported.

Keep the tag in the **feature key** and the **gene symbol** — they are internal
identifiers and have to be unique, and S1 Table namespaces symbols for exactly
this reason. It is only the published annotation string that must stay honest.

The consequence is that many features share
`Uncharacterized lineage-specific protein` — 33 of them here. **That is the
correct outcome, not a failure.** They are uncharacterized lineage-specific
proteins, the Taxon column says which lineage, and the feature key still tells
you which ORF. Inventing distinct names to avoid the repetition would be
inventing information.

Drive it from an explicit allowlist so the rule is auditable:

```python
BACKED_DESIGNATION = {"HAPA_U1", "HAPA_U2", "HAPA_U3", "TIBRO_U1", ...}
```

A designation earns its place by appearing in the source data or in a DLIT that
names it. Anything else gets stripped at generation time.

## A functional claim needs a DLIT someone has read

The same principle one level up. `RNA silencing suppressor protein` asserts a
function. If the only citation behind it is model-proposed — real paper, right
topic, nobody read it — the claim is not established and the annotation falls
back to the uncharacterized form. The citation stays on the feature so it can
be checked and the name restored.

Two exemptions, or the rule eats the vocabulary:

- **Established vocabulary is exempt.** Calling a 2,100 aa rhabdovirus ORF
  `RNA-dependent RNA polymerase` rests on no particular paper. Any string
  already in S1 Table, or in the module's own controlled set, passes.
- **A backed designation is exempt**, since the source data already makes the
  claim.

Applied here it downgraded exactly two features, which is the right size for a
rule like this: one viroporin claim and one silencing-suppressor claim, both
resting solely on citations a model proposed.

## Gene symbols, and when to namespace them

The symbol is the community's short name, in its own column: `L`, `N`, `SH`,
`NSP13`. Not a description, and not the annotation repeated.

When one module collapses several subgenera that each own a same-numbered ORF,
S1 Table prefixes the symbol with the lineage to keep them apart. 41 of its 442
symbols do this:

```
Embeco_NS2a   Merbeco_ORF4a   Hibeco_ORF3   B_COV_ORF6   G_COV_ORF4a   D_COV_NS7b
```

So a collapsed module should namespace exactly where it needs to disambiguate
and not otherwise. Alpharhabdovirinae holds ~38 genera and several of them have
their own SH or U1, so `Sunrha_SH` and `Tupa_SH` are correct; a bare `SH` would
collide. `check_annotations.py` recognises the `<Lineage>_<SYM>` form, strips
the prefix before comparing, and reports those separately from a genuine
disagreement.

A genuine disagreement is when the base symbols differ — your `P7` against S1
Table's `M` for `Matrix protein`. Either is defensible: the manuscript says
"historically accepted gene symbols and names commonly used by the viral
research community were retained", so a positional symbol that the community
actually uses for that virus can stay. Decide it deliberately rather than by
accident.

## Feature Type, Segment, quality flag

`Feature Type` is `CDS` for full-length proteins and `mat_peptide` for cleavage
products — the manuscript states this explicitly, and says that when multiple
cleavage products are possible, all are annotated.

`Segment` is `Single RNA Segment` for nonsegmented viruses, otherwise the
segment name exactly as used in the JSON `segments` block.

`Used for Genome Quality` is 1 if the feature's presence and copy number should
count toward the genome quality score. Set 0 for accessories genuinely absent
from many members of the module, or every genome lacking an optional protein is
marked poor quality.

## Adding to the Excel

New strings go into `S1-Table.xlsx` as new rows, one per taxon that uses the
string, with all seven columns filled. Put the DLITs — the papers defining the
protein's function or sequence — in `PubMed IDs*`, and the same IDs in the JSON
`PMID` field. A new string with no citation is a liability for whoever maintains
the module next.

**Only human-verified citations belong in the Excel.** S1 Table has no
provenance column, so a model-proposed PMID entering it silently becomes a
curated one. The JSON can carry the distinction (`PMID_claude_generated`, see
`references/json-schema.md`); the Excel cannot. Verify it or leave the cell
empty.

## What `check_annotations.py` reports

- **REUSED** — already in S1 Table. The good outcome.
- **GENE SYMBOL namespaced** — your symbol is `<Lineage>_<SYM>` and the base
  matches S1 Table. That is S1 Table's own convention for a collapsed module;
  no action.
- **GENE SYMBOL differs** — the base symbols genuinely disagree (`P7` against
  S1 Table's `M`). A real decision, not a typo; see the section above.
- **VARIANT** — the same string but for case, punctuation or spacing. Always a
  defect. The check is deliberately narrow and does not use edit distance: this
  vocabulary distinguishes proteins by a single character (`P5` / `P6` / `P7`,
  `NS7a` / `NS7b`), so any tolerance produces only false positives.
- **NEW** — needs adding to the Excel, and gets a style check.
