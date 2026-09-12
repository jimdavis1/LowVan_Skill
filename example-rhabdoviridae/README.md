# Worked example: Rhabdoviridae

The two programs that do the central work of a module build. Both are
**taxon-specific by design** — copy them, replace the rule tables at the top,
and leave the machinery below alone.

```
build_collections.py          source annotation strings -> feature collections
build_rhabdoviridae_json.py   collections + curated alignments -> the module JSON
query_PATRIC_bob.pl           the BV-BRC export used to produce the inputs
```

## What to replace in `build_collections.py`

Everything above `def compile_rules` is Rhabdoviridae. Everything below it is
general and should be kept as-is.

| table | what it is | replace |
|---|---|---|
| `ALPHA_GENERA`, `BETA_GENERA`, `GENUS_MODULE` | which genus routes to which module | yes |
| `CORE` | family-wide regexes for N, P, M, G, L | mostly — the shape is reusable |
| `GENUS_RULES` | genus-specific accessory proteins | yes |
| `EXPECTED` | length range per feature, used to split outliers | yes |
| `MISBINNED` | per-sequence overrides for annotations that are simply wrong | yes |
| `NOT_MODELLED` | genera deliberately excluded, with the reason | yes |
| `NTERM_TRIM` | start-codon misassignments in the source | yes |
| `UNCHAR_RULE` | what counts as an unnamed protein | probably reusable |
| `UNCHAR_L_RESCUE` | length above which an unnamed protein must be the polymerase | retune |
| `UNC_PIDENT`, `RESCUE_PIDENT` | 0.80 — the line at which two proteins are the same protein | keep |

The general machinery below, which you should not need to touch:
`homology_rescue`, `_drop_named_duplicates`, `split_unchar`,
`reroute_length_outliers`, `merge_extracted`, `classify`, `compile_rules`.

## What to replace in `build_rhabdoviridae_json.py`

The `TAXA[...]` blocks and the `ANNO_*` / `PMID_*` constants. Keep `feat`,
`unchar`, `seg`, `bounds`, `bit_for`, `derive_bounds`, `enforce_designations`,
`declare_unchar_groups`, `canonical` and `sort_deep`.

## Read the comments

Around 1,000 lines of these two files are comments recording why a rule exists,
usually with the measurement that prompted it — why `outer coat protein` is the
glycoprotein and not the nucleocapsid, why `M1` is the phosphoprotein, why
`polyprotein` is safe as an L synonym but `P6` must never be a regex at all.
Those are the decisions that are expensive to rediscover, and most of them
generalise even where the strings do not.
