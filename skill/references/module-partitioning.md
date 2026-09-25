# Deciding what the modules are

Settle this before writing a single classification rule. It is the one decision
that is expensive to reverse, because collections, alignments, PSSMs, JSON keys
and installed directory names all carry the module name.

## The rule

**A module is a group of taxa that share a genome organisation.** Not a clade.

The annotator picks a module by BLASTn against rep contigs, then assumes the
incoming genome has *the same set of proteins* as that module declares. So the
grouping has to be by gene layout, because that assumption is what the module
encodes.

Consequences:

- Genera from **different subfamilies** belong in one module if they carry the
  same core gene set. Rhabdoviridae's Alpharhabdovirinae module holds ~38 genera
  spanning two subfamilies, all of them N-P-M-G-L single-segment.
- A genus with an **extra segment** gets its own module even if it is closely
  related. Dichorhavirus is two segments and Trirhavirus is three, so neither
  can share a module with the single-segment plant rhabdoviruses.
- A genus with a **different accessory set** gets its own module. Novirhabdovirus
  carries NV, which nothing else does.

## Why collapsing is safe

The obvious worry about putting 38 genera in one module is that a genus-specific
accessory PSSM will fire on a sibling genus. It does not, because `bit_cutoff`
does that work, not the routing. A Tibrovirus U1 profile sitting in the same
directory as a Lyssavirus genome simply scores below its cutoff.

What collapsing buys you: one rep contig per genus routes to one shared set of
core PSSMs, so N/P/M/G/L are built once from the deepest possible collection
instead of 38 times from thin ones.

## Cutoffs and judgement

A working heuristic: a genus with **fewer than ~5 genomes** in the dump does not
justify its own module. Fold it into the nearest module with the same layout and
accept that its accessories will not be called.

Check the genome organisation against ICTV rather than against the dump. The
dump tells you what was submitted; ICTV tells you what the genome is. In this
project Varicosavirus was declared as a module, built, and then dropped after
ICTV showed a set of accessory genes that could not be supported from the
available sequences — better to have no module than one that misses half a
genome.

## Declaring a module you cannot build

It is legitimate to declare a module in the JSON with no PSSMs behind it, as a
marker that the taxon is known and unbuilt. Two Rhabdoviridae modules are in
that state: the BV-BRC dump holds no protein sequences at all for
Alphagymnorhavirus or Trirhavirus.

If you do this, `install_module.py` will skip the module and say why, and the
registry artifact will show it as `no-collection`. That is the intended
behaviour — the declaration documents the gap instead of hiding it. What you
must not do is declare *features* inside a working module with nothing behind
them and leave that unrecorded, because at runtime they are silently never
called.

## Declared is a promise; retire what you cannot keep

A feature in the JSON is a claim that the module can call that protein. When
there is no profile behind it the claim is false, and false in a way nothing
downstream detects: the annotator looks for a PSSM directory, finds none, and
moves on silently. The count only surfaces if you go looking for it.

Retire a feature when the sequences will not support a profile. Not enough data
is not a reason to ship a call. In this build that removed 33 features at once
-- 21 with no sequences at all, 12 with a handful that never clustered -- and
one whole module.

Two of those had real literature behind them, which is the case worth thinking
about. `ANUCLEO_P6` is a described RNA silencing suppressor with a PMID, and its
two sequences are genuinely homologous: 44% identity over 54 residues, e-value
2e-09, sharing a C-terminal `[VL]YELIDEC` motif, and hitting nothing else in the
subfamily. The declaration was right. It still goes, because **a citation does
not make a profile.** Two sequences cannot train one at any threshold.

So retire on the data, not on the strength of the naming. Comment the
declaration out rather than deleting it, with the sequence count and what would
unblock it, so reviving it later is uncommenting rather than re-deriving.

### Retire the designation with the feature

If you keep a list of designations considered literature-backed, prune it in the
same pass. Nine retired keys were still listed as backed here, and because that
list is consulted when deciding whether an annotation may carry a positional
label, it went on shielding *surviving* features: `BCYTO_P4`, `BCYTO_P5`,
`SAWGR_GY`, `CURIO_PMIP`, `BNUCLEO_U` and `ACYTO_VIROPORIN` were all still
annotated "... P4 protein", "... Gy protein", "Class Ia viroporin P9 protein"
on the authority of entries for features that no longer existed. Clearing the
dead entries let the designation check do its job and collapse all six to the
plain uncharacterized form.

A designation is only backed while the feature it names exists.

## Recording the decision

Write the partitioning down with reasons, in the working directory, before you
build. Something a stranger can audit:

```
Module               Genera                        Layout            Why
Alpharhabdovirinae   Lyssavirus, Vesiculovirus...  N-P-M-G-L         shared core, 38 genera
Betarhabdovirinae    Alphacytorhabdovirus...       N-P-P3-M-G-L      plant, single segment
Dichorhavirus        Dichorhavirus                 2 segments        layout differs
Novirhabdovirus      Novirhabdovirus               N-P-M-G-NV-L      carries NV
```

Then map each module to its rep contigs, one per genus where you can get an
exemplar, and keep the accessions with provenance. Never guess an accession by
pattern-matching `NC_` numbers — look each one up and record what confirmed it.

## Name the module after a real taxon

The module name is not a label. It is the key in `Viral_PSSM.json`, the
`Viral-PSSMs/<Module>.pssms/` and `PSSM-Alignments/<Module>/` directory names,
the `Splice-Variants/<Module>/` lookup that `get_splice_variant_features.pl`
does, and the `<Module>.<n>.dna` rep-contig filenames. It is also written into
every output GTO as **`viral_family`**, so anything downstream that reads that
field as taxonomy gets whatever was invented here.

**Partition on genome organisation, then find the taxon that names the
result.** The two are usually the same, because gene layout is what the
taxonomy was built from. Betaflexiviridae splits cleanly into a
single-movement-protein half and a triple-gene-block half — worth 16 points of
routing — and those halves are exactly the ICTV subfamilies **Trivirinae** and
**Quinvirinae**. Shipping them as `Betaflexiviridae_MP` and
`Betaflexiviridae_TGB` would have invented two names for groups that already
had them.

When a member has to be handled separately, the parent keeps the parent's
name and the exception gets its own real name:

| module | split out | both valid |
|---|---|---|
| `Orthopneumovirus` | `Orthopneumovirus_muris` | ICTV binomial |
| `Alpharhabdovirinae` | `Merhavirus` | genus |
| `Trivirinae` | `Capillovirus` | genus; its CP is inside the ORF1 polyprotein |

**If no single taxon covers the group, name every member.** Botrexvirus,
Platypuvirus and Sclerodarnavirus share no subfamily, so that module is
`Botrexvirus_Platypuvirus_Sclerodarnavirus`. It is long and it does not
resolve anywhere, and both are fine: every component is a real genus and the
membership is on the label. What is not fine is `Alphaflexiviridae_noTGB`,
where `noTGB` is a property of the group rather than a name for it and you
cannot tell from the string what is inside.

### The test

Every component of the name is a real taxon, and a reader can tell what the
module contains without opening it.

| | |
|---|---|
| `Trivirinae` | best: one taxon covers the group exactly |
| `Merhavirus`, `Orthopneumovirus_muris` | a real taxon split out of its parent |
| `Botrexvirus_Platypuvirus_Sclerodarnavirus` | fine: no taxon covers them, so all three are named |
| `Betaflexiviridae_MP`, `Alphaflexiviridae_noTGB` | wrong: `MP` and `noTGB` are descriptions, not names |

Resolving in NCBI taxonomy is a useful check for a single-taxon name and
nothing more — a concatenation will not resolve and does not need to. The
requirement is that nothing in the name was invented.

**Reach for a real taxon first. The concatenation is the fallback.**

Before naming a module anything else, look for the subfamily, genus or
species that already covers the group. Look properly: `Betaflexiviridae_MP`
shipped because the split was made on gene layout and never checked against
the taxonomy, and **Trivirinae** had covered those exact eight genera the
whole time. One NCBI lookup would have found it.

Finding a real taxon is also evidence the partition is right. Gene layout is
what the taxonomy was built from, so a split that lands on an existing
subfamily has landed where the biology already agrees. A split that lands on
nothing is worth a second look before it is worth a compound name — the
group may be wrong, not merely unnamed.

Only once that search comes up empty does the concatenation apply.

## Splitting on measurement, when the organisation rule says one module

The rule above is that a module is a group of taxa sharing a genome
organisation. Sometimes every genus in a family shares the organisation and the
module still should not be one. **Gate the split before building and let the
measurement decide.**

The mechanism is budget competition. 25 references are chosen greedily,
largest cluster first, so a divergent genus absorbs references a tight one
needed and neither ends up well covered.

### Quinvirinae, 24 September 2026

Six genera, all with replicase / TGB1-3 / CP / NABP. Lumped they route
**81.0%**, and the loss is not confined to the divergent genus:

| genus | total | uncovered lumped | lumped | split |
|---|---|---|---|---|
| Foveavirus | 794 | 43 | 94.6% | **100%** |
| Carlavirus | 714 | 246 | 65.5% | **76.3%** |
| Robigovirus | 121 | 22 | 81.8% | **100%** |
| Ravavirus, Sustrivirus | 4 | 4 | **0%** | covered |

13 of the 25 references go to Carlavirus, which has 102 clusters over 714
genomes and still reaches 65.5%; Foveavirus needs five references for 92%, gets
seven, and falls from 100%. **Splitting improved every genus, the divergent one
included.**

The collections agreed independently: coat protein is 299 aa in Carlavirus and
379 in Quinvirinae — a shared window spans both and admits junk — and NABP is a
Carlavirus feature that a lumped module would declare on 1,983 genomes that do
not have it.

### Closteroviridae, 24 September 2026

Six genera sharing ORF1a / ORF1b / p6 / HSP70h / HSP90h / CPm / CP.

| grouping | 25 refs |
|---|---|
| lumped | **77.0%** |
| two-way | 88.6% / 89.5% |
| **three-way** | **95.8% / 89.5% / 99.4%** |

Three beats two for every group; a fourth module adds nothing over the
three-way. Test the groupings rather than assuming one split point.

### Tymoviridae — the precedent

Tabled at **74.6%** lumped. Gated per genus: Marafivirus 100%, Tymovirus 91.2%,
Maculavirus 100%. Same genomes, same budget. **The lumped figure was measuring
the module boundary, not the taxon.**

### How to run the test

```bash
#  monopartite: records at >=6,000 nt, no merging
python3 scripts/repcontig_budget.py --meta <class>.full.tsv --any-family \
    --family X --genera "GenusA" --min-len 6000 --budget 25 --label GenusA

#  multipartite: merge records sharing a genome_name FIRST
python3 scripts/repcontig_budget.py ... --merge-by-name
```

**Do not guess a per-taxon length floor.** Two attempts failed differently: a
floor from the median record length admits fragments and reported Quinvirinae
at 89.6% against a true 81.1%; a floor from a p90 of merged lengths broke the
other way and gated Tymovirus on 5 genomes of 80. One fixed 6 kb criterion,
merged only where the genome is genuinely segmented, reproduces independent
measurements exactly.

**And merge only what is segmented.** For a monopartite virus, records sharing
a genome_name are separate submissions of the same isolate, and concatenating
them manufactures coverage — Closteroviridae reads 97.5% merged against 77.0%
unmerged.
