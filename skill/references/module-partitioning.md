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

**Prefer the single taxon when one exists.** Reach for the concatenation only
after looking for a subfamily or genus that already covers the group, because
finding one usually means the partition landed where the taxonomy already
is — which is a good sign about the partition.

