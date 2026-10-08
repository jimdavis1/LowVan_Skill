# Triaging annotation strings into feature collections

This step decides which sequences train which profile. Everything downstream
inherits its mistakes, and the mistakes are quiet — a misbinned cluster produces
a perfectly healthy-looking alignment and PSSM that calls the wrong protein.

## Start with the histogram

```bash
paste <(cut -f2 <T>.uniq.id_ann) <(cut -f1 <T>.uniq.id_ann) \
  | sort | uniq -c | sort -rn
```

Better, cross-tabulate string × genus, because the same string means different
things in different genera. Do this before writing any rule.

## The U-number trap

**Positional labels are not homology groups.** `U1`, `U2`, `ORF3`, `VP2`,
`protein 4`, `P6` name a slot in a genome, not a protein family. Two viruses'
`U1` are usually unrelated, and binning them together produces a profile that
matches nothing.

The test is cheap and mandatory: before building a feature from a positional
label, blast its members against each other. If they are not homologous, they
are not a feature.

Worked example from this project: `SRIPU_GX` (2 sequences), `SRIPU_MX` (1) and
`SRIPU_U1` (1) were checked pairwise and are not homologous to one another —
three positional labels covering three different proteins. Forcing a profile
from any of them would generate false positives, so all three stayed unbuilt.

The same trap in reverse: Stangrhavirus labels its first three genes
`VP1/VP2/VP3` while a sister genus labels the same three `ORF1/ORF2/ORF3`, and
gene order plus a homology hit to the N collection is what identifies them as
N/P/M. Positional labels can be *resolved*, but only with evidence — gene order,
length, and a homology hit — never by the number in the name.

## Every rule needs a genus column

The same product string maps to different features depending on genus. Real
example, and a bug this project shipped before finding it:

| String | Genus | Feature |
|---|---|---|
| `coat protein` | Varicosavirus, Lyssavirus | N — it is the nucleocapsid |
| `outer coat protein` | Sigmavirus | **G** — "outer coat" is the envelope |

`outer coat protein` was originally in the N rule, which sits ahead of G in the
rule order, so 176 Sigmavirus glycoproteins were binned as nucleocapsid. The
cross-feature QC found it when a G cluster's master turned out to be 99.8%
identical to a member of the N collection.

So: rules are `(feature_key, genus_or_None, regex)`, matched **in order**, and
order matters. Put genus-specific rules ahead of general ones.

## Match twice: raw, then normalised

BV-BRC punctuates the same product string inconsistently. Match the annotation
exactly as supplied, then match a normalised form with hyphens and underscores
turned into spaces and whitespace collapsed:

```python
norm = re.sub(r"[-_]", " ", ann)
norm = re.sub(r"\s+", " ", norm).strip()
if norm != ann:
    for key, genus, prx in rules:
        if prx.search(norm):
            return key
```

Strictly additive — anything that matched before still matches on pass one. It
exists because `L-protein` missed `^L (protein|polymerase)` on a single hyphen
and sat unbinned, which is how Varicosavirus lost a third of its L sequences.

## Discard anything the source annotates as partial

**A record whose annotation says it is partial or incomplete is DISCARDED.**
Not binned, not offered to the homology rescue, not parked in
`<FEAT>.outliers.fasta`. The submitter is telling you the sequence is not the
whole protein, and a profile has no business being trained on it.

This is a different test from the length gate below, and it runs first. The
length gate asks *is this sequence the wrong size*; this asks *does the source
say it is a fragment*. A record can fail either independently.

The qualifier vocabulary is wider than it looks. Measured across one
Sedoreoviridae dump:

```
, fragment     43,026 features
, N-terminal   17,743
, C-terminal    6,135
truncated X        25     <- a PREFIX, not a comma suffix
```

so match both shapes:

```python
PARTIAL = re.compile(r",\s*(fragment|partial|incomplete|truncated|[NC]-?terminal)\s*$"
                     r"|^\s*truncated\b", re.I)
```

The 25 `truncated VP7` / `truncated NSP3` records had been binned as clean
because only the comma form was being matched.

### N-terminal and C-terminal are not always partiality markers

A trailing bare positional qualifier is the submitter saying the record is only
a piece. The same words used *adjectivally inside a full product name* are
naming a genuine cleaved subunit, and those are kept:

```
putative capsid glycoprotein VP7, N-terminal          -> DISCARD (a piece of VP7)
Mature N-terminal envelope glycoprotein Gn            -> KEEP    (a real subunit)
Mature hemagglutinin C-terminal membrane fusion subunit -> KEEP  (HA2)
Mature N-terminal spike protein receptor binding domain S1 -> KEEP
```

32 strings in S1 Table are of the second kind -- Gn/Gc across five bunyavirus
families and both filoviruses, HA1/HA2, HEF1/HEF2, spike S1/S2, fusion F1/F2.
Discarding those would delete a cleaved product from every module that declares
one.

The anchoring is what separates them: the qualifier must be **comma-introduced
and at the end of the string**, with nothing after it. `check_annotations.py`
asserts that no string in the published vocabulary is classified partial, so a
future loosening of this regex fails the QC run instead of silently removing
subunits.

**`truncated` is always discarded**, in either shape -- as a comma suffix or as
the `truncated VP7` prefix. There is no homology test that rescues it, because
the objection is not to the identity but to the completeness.

**`putative` is a different thing, and it is kept -- but it has to be earned.**
It is not a partiality claim; it is the submitter saying they are unsure *which*
protein this is. So the string is not taken at face value: the record goes to
UNASSIGNED and the homology rescue adopts it only at >=80% identity over >=60%
of the query against a named collection. `probable`, `possible` and `predicted`
are treated the same way.

Binning `putative X` by name instead would let one uncertain submitter's guess
into the training set for X, where it trains the profile that then propagates
the same guess. Routing it through homology means the collection decides whether
the guess was right.

### Why parking them in outliers is not good enough

Rotavirus VP7 is the worked example. It is full length at 326 aa in Rotavirus A
and at 244-251 in RVB/RVG/RVJ, so the module window has to reach down to ~195 to
keep the short species -- and a *truncated RVA* VP7 at 276 then passes it.

    RVA VP7, clean label     5,461 records, median 326   <- true full length
    RVA VP7, partial label  10,105 records, median 276   <- truncated

The truncated records outnumbered the full-length ones two to one. No single
length window separates them, because the truncated RVA sequences are the same
size as legitimate RVB ones. Discarding on the annotation does separate them:
RVA VP7 goes straight back to a 326 median.

The homology rescue makes this worse if it is allowed to see them, because they
*are* genuine VP7 by homology -- it adopted 36,800 such records at >=80%
identity and handed them all back. They have to be removed before the rescue,
not after.

Counts still go to `NOT_MODELLED.tsv`, and an audit file naming every discarded
string keeps the decision reviewable. What is not kept is the sequence.

## Length gates catch some errors, not this one

Keep an `EXPECTED = {feature: (min_aa, max_aa)}` table and split anything outside
the range into `<FEAT>.outliers.fasta` rather than discarding it. Most outliers
are partial single-gene submissions; some are real misannotations — a 2,128 aa
"glycoprotein" is an L protein.

Be aware of the limit. A 524 aa protein labelled "nucleoprotein" sits inside any
sane nucleocapsid range (a genuine Betaricinrhavirus N is 557 aa), so the length
gate passes it and the cross-feature QC is what catches it. Do not tighten the
ranges to compensate; you will lose real proteins. Run
`scripts/qc_cross_feature.py` instead.

## Length also rescues, not just rejects

The `EXPECTED` table above is a filter, but the same measurement runs the other
way. Every rule in this file reads annotation text, so a protein the source
labelled only "hypothetical protein" carries no signal at all and lands in the
grab-bag no matter what it is. Length still knows.

Nine Alpharhabdovirinae `UNCHAR` members were 1452-2111 aa. They are L
polymerases, and **none of those nine genomes had an L binned at all** — the
collection was simply missing them, and nothing in the annotation text could
have said so. The tell was a derived bound: `UNCHAR` came out at 47-1947 aa,
which is not a protein family, it is a bag with a polymerase in it.

The rescue is safe when the sizes do not overlap. Excluding L, the longest
protein in any Rhabdoviridae collection is 814 aa; L runs 1804-2295. Nothing
sits between, so a floor at 1400 cannot take anything else:

```python
UNCHAR_L_RESCUE = 1400
...
if key == "UNCHAR" and len(s) >= UNCHAR_L_RESCUE:
    key = "L"
```

Route to the feature key, never straight to the fasta. Routing means `EXPECTED`
still applies afterwards, so full-length rescues join the collection and the
partials (1452-1637 aa here) go to `L.outliers.fasta` where partials belong. A
rescued sequence should face exactly the same scrutiny as one that was named
correctly — appending it directly would smuggle four truncated polymerases into
the profile that the length gate exists to keep out.

Check for this whenever a grab-bag's derived length range spans more than about
threefold. Look for the module's largest feature first: the polymerase is the
one that hides in plain sight, because nobody expects the unnamed protein to be
the biggest one in the genome.

## Homology rescue, and the line at 80%

Text rules miss two kinds of protein: the ones whose annotation is a typo or a
synonym you did not think of, and the ones whose annotation is not a name at
all. Blast every protein no rule matched against the named collections **of its
own module** and the difference is stark — in this dump 485 of them were a
named feature at 80% identity or better.

The recoveries are humbling, because each is a real protein lost to a string:

```
"nulcleoprotein"   99.5% -> N     a typo
"non-viron"       100.0% -> NV    a spelling
"P(M1) protein"    95.2% -> P     a synonym nobody listed
(blank)            90.9% -> N     no annotation at all
```

**Adopt at >= 80% identity over >= 85% of the query. Below that, leave the
protein out.**

> The coverage floor was **0.60 until 1 October 2026** and is now **0.85**, on
> the curator's instruction after the Alsuviricetes backfill: *"bring the
> coverage up to 85% if we ever do the rescue again, i don't want to slurp back
> in all those short proteins."* Every module rescued before that date was built
> at 0.60 and its audit says so. `filter_unchar.py --cov` and
> `reroute_outliers.py --min-qcov` were raised with it and must stay equal.
>
> Know what the floor measures. It is coverage of the **query**, so it rejects
> a partial *alignment*, not a short *protein* — a truncated sequence aligns
> over nearly all of itself and passes. The length window is what bounds
> truncation after adoption. Measured at 0.60 over six modules: of 2,196
> adoptions, 20 (0.9%) fell below their collection's previous minimum length.
> To exclude short proteins outright, use a floor at the collection's existing
> minimum; this flag will not do it. Not adopted, and not swept into the uncharacterized bag either.
Below 80% it is its own lineage-specific protein, and absent a literature-based
function saying what that protein is, there is nothing to build. Filing it under
"uncharacterized" only dilutes that collection with proteins you declined to
identify.

Two things make leaving them out the right call rather than a loss. A collection
is training data, not an inventory — the annotator blasts each finished profile
back against the incoming genome, so a protein omitted here can still be called
by the profile its relatives built. And the alternative asserts a homology the
number does not support: Rice yellow stunt gene 3 is 69% identical over 99% of
its length to a movement protein, which makes it movement-protein-*like* and not
a movement protein.

Whatever threshold you choose, **the UNCHAR filter must use the same one**. If
adoption says 80% means "same protein" and `filter_unchar.py` says 90%, a
protein 85% identical to a named feature is adopted into that feature *and* left
in the grab-bag, and both profiles fire on one locus.

### What must never be gathered by regex

`P6`, `ORF3 protein`, `hypothetical protein`, `gp3`, `protein 3`. These are
positions in a genome or a shrug, not proteins. A rule keyed on one of them
manufactures a feature out of sequences with nothing in common — the U-number
trap with a different label. If such proteins are to be grouped at all, only
homology may do it.

The distinction is not about how vague a word looks, it is about whether the
word names a *protein*. `replicase`, `viral polymerase`, `polyprotein` and
`large polymerase protein` are all broad, and all fine: nothing in a rhabdovirus
except L runs 1800-2300 aa, so `EXPECTED` catches any mistake and routes it to
`L.outliers.fasta`. `P6` names the sixth reading frame and says nothing about
what is in it.

Prefer the submitter's own word to a homology score wherever the word names a
protein. Fifty polymerases in this dump were being recovered by blast when they
could have been bound on annotation, simply because the L rule lacked
`rna-directed` (as against `rna-dependent`) and four other synonyms. Widening
the regex is better evidence and cheaper than blasting.

## When no rule can fix it

Some source annotations are simply wrong about which protein the sequence is,
and the offending string is a legitimate name for the feature it landed in. No
regex helps. Keep a per-feature override table alongside the rules:

```python
MISBINNED = {
    # feature_id -> (correct key or None to drop, the evidence that settled it)
    "fig|57482.650.peg.2": ("P", "298 aa 'matrix' in the P slot; real M is peg.3 at 202 aa"),
}
```

Apply it right after `classify()` and print every hit, plus a warning for
entries that matched nothing so stale ids surface on the next rebuild. Listing
them individually is the point: a rebuild must not silently reintroduce them,
and each one should carry the evidence rather than just an assertion.

Two real cases from Rhabdoviridae, both found by `qc_cross_feature.py`:

- a 298 aa protein annotated `matrix` in a Lyssavirus whose real M is 202 aa in
  the next slot — the source data mislabelled a phosphoprotein
- nine rabies glycoproteins annotated `nucleoprotein`, which had to be dropped
  at the cluster level rather than the rule level because the string is right
  for N in every other genome

Prefer a rule fix when the string is systematically wrong (`outer coat protein`
always means G). Use the override only when the same string is correct
elsewhere.

## Nothing may vanish

Emit three side files so every input sequence is accounted for:

- `synonyms.tsv` — `count \t annotation \t genus \t bins_to`, one row per
  string×genus. This is the input to the string-collapse artifact and the
  auditable record of every binning decision.
- `UNASSIGNED_TRACKING.tsv` — strings that matched no rule, with genus and
  count, so a reader can see what was left out and decide whether it matters.
- `NOT_MODELLED.tsv` — features deliberately excluded, with the reason.

A string that matches no rule is a finding, not a leftover. Read the unassigned
list before building; that is where the missing accessories are.

## The uncharacterised grab-bag

Real proteins with no name are common — 152 in Hapavirus alone in this dump.
Collect them under an `UNCHAR` key so they are not lost, and expect the profile
to be weak: Alpharhabdovirinae `UNCHAR` recovers 4% of its own collection,
because a single profile cannot cover a bag of unrelated proteins. That number
is honest, not a failure to fix. Do not tune it; either split the bag into real
features when the literature names them, or leave it as a documented floor.

## A mass is not a homology group, any more than a position is

The U-number trap above covers positional labels — `U1`, `ORF3`, `VP2`. The
same applies to **masses**, and closterovirids are named almost entirely that
way: `p13`, `p20`, `p23`, `19 kDa protein`, `5 kDa protein`.

Two consequences, both measured on the Closteroviridae split:

**String rules miss core features.** First-pass binning reached only 37–46% of
occurrences because the source calls things by mass:

| module | what the rule expected | what the source says |
|---|---|---|
| Ampelovirus | polyprotein | `methyltransferase/helicase` |
| Ampelovirus | HSP90 | `p55`, `60 kDa protein`, `Hsp90-like` |
| Closterovirus | HSP70h / HSP90h | `p65` / `p63` |
| Closterovirus | CP / CPm | `p25` / `p27` |
| Ampelovirus | CP | also `divergent coat protein`, `coat protein duplicate` |

Clustering the unbinned pool found all of them and lifted binning to
**83.0% / 74.3%**.

**Similar masses are not the same protein.** Ampelovirus `19.6 kDa protein`
and `19.7 kDa protein` cluster **separately** — p20A at 177 aa and p20B at 179.
Any rule keyed on mass or length merges two distinct proteins.

The same held for Quinvirinae's accessories: 471 positionally-labelled
sequences clustered into three groups sharing no homology, of which p14
(median 119 aa) and ORF2A (median 123) are inseparable by length, and 295 of
the records are labelled only `hypothetical protein`.

### The method

1. Bin what the strings clearly name.
2. Cluster the **unbinned pool** at ~30% identity, 60% coverage.
3. Read each cluster's dominant label as a *description* of the group, not as
   its definition. Groups above ~25 members are worth declaring.
4. Revise the rules from what the clustering found, then rebuild.
5. For records carrying no informative label at all, seed each collection from
   the explicitly-labelled members and assign the rest by blastp — and **leave
   the ones that match nothing unassigned rather than guessing**. On
   Quinvirinae that placed 298 of 309 and left 11.

## Set length windows per module, and check them per genus

A relative floor — a fraction of the feature's own median — fails on a bimodal
or genus-split distribution, because the median sits in the larger mode.

Quinvirinae coat protein, measured by genus:

| genus | n | p10 | median |
|---|---|---|---|
| Foveavirus | 858 | 256 | **395** — bimodal, 302 in a short class |
| Robigovirus | 132 | 267 | 268 |
| Banmivirus | 49 | 238 | **238** |

Pooled median 376 puts a 70% floor at 263, which deleted **46 of Banmivirus's
49**, Sustrivirus's only one, and all 302 short Foveavirus. The short and long
Foveavirus classes are homologous at ~41% identity over 74–80% coverage, so
both are real. Explicit bounds of 200–460: **679 sequences → 987**.

This is the Betaflexiviridae/Citrivirus lesson restated — and note it bit again
even with that lesson written down, because the failure was *within* a genus,
not across genera.


## Use the engine, not a copy of the last taxon's script

`scripts/collection_engine.py` holds the binning logic -- rules, length
windows, ambiguity, cleavage chemistry, every tracking file. A taxon's
`build_collections.py` should declare only what is specific to it and call
`run()`. It was extracted after the fourth taxon reproduced the same 180
lines with a different `RULES` table, and the copies had already started
drifting.

```python
import os, sys
sys.path.insert(0, os.environ["LOWVAN_KIT"] + "/skill/scripts")
from collection_engine import run

run(here=..., T="Crini", MODULE="Crinivirus",
    RULES=[...], NOT_MODELLED=[...], WINDOW={...},
    PRE_FLOOR={...}, JUNCTION={...}, ACCESSORY={...})
```

## A pooled median cannot arbitrate between two modes

`RUNT`/`GIANT` (0.7x-1.3x the feature's median) assume the feature is
unimodal. When it is not, the median lands in the larger mode and the smaller
mode is deleted wholesale as "runt" or "overlong".

Capillovirus MP has a 320 aa class and a 463 aa class. They are 42% identical
over 212 residues (E=6e-61) -- divergent but unambiguously the same ORF2
homology group, and both are real. The pooled median is 320, so the default
window is 224-416 and **all 119 sequences of the 463 class are discarded**,
with every number in the report self-consistent.

This is the Citrivirus CP failure arriving through the length filter instead
of through the string rules, and it needs the same answer: an explicit window.

```python
WINDOW={"MP": (260, 480)}     # replaces the median-relative test entirely
```

**Plot the length histogram of every feature before accepting its window.**
Two modes with a gap between them means either two features or one feature
needing an explicit window; it never means take the median.

## A mass is not a homology group -- so do not let a mass form one

Already recorded for U-numbers and positional labels. Closterovirid
accessories make the point quantitatively. Measured on the Crinivirus dump:

    "p22"        FOUR non-homologous groups, at 188, 191, 192 and 193 aa
    "p26/27/28"  one group holding all three names, plus two further
                 separate groups that also call themselves p26/p27/p28
    "p59"/"p60"  ONE group of 78 sequences -- and it is the HSP90 homolog,
                 which Closterovirus, Ampelovirus and Velarivirus all name
                 HSP90 rather than by its mass
    "p6/7/8/9"   at least seven groups, fragmenting by size

Binning any of those on the string invents proteins that do not exist and
splits ones that do. The mass is a *name*, and a name is only safe once the
group exists.

So: **discover accessory features by clustering, then name them.**

```python
ACCESSORY={"pool": r"^p\d{1,3}(\.\d)?( protein)?$|^\d{1,3} kda protein$|"
                   r"suppressor of rna silencing",
           "min_members": 8,      # enough to build a profile at -m 2
           "min_med_len": 45,     # enough to give it a usable bit_cutoff
           "min_seq_id": 0.25, "cov": 0.5}
```

The engine clusters the pool, keeps groups that clear both floors, names each
after its dominant label, and suffixes A/B/C where one name covers several
groups -- which is the Ampelovirus `P20A`/`P20B` precedent made general.
`ACCESSORY_CLUSTERS.tsv` records what became what.

Two floors, both load-bearing. `min_members` asks whether a profile can be
built at all. `min_med_len` asks whether it can be given a cutoff: the
discovered Crinivirus "p4" group is 33 aa, and a 33-residue profile cannot
both fire on its own members and stay silent across a 9 kb genome. Ship it and
it can only be wrong in one of two directions.

### The free validation this buys you

For a multipartite taxon, check which segment each discovered group sits on.
Anchor two record sets -- the ones carrying POLY/RDRP are RNA1, the ones
carrying HSP70/CPM are RNA2 -- and ask where each group's members fall. On
Crinivirus the two anchor sets did not overlap at all, and **every one of the
23 features landed entirely on one segment or the other**. A cluster that
straddled both would have been an artefact; none did. It costs one query and
tells you the clustering found real genes.

## Re-running the binning leaves the previous rule set on disk

`build_collections.py` writes `collections/<Module>/<KEY>.fasta` per feature.
It does not remove files it no longer produces, and **every downstream step
globs that directory**. Rename a feature and both names persist; change an
accessory suffix and the old cluster stays. `SUPPRESSOROF` and `SUPPRESSOR`
sat side by side, and a PSSM would have been built for both.

The engine now deletes any `*.fasta` the current run did not produce and says
which. That interacts with the rescue, which appends: the order is always
**build collections, then rescue, then compare every count against the
previous build** -- never one without the other.
