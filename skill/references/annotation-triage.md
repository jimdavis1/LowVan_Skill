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

**Adopt at >= 80% identity over >= 60% of the query. Below that, leave the
protein out.** Not adopted, and not swept into the uncharacterized bag either.
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
