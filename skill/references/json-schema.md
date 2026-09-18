# The Viral_PSSM.json module block

One top-level key per module. The key is the module name and it must match the
PSSM directory (`Viral-PSSMs/<key>.pssms/`) and the rep-contig filenames
(`Viral-Rep-Contigs/<key>.<n>.dna`) exactly.

```json
"Alpharhabdovirinae": {
  "close_genomes": {
    "Alpharhabdovirinae.1.dna": {
      "genome_ids": "NC_001542.1",
      "genome_name": "Rabies lyssavirus Pasteur virus"
    },
    "Dichorhavirus.1.dna": {
      "genome_ids": "NC_009608.1,NC_009609.1",
      "genome_name": "Orchid fleck virus"
    }
  },
  "segments": {
    "Single RNA Segment": { "max_len": 16000, "min_len": 10000,
                            "replicon_geometry": "linear" }
  },
  "features": {
    "N": {
      "anno": "Nucleocapsid protein",
      "gene_symbol": "N",
      "feature_type": "CDS",
      "segment": "Single RNA Segment",
      "bit_cutoff": 100,
      "coverage_cutoff": 0.65,
      "upstream_ext": 1,
      "downstream_ext": 1,
      "copy_num": 1,
      "kmers": 1,
      "min_len": 330,
      "max_len": 560,
      "PMID": ["1404364"]
    },
    "NV": {
      "anno": "Non-virion protein",
      "PMID": ["9010293"],
      "PMID_claude_generated": ["4038520", "8683214"]
    }
  }
}
```

## Fields

| Field | Meaning |
|---|---|
| `anno` | the annotation string emitted; a functional description. Must follow the controlled vocabulary. |
| `gene_symbol` | the community's short name (`L`, `N`, `SH`). Must match the vocabulary's Symbol column where the annotation already exists there; do not namespace it per genus. |
| `feature_type` | `CDS`, `mat_peptide`, or `RNA` |
| `segment` | which segment; used by the quality tool |
| `bit_cutoff` | tBLASTn bitscore floor for calling this feature |
| `coverage_cutoff` | subject coverage floor, typically 0.65 |
| `upstream_ext` | may scan upstream for a Met start |
| `downstream_ext` | may scan downstream for a stop codon |
| `copy_num` | expected copies; quality tool flags departures |
| `min_len` / `max_len` | expected feature length in **amino acids** |
| `PMID` | DLITs — papers defining the function or sequence. Bare numeric PubMed ids only, never DOIs. |
| `PMID_claude_generated` | citations a language model proposed rather than a curator choosing. **Disjoint from `PMID`** — never list an id in both. See below. |
| `internal_stop` | `1` if the match may legitimately contain an in-frame stop (readthrough). Omit otherwise — see below |
| `special` | `transcript_edit` or `splice`; called by an external program, no PSSM expected |
| `non_pssm_partner` | places a location-based feature relative to another |

`segments` entries carry `min_len`/`max_len` in **nucleotides** for the whole
replicon, plus `replicon_geometry`.

### `internal_stop`, and why there are no zeros

Set it to `1` on a feature whose profile match may legitimately span an in-frame
stop — a readthrough codon such as the alphavirus opal. Without it, tblastn's
match is cropped at the stop and the protein comes out short by a constant number
of residues across most of the taxon.

**Omit the flag rather than writing `internal_stop: 0`.** Then
`grep internal_stop` returns exactly the features that permit one, which is the
question an auditor actually asks; writing the zeros makes that grep return
everything and forces a read of values instead.

The annotator caps the number of tolerated stops with `-mis N` (default 1): one
readthrough codon is expected, a run of stops is a broken genome and is cropped.
The flag is necessary but not sufficient — the feature's profiles must also be
built from alignments that span the stop, or they have never seen what follows
it. Full treatment in `references/special-features.md`.

## Setting the numbers

**`min_len` / `max_len`** from the collection's own length distribution, not
from a textbook. Take the observed range and widen it a little; too tight and
the quality tool flags healthy genomes.

**`bit_cutoff`** is what keeps a genus-specific accessory from firing on a
sibling genus inside a collapsed module. Raise it for short accessories and for
anything whose profile is close to another feature's. If a feature is calling
the wrong locus, the cutoff is not the fix — run `qc_cross_feature.py`, because
the cluster is probably misbinned.

Set it from two measurements, not by feel. It must sit **above the noise
ceiling** — the best score the profile reaches on genomes that do not contain
the feature — and **below the weakest true signal**. Short profiles are where
this matters most, because the two can converge:

```
Alpharhabdovirinae G_SP, ten 19-residue profiles, seven-genome decoy panel

  weakest true signal   39 bits   (lowest median self-recall across clusters)
  noise ceiling         18 bits   (worst off-target hit any profile reached)
  -> any cutoff in 19..39 works;  30 chosen, erring toward specificity
```

`build_sp_pssms.py --decoys <panel.fna>` reports exactly this window. If the
noise ceiling reaches the true signal there is **no safe cutoff**, and the
answer is to split the clusters further rather than pick a number and hope — a
19-residue profile that has to run at 20 bits will fire on anything.

The same reasoning applies to any short feature: viroporins, small hydrophobic
proteins, accessory ORFs under ~60 aa.

**`copy_num`** only where you are confident. Note that a module collapsing many
genera may legitimately see a feature absent in most of them; omit `copy_num`
rather than guarantee a copy that some members do not have.

## Flagging model-proposed citations

A DLIT means a paper a curator read, which establishes the function or the
coordinates the feature claims. A citation a language model proposed is not
that, even when the paper is real and topically correct — nobody has confirmed
it supports the claim.

So any citation a model supplied is listed in `PMID_claude_generated` **and not
in `PMID`**. The two lists are disjoint:

```json
"PMID":                  ["9010293"],
"PMID_claude_generated": ["4038520", "8683214", "28276468"]
```

One curated citation, three proposed. A feature whose citations are all
model-proposed carries **no `PMID` key at all** — which is the normal state for a
newly built module, and is what Matonaviridae and Hepeviridae both ship.

`check_dlits.py` treats an id appearing in both lists as an error
(`CITATION IN BOTH FIELDS ... they must be disjoint`), and every one of the
installed modules follows the disjoint form. An earlier revision of this file
described the nested form, where the flag named a subset of `PMID`; that was
wrong, disagreed with the checker and with the shipped data, and is the reason
to state it plainly here.

**The test is who supplied the id, not how much work went into it.** A model
that searched PubMed, matched the title and read the abstract has proposed a
citation, not curated one. Reading does not clear the flag. The only thing that
clears it is a human curator deleting the entry from the generator's registry,
which is why the registry is the single place the provenance lives.

Drive it from a registry in the generator rather than hand-editing the JSON, so
the flag cannot drift from reality:

```python
CLAUDE_PMIDS = {
    "4038520": "IHNV mRNA species (Novirhabdovirus NV)",
    ...
}

if pmid:
    ids      = [str(p) for p in pmid]
    proposed = [p for p in ids if p in CLAUDE_PMIDS]
    curated  = [p for p in ids if p not in CLAUDE_PMIDS]
    if curated:                                  # omit the key entirely when
        d["PMID"] = curated                      # every id is model-proposed
    if proposed:
        d["PMID_claude_generated"] = proposed
```

When a curator has read a paper and confirmed it, delete that entry from the
registry; the flag disappears on the next regeneration. That makes clearing the
flag an explicit act rather than something that happens by editing away an
inconvenient line.

Verify with:

```bash
python3 scripts/check_dlits.py --json <T>_Viral_PSSM.json
```

which catches DOIs in the PMID field, ids that do not resolve in PubMed (for a
model-proposed citation, that is the fabrication case), and disagreement between
the two lists — then prints every citation with its title so relevance can be
eyeballed.

**Do not put a model-proposed citation into the vocabulary.** The table has no
provenance column, so a flagged citation entering it silently becomes a curated
one. Verify it first, or leave the `PubMed IDs*` cell empty.

## close_genomes: name every rep contig, and say what is in it

One entry per `Viral-Rep-Contigs/<Module>.<n>.dna`, keyed by the filename.

`genome_ids` identifies the sequence the file holds. Both conventions are in
use upstream — BV-BRC genome ids (62 entries) and **versioned GenBank
accessions** (15 entries) — and either is fine as long as it is filled in.
A multi-segment exemplar lists one accession per segment, comma-separated and
in the same order as the records in the file:

```json
"Dichorhavirus.1.dna":  { "genome_ids": "NC_009608.1,NC_009609.1", ... }
"Trirhavirus.1.dna":    { "genome_ids": "PQ653982.1,PQ653983.1,PQ653984.1", ... }
```

Leaving it blank costs you provenance you cannot reconstruct later: the `.dna`
header is free text, and nothing else records which release of which record the
routing was built against.

`genome_name` is what the annotator prints in the feature table as the matched
reference, so a blank one produces a nameless match column.

**The declaration and the directory must agree in both directions.** A `.dna`
file that no entry declares is still BLASTed — the annotator globs the
directory — so it silently routes genomes while carrying no name, and the next
regeneration of the JSON does not know it exists. That happened here: a rep
contig added to close a routing gap stayed installed while a later regeneration
dropped its entry. `install_module.py --check` now reports both directions plus
any entry with a blank `genome_ids` or `genome_name`.

## Extension flags on mature peptides

`upstream_ext` and `downstream_ext` say the annotator may scan outward for a
start codon or a stop codon. **A cleavage site is neither.**

A precursor rarely yields just one product. Glycoproteins in particular are
cleaved more than once — GPC into a stable signal peptide plus Gn and Gc,
HA0 into HA1 and HA2, S into S1 and S2, filovirus GP into GP1 and GP2 — and
The vocabulary carries all of those as separate `mat_peptide` rows. Work out where
each product sits in the chain, then set the flags from its two termini:

| Product | N-terminus | C-terminus | `upstream_ext` | `downstream_ext` |
|---|---|---|---|---|
| first (signal peptide) | initiating Met | cleavage | **1** | **0** |
| internal (e.g. Gn, HA1) | cleavage | cleavage | **0** | **0** |
| last (e.g. Gc, HA2) | cleavage | protein's true end | **0** | **0** |

Which collapses to a rule worth remembering: **only the first product of a
precursor ever gets `upstream_ext: 1`; every mature peptide gets
`downstream_ext: 0`.** The last product is 0 as well, because a mature peptide
never includes the stop codon even when it runs to the protein's real end.

Sanity check on any multi-product precursor: the products should **tile** it.
Rabies G, which has a single cleavage:

```
G         3318..4892      524 aa + stop
SP        3318..3374       19 aa      starts where G starts
G_mature  3375..4889      505 aa      starts one base after SP ends
                                      ends 3 nt short of G -- the stop codon
```

Gaps or overlaps between adjacent products mean a flag is wrong or a cleavage
coordinate is off.

Leaving the default in place is not a cosmetic error. `G_MAT` kept
`upstream_ext: 1`, so the annotator scanned upstream from the mature
N-terminus, found the *precursor's* start codon, and emitted the mature peptide
on the same coordinates as the precursor:

```
CDS          G         3317  4891   Envelope glycoprotein precursor
mat_peptide  G_mature  3317  4888   Mature envelope glycoprotein     <- wrong
```

With `upstream_ext: 0` the same genome gives the real cleavage:

```
CDS          G         3317  4891   MVPQALLFVPLLVFPLCFG...  524 aa
mat_peptide  G_mature  3374  4888   KFPIYTIPDKLGPWSPID...   505 aa
```

19-residue signal peptide, cleaved before `KFPIY` — the documented rabies G
processing. Nothing errored in the broken state; the feature was simply wrong
in every genome the module annotated.

Two consequences elsewhere, both worth copying into a new module:

- Do not expect Met at a cleaved N-terminus. `validate_calls.py` keys its
  `NO_MET` check off `upstream_ext == 0` for exactly this reason.
- A mat_peptide sharing a boundary with its parent CDS is the symptom to watch
  for. If mature and precursor start at the same coordinate, the extension flag
  is wrong.

## The annotation string is a key, not a label

**Within one module, no two features may share an `anno` string** if both
declare `copy_num` and a `CDS` or `mat_peptide` `feature_type`.
`viral_genome_quality.pl` builds its `%essential` hash keyed by the annotation
string, not by the feature key:

```perl
$essential->{$anno}->{max_len}  = ...->{max_len};
$essential->{$anno}->{min_len}  = ...->{min_len};
$essential->{$anno}->{copy_num} = ...->{copy_num};
```

Two features sharing a string means the second one read **silently overwrites
the first one's length window**, and `%anno_count` then sums both features'
calls against the single surviving `copy_num`.

The symptom is distinctive, and it looks like a biology problem when it is not:

- every genome flagged `Genome has too many HSPs for: <string>; Count = 2`
- length flags that **flip between "too short" and "too long"** across the run,
  because Perl randomises hash order and the script runs once per genome, so a
  different feature's window wins each time

Tobamovirus shipped `REP126` and `REP183` both as `Nonstructural polyprotein`
and scored **0% good on 89 genomes** for this reason alone — 52 genomes where
the 183K window won and the 126K call read too short, 38 where the 126K window
won and the 183K call read too long, one flag per genome. Nothing was wrong
with the profiles, the cutoffs or the genomes.

A readthrough or nested pair is still two proteins and needs two strings.
Togaviridae is the worked precedent: `nsP1234` is `Nonstructural polyprotein`
and `nsP123` is `Nonstructural polyprotein P123`. Where the two products have
genuinely different functions, say so instead of appending a name — Tobamovirus
now uses `Methyltransferase and helicase replication protein` (126K) and
`RNA-dependent RNA polymerase` (183K), which also reuses an existing
vocabulary string for the readthrough product.

`scripts/check_annotations.py` fails with a non-zero exit on this. Run it
before installing.

## Two failure modes the schema does not prevent

**A feature declared with no PSSM is never called.** Nothing errors; the
annotator simply never emits it. This is fine when deliberate (a documented gap)
and invisible when accidental.

**A PSSM with no JSON entry is never loaded.** Also silent.

`scripts/install_module.py --check` reports both and refuses to install. Run it
every time.

## Generating the block

Write a generator script rather than editing the JSON by hand, and keep it in
the working directory so the JSON can be regenerated byte-identically. Lengths
come from the collections, `anno` from your vocabulary decisions, PMIDs from the
literature you consulted. Hand edits to a generated file are lost on the next
regeneration and are the usual reason a JSON and its collections drift apart.

Validate before installing:

```bash
python3 -c "import json;d=json.load(open('<T>_Viral_PSSM.json'));print(len(d),'modules')"
python3 scripts/check_annotations.py --json <T>_Viral_PSSM.json
python3 scripts/install_module.py --workdir . --repo $LOWVAN_DATA_DIR --check
```

## One canonical format, always

`Viral_PSSM.json` is written by tools in two languages and lives in git, so it
needs exactly one serialisation or every write reformats the file and buries the
real change.

The canonical form is what Perl produces with:

```perl
JSON::XS->new->pretty->canonical->encode($data)
```

3-space indent, `" : "` between key and value, keys sorted at every level,
trailing newline. The indent width is fixed at 3 inside JSON::XS and cannot be
configured, so that is a fact about the tool, not a style choice.

**`canonical` is the part that matters.** Perl randomises hash order, so a plain
`decode`/`encode` round-trip of the master reorders all 7,373 lines
*differently on every run* — verified here, three runs, three different files.
Sorting keys is what makes the file diffable at all.

Python matches byte for byte:

```python
json.dumps(normalize(obj), indent=3, separators=(",", " : "),
           sort_keys=True, ensure_ascii=False) + "\n"
```

where `normalize` folds an integral float to an int, because Perl has no
separate float type for one and writes `1` where Python writes `1.0`. That is
the only place the two languages disagree on this schema.

Use `scripts/json_canon.py` rather than reimplementing it:

```bash
python3 scripts/json_canon.py --check  Viral_PSSM.json <T>_Viral_PSSM.json
python3 scripts/json_canon.py --write  <T>_Viral_PSSM.json
```

`install_module.py` writes the master canonically and says so when a write
reformats a file that was not already canonical — that happens once, and it
looks alarming in a diff if nobody warned you. A module generator should emit
the same format directly; keep the serialiser inline there so the generator has
no dependency outside its own working directory.

## Merging into the master

`install_module.py` merges module blocks into the repo's `Viral_PSSM.json`,
preserving existing key order and appending new modules, and writes a `.bak`
first. It replaces a module wholesale if the key already exists — so a
re-install after a rebuild is safe and idempotent.
