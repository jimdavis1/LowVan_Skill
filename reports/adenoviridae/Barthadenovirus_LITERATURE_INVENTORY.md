# Barthadenovirus (= Atadenovirus) — literature inventory

**Genus name.** ICTV renamed *Atadenovirus* to **Barthadenovirus** in the 2024
ratification (proposal 2023.008D). BV-BRC files the genomes under the new name;
every paper cited below predates the rename and says *Atadenovirus*. The module
is named `Barthadenovirus` after the current name, and `Atadenovirus` is
recorded here as the synonym a reader will meet in the literature.

**What this file is.** A per-protein check of three things the curator asked
for, in this order: (1) is every protein the literature describes present in
the module, (2) are its coordinates right, (3) does it have a real function
rather than a position label.

**Coordinates below are DATABASE annotations, not literature boundaries.**
See §0 for what that means and how far each one can be trusted.

| Reference | Accession | Length | Role here |
|---|---|---|---|
| Ovine adenovirus 7 (OAdV-7 / OAV287) | `NC_004037` | 29,576 bp | genus prototype; the RH and LH genes are defined on it |
| Snake adenovirus 1 (SnAdV-1) | `NC_009989` | 27,751 bp | reptilian branch; carries a *different* right-hand set |

Citations are in `PMID_claude_generated`, never `PMID`. Every one below was
resolved against PubMed and its title read before it was written here.

---


## 0. Where these coordinates come from, and how far they can be trusted

**GenBank and RefSeq records carry the same kind of error BV-BRC does.** A
`/product=` or a `join()` is a submitter's claim with at most a curation pass
over it. This family proves it twice over:

- `NC_004037` annotates the LH3 locus as **"E1B 55K"**, and atadenoviruses
  have no E1B region at all.
- `NC_001958` is a **PROVISIONAL RefSeq** — its own header says it "has not
  yet been subject to final NCBI review".

So the tables in this file are a **consistency check**, not an authority. They
answer "does this collection hold one gene or two mixed together", which a
collection cannot answer about itself. Where the literature establishes a
boundary, the literature wins and any disagreement is written down.

### Evidence grade, per source

| Grade | What it is | Sources for this genus |
|---|---|---|
| **Experimental** | transcript mapping, structures, proteomics of real virions | see table below |
| **Computational** | genes called by homology to other adenoviruses | the genome papers' own gene tables |
| **Database only** | a feature table with no paper behind the specific boundary | RefSeq re-annotation, BV-BRC |

| Source | PMID | Grade | What it actually establishes |
|---|---|---|---|
| Vrati *et al.* 1996, the OAdV-7 genome paper | 8659111 | **Experimental (proteomics)** + computational | The primary genome description. Reports **proteinase cleavage sites measured on purified virions**, and notes that some "**did not conform to the known specificity of the human Ad proteinase**". Also states that structural proteins **V and IX are absent** from the genome, that a **new processed 28 kDa virion polypeptide** is encoded opposite the proposed E1A region, that two structural coding sequences were then **unidentified**, and that the OAV penton **lacks the RGD motif** human adenoviruses use to engage integrins. |
| Khatri & Both 1998, transcripts and promoters | 9614874 | **Experimental** | Transcription start sites and **RNA splice junctions** mapped for the left-hand end, IVa2, E2, P32K and E4 regions. Left-hand end has three temporally distinct promoters, E4 two; E2 has a single promoter, unlike mastadenoviruses; p32K and IVa2 each have their own. **This is the real boundary authority for those regions.** |
| Both 2002, right-end transcription map + F-box | 12504581 | **Experimental** | Eight transcripts from a terminal-repeat promoter, five potentially bicistronic, covering RH1/RH2/RH4/RH6; and the F-box motif itself. |
| Marabini 2021 / Menéndez-Conejero 2017 | 33789897, 28943338 | **Experimental** | LH3, p32K and LH2 as virion components; LH3's fold and capsid-binding motif. These **supersede the 1996 paper's tentative "E1A/E1B" reading of the left-hand end**, which is the origin of the "E1B 55K" name still in the database. |

**The provenance of the "E1B 55K" error, traced:** the 1996 paper proposed
only *tentatively* that the left-hand group "probably represents the E1A/E1B
regions". That tentative reading hardened into a `/product="E1B 55K"` in the
database, and the 2017/2021 structural work identified the same protein as
LH3, a capsid cementing protein. The literature corrected itself; the database
did not follow.

**Not yet obtained.** PMIDs 8659111, 9614874, 12504581 and 18166240 are **not
open access** and only their abstracts were read. The cleavage sites in
8659111 and the splice junctions in 9614874 are the authoritative boundaries
for this genus and they are **not yet in this module** — the figures here
still come from the database record. This is a known, named gap, not a
silently accepted one.


## 1. Completeness — what the module was missing

Three right-hand genes were **absent from the feature set entirely**, and the
reference genome has all of them:

| Gene | OAdV-7 coordinates | aa | Was | Now |
|---|---|---|---|---|
| RH1 | `complement(29130..29459)` | 109 | 7 deposits → UNASSIGNED, 3 → mixed UNC | declared `RH1` |
| RH6 | `complement(26688..27161)` | 157 | 3 deposits → mixed UNC | declared `RH6` |
| RH0 | not in OAdV-7; 3 deposits in the dump | ~187 | 1 → mixed UNC, 1 → UNASSIGNED | declared `RH0` |

The right-hand region is numbered **RH0–RH6** and the module declared only
RH2–RH5. The missing three were not reported as a gap anywhere: with no key to
bind to, they were swept into one accessory cluster (`UNC3`) that mixed RH6,
RH1 and RH0 — three different genes in one collection, which would have trained
one profile on all three.

## 2. Coordinates — the two that were wrong

**`LH3` is the gene the reference genomes call "E1B 55K".** This is the single
largest correction in the module.

- OAdV-7 annotates `1916..3064`, 382 aa, as **E1B 55K**. It sits in the
  *left-hand* region, immediately after LH1 (`1143..1505`) and LH2
  (`1490..1873`) and before IVa2 — i.e. exactly where LH3 belongs. SnAdV-1
  annotates the same locus `2169..3290`, 373 aa, the same way.
- Measured: all 13 sequences deposited as "E1B 55K" / "E1B-2" / "E1B protein,
  large T-antigen" hit the `LH3` collection at **42–77 % identity over the full
  length**, at 370–393 aa against LH3's 361–383 aa.
- The structural literature predicts precisely this confusion: LH3 has a
  trimeric β-helix phage-tailspike fold and binds the capsid through a
  triskelion *identical* to mastadenovirus protein IX, and mastadenovirus
  **E1B-55K was exapted from the atadenovirus LH3** — an ancient duplication
  (PMID_claude_generated 33789897; fold first solved in SnAdV-1,
  PMID_claude_generated 28943338).

So the deposited name is a claim, and here it is a misleading one:
atadenoviruses have **no E1B region at all**. Without the fix, `LH3` would have
shipped a profile built from 4 sequences while 13 real LH3s sat in an unnamed
cluster.

**`RH3` shipped empty, and it is not an absent gene.** Four RH3 deposits
(160 aa) exist in the dump, and consolidation moved every one of them into
`RH4`, along with 7 of RH2's 22. RH1/RH2/RH4/RH6 are *related paralogues* —
the F-box paper says so in as many words — so homology cannot be allowed to
merge them. They are now a declared paralogue set, as FIBER1/FIBER2 already
were.

**Splicing.** Three features in OAdV-7 carry exactly one intron each, which is
what the splice caller models:

| Feature | OAdV-7 | Note |
|---|---|---|
| 33K | `join(20891..21062,21168..21397)` | 133 aa |
| pTP | `join(7349..9121,11883..)` | 595 aa |
| U exon | `complement(<22106..22266)` | 54 aa; "continues upstream, probably spliced" |

IVa2 is annotated "protein-coding region continues upstream and is probably
spliced" in **both** references — flagged, not yet modelled.

## 3. Functions — what is actually known

| Feature | Function | Evidence |
|---|---|---|
| **RH1, RH2, RH4, RH6** | **F-box proteins.** They carry an F-box motif and are expected to drive degradation of cellular regulatory proteins through SCF ubiquitin-ligase complexes — a strategy that distinguishes this genus from every other adenovirus genus and from most viruses generally. At least three are non-essential for replication in vitro. | PMID_claude_generated 12504581 |
| **LH3** | **Capsid cementing protein.** Trimeric β-helix, bacteriophage-tailspike fold; binds the capsid surface via the same triskelion motif mastadenovirus pIX uses. Mastadenovirus E1B-55K was exapted from it. | PMID_claude_generated 33789897, 28943338 |
| **p32K** | **Genus-specific minor coat protein.** α-helical clusters beneath the capsid facet. `complement(238..1098)`, 286 aa in OAdV-7. | PMID_claude_generated 33789897, 28943338 |
| **LH2** | **Structural component**, previously unrecognised as one; an extended polypeptide that drives a large conformational change in the internal vertex protein IIIa relative to mastadenoviruses. | PMID_claude_generated 33789897, 28943338 |
| **E4.2, E4.3** | Share a conserved motif with each other **and with mastadenovirus E4 34K**, per the reference annotation's own note. Function not established. | `NC_004037` annotation |
| pVII, pX | Deposited as "core protein 1" and "core protein 2". Confirmed by position in both references: OAdV-7 pVII `13313..13648`, pX `13669..13884`. Bound by homology, never by that name. | `NC_004037`, `NC_009989` |
| **RH0** | One deposit is annotated "F-box related protein RH0", which fits the RH family above, but RH0 is in neither reference genome and no paper covers it. **Keeping the number.** | — |
| RH5, E4.1, LH1 | No function established in the literature. RH5 carries no F-box note in the reference annotation, unlike RH1/2/4/6. **Keeping the numbers** — a functional string the literature does not support is worse than a position label. | — |

**SnAdV-1 carries a different right-hand set** — `ORF1` (112 aa) and `105R`
(178 aa, "similar to 105R of TSAdV-1") in place of most RH genes, with only RH1
shared. The reptilian atadenoviruses are not simply OAdV-7 with substitutions,
and a module built only on the ovine/bovine layout will not call them. Measured
against the dump before the module ships.

## What was probed and not found

Honest edge of this check, so the next reader can see it:

- **Protein IX.** An earlier note in `probes.json` said Barthadenovirus has pIX
  because "the atadenovirus capsid structure resolved it together with p32K".
  **That was wrong.** The structure paper resolves the genus-specific proteins
  LH2, p32K and LH3, and its point is that *LH3 occupies pIX's structural role*
  — not that pIX is present. Neither reference genome annotates a pIX. The
  probe is corrected to ABSENT.
- No literature was found assigning a function to RH5, RH3, LH1, E4.1, or to
  the SnAdV-1-specific ORF1 and 105R.

---

## Final state (pass 10, 2026-11-03)

Every probe scores **PRESENT** except protein IX (decided ABSENT, above).
Lengths checked against `NC_004037`: **29 of 30 features match.**

The paralogue audit moved 3 sequences, each needing reference length *and*
homology to agree:

| Move | Length | Own ref | Target ref |
|---|---|---|---|
| LH1 → LH3 | 382 aa | 120 | 382 |
| RH1 → RH4 (×2) | 138 aa | 109 | 157 |

`LH1` is now 117–130 aa with a median of exactly 120, `LH3` 361–393 with a
median of exactly 382, and `RH1` a median of exactly 109. `RH3` and `RH0`
carry sequences for the first time.

### The one feature still flagged

**`RH2` is off by 30 %** — reference 176 aa, collection median 124, range
103–186 across 7 sequences. Genome `130499.16` deposits **two different
proteins both called RH2**, at 103 aa and 186 aa, and genome `114430` deposits
two at 158 and 176.

The audit deliberately did **not** move them. The length signal points at RH1
(109 aa) for the short ones, but the homology signal does not agree, and this
module only moves a sequence when both do — inside an F-box paralogue family,
homology is precisely what cannot separate the members, and length alone
cannot tell a short paralogue from a truncated deposit. Recorded as an open
item rather than resolved by guesswork. It affects 7 sequences in a 54-genome
module and should be settled per species before `RH2`'s profile is trusted.

### Module partitioning — not yet decided

55 genomes: **24 ruminant, 23 bird, 7 reptile, 1 marsupial.** SnAdV-1
(`NC_009989`) carries a different right-hand set — `ORF1` (112 aa) and `105R`
(178 aa) in place of most RH genes, with only RH1 shared. Genome organisation
is what defines a module, so this is a real question and not a phylogenetic
one. Building as a single module and measuring routing first; a three-way
split would leave the reptilian module at 7 genomes, which is below what the
greedy reference picker can serve.
