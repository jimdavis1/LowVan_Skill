# Siadenovirus — literature inventory

A per-protein check of the three things the curator asked for: (1) is every
protein the literature describes present, (2) are its coordinates right,
(3) does it have a real function rather than a position label.

Coordinates are read from the annotated reference genome, not assumed:

| Reference | Accession | Length |
|---|---|---|
| Turkey adenovirus 3 (TAdV-3, *Siadenovirus gallopavotertii*; formerly turkey haemorrhagic enteritis virus, THEV) | `NC_001958` | 26,263 bp |

Citations are in `PMID_claude_generated`, never `PMID`. Every PMID below was
resolved against PubMed and its title read before it was written here.

---

## 1. Completeness — the genus-specific set

The genus is defined by five genus-specific genes: **ORF1 (sialidase), hyd, E3,
ORF7 and ORF8**. All five are in the module:

| Gene | TAdV-3 coordinates | aa | Module status |
|---|---|---|---|
| ORF1 / sialidase | `399..1952` | 517 | `SIALIDASE`, 50 features |
| hyd | not annotated in `NC_001958` | 68–75 in deposits | `HYD`, 15 features |
| E3 | `21213..22115` | 300 | `E3`, 47 features |
| ORF7 | `24509..25165` | 218 | was `UNC1` (16 seqs, median 225 aa) → **name it ORF7** |
| ORF8 | `complement(25201..25698)` | 165 | was `UNC2` (11 seqs, median 174 aa) → **name it ORF8** |

The two largest accessory clusters are therefore not unknowns at all — they are
two of the five genus-defining genes, and the median lengths match the
reference to within a few residues. The remaining two clusters map to the
reference's other numbered ORFs: `UNC3` (8 seqs, median 111 aa) against
ORF4 `1964..2311` (115 aa), and `UNC4` (4 seqs, median 142 aa).

TAdV-3 also annotates ORF2 (100 aa), ORF3 (124 aa), ORF5 (165 aa) and ORF6
(120 aa) as hypothetical proteins. These are below the module's accessory
floor in the current dump; recorded here so that they are not rediscovered as
novelties later.

## 2. Coordinates

**`core protein I` and `core protein II` are pVII and pX.** TAdV-3 annotates
`12346..12708` (120 aa) as "core protein I" and `12711..12887` (58 aa) as
"core protein II", positioned between the penton base and pVI — which is pVII's
and pX's position. They are bound by homology, never on those strings.

**The reference's own E3 call is weak.** `NC_001958` annotates E3 at
`21213..22115` with the note *"identification based on very limited
similarity"*. The database is flagging its own uncertainty, and the module
inherits it: `E3` holds 47 features but the identification rests on the
deposits, not on a paper. Recorded as a known soft spot rather than laundered
into a fact.

**Splicing.** A 2026 RACE/Sanger transcriptional map of TAdV-3 gives
full-length annotation for 14 of 23 predicted transcripts and reports a
**conserved tripartite leader shared across all late transcripts**, alternative
leader exons (DBP mRNA carries a bipartite leader while Hyd mRNA uses the
tripartite leader), and shared polyadenylation signals between ORF1 and Hyd
(PMID_claude_generated 41960683). As in Mastadenovirus, the tripartite leader
is 5′ UTR rather than coding, so it does not force multi-exon CDS models.

## 3. Functions

| Feature | Function | Evidence |
|---|---|---|
| **SIALIDASE (ORF1)** | The genus-defining gene — the *sia* in Siadenovirus. Its product carries a region of high similarity to **bacterial sialidases**, and it is a confirmed structural component of the TAdV-3 virion. Whether the conserved sialidase domain is catalytically important is **not known**. | PMID_claude_generated 36016389 |
| **FIBER** | Receptor binding. TAdV-3 uses **sialic acid on N-linked membrane glycoproteins of B cells** as its cellular receptor, and this is required for infection. The fibre head is a twisted β-sandwich and binds sialyllactose. | PMID_claude_generated 32459612 |
| **K100** | Hexon assembly chaperone / late host shut-off protein. See the correction below. | `NC_001958` `18184..20226`, 680 aa |
| **HYD** | A hydrophobic protein of **unknown function**. Named for the property, not a function. Its mRNA uses the tripartite leader and shares a polyadenylation signal with ORF1. | PMID_claude_generated 41960683 |
| **E3** | Possible transcriptional regulator, but no functional study was found and the reference's own identification is flagged as weakly supported. **Not asserting a function.** | — |
| **UXP** | U exon protein; late l-strand product. | PMID_claude_generated 17881437 |
| **ORF7, ORF8** | Genus-specific, **function unknown**. **Keeping the numbers** — these are two of the five defining genes and their existence is well established, but no paper assigns either a function, and inventing one would be worse than the position label. | — |

## 4. A feature that should not have existed

**`SHUTOFF` was declared as a separate feature and it was wrong.** It shipped
with zero sequences, which looked like a missing protein. It is not:

- The 8 deposits reading "shutoff protein" are **703 aa**, and consolidation
  placed every one of them into `K100` — correctly.
- The late **100K protein is** the host shut-off protein. The string is a
  synonym for a feature the module already declares, not a gene of its own.

Declaring both would have put two profiles on one locus. `SHUTOFF` is removed
and the strings are now listed as `K100` synonyms. This is the opposite of the
LH3 case in Barthadenovirus, and the two together are the argument for checking
*every* empty collection rather than assuming either outcome: one empty
collection was a real protein hiding under a wrong name, the other was a
duplicate of a protein already present.

## 5. Confirmed absent

**Protein IX.** Siadenovirus is the one genus with no pIX. `NC_001958`
annotates none, and the cross-genus capsid-cement comparison says so directly
(PMID_claude_generated 33789897). Not a gap — recorded as a decided ABSENT.

## What was probed and not found

- No functional study for E3, hyd, ORF2–ORF8.
- No evidence that the sialidase domain of ORF1 is enzymatically active; the
  literature is explicit that this is open.

---

## Final state (pass 10, 2026-11-03)

Every probe scores **PRESENT** except protein IX, which is the decided ABSENT
above. Lengths checked against `NC_001958`, **17 of 17 features match**:

| Feature | ref aa | n | median | range |
|---|---|---|---|---|
| POL | 1112 | 23 | 1112 | 1107–1112 |
| HEXON | 906 | 13 | 912 | 906–922 |
| PTP | 597 | 23 | 600 | 577–625 |
| K100 | 680 | 21 | 665 | 630–706 |
| SIALIDASE | 517 | 27 | 565 | 517–620 |
| PIIIA | 505 | 11 | 505 | 482–513 |
| FIBER | 454 | 28 | 454 | 165–664 |
| PENTON | 448 | 10 | 448 | 445–449 |
| IVA2 | 367 | 15 | 366 | 315–385 |
| DBP | 346 | 16 | 372 | 302–392 |
| E3 | 300 | 19 | 296 | 163–313 |
| K52 | 300 | 12 | 299 | 282–348 |
| PVI | 231 | 12 | 222 | 216–246 |
| **ORF7** | 218 | 16 | 225 | 218–238 |
| PROTEASE | 214 | 12 | 204 | 202–214 |
| PVIII | 200 | 14 | 200 | 153–211 |
| **ORF8** | 165 | 11 | 174 | 164–185 |
| PVII | 120 | 11 | 121 | 120–149 |
| PX | 58 | 10 | 57 | 57–58 |

Accessory clusters are down from 4 to **2**, both of which are the reference's
own numbered hypotheticals (ORF4 at ~111 aa, ORF9 at ~142 aa).

### "ORF7" is a position, not a protein — measured

21 unique sequences in this dump carry the label `ORF7`/`ORF 7`/`ORF07`, at
**46, 106, 107, 131, 160, 172, 186, 187, 198, 218, 220, 223, 225, 377 and
977 aa**. Only the 218–225 band is the Siadenovirus ORF7; the rest are other
genomes' seventh reading frames, and a name rule on the string would have
built one profile out of all of them.

The collection was formed by **homology clustering first**, and the label was
read off afterwards only to name what had already been found. The other
`ORF7`-labelled sequences remain unplaced, which is the correct outcome and
not a gap. This is the U-number trap behaving exactly as the triage reference
describes it.
