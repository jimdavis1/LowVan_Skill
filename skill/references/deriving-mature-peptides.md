# Deriving mature peptides from published cleavage sites

For a polyprotein virus, the source annotation is the weakest input in the
chain. Where a cleavage site is established in the literature, **cut there
instead of modelling the product** — the site is better evidence than the
records, and it reaches lineages that have no records at all.

Curator's instruction, 24 September 2026: *"you should be doing all of the
matpeptides this way when possible"*, qualified immediately afterwards —
*"when linked to the literature"*. Conservation measured in your own
collections is **not** sufficient. A published site is.

## When it is worth doing

The case for Pestiviridae NS2 is the clearest in the project:

| | |
|---|---|
| records under an NS2 string | 69 |
| of those, the right length | 34 |
| of those 34, ending in the right residue | **10** |
| pestivirus genomes that cleave NS2-3 at all | a minority — uncleaved NS2-3 is the normal product of non-cytopathic BVDV |
| NS2 features in the CSFV reference `NC_002657` | **none** |

Against that, the site: **Arg1589–Gly1590**, fixed by protein sequencing of two
cytopathic BVDV strains (Tautz & Thiel, *Arch Virol* 2003, PMID 12827468). The
motif `R|GPAVC` is present in **358 of 400** pestivirus polyproteins and
**534 of 556** NS3 records begin `GPAVCK`.

Result: **110 of 111** derived NS2 end in the correct Arg, zero internal stops,
length 394–547 against a 457 aa reference — for a feature with no usable
training data.

## When it fails, and why

The same method was tried on Npro and **rejected on measurement**. Its site is
equally well established — Cys168–Ser169, Cys168 absolutely conserved
(Rümenapf 1998 PMID 9499122; Gottipati 2013 PMID 24146623):

| | calls | ending in the conserved Cys168 |
|---|---|---|
| PSSM | 122 | **96%** |
| derived | 120 | **58%** |

No coverage gain, and the C-terminus collapses. Only 43 of the 100 genomes
called by both agreed on coordinates.

**A derivation is only as good as its anchor's boundary.** NS2 works because
its downstream anchor is NS3's N-terminal Gly — the literature site itself,
called precisely. Npro's anchor is the capsid profile, which does not place
that junction reliably. Npro kept its PSSMs.

## How to decide

1. **Measure terminal-residue conservation per feature first**, from the
   unfiltered dump. It shows which junctions are enzyme-determined and
   therefore derivable. Pestivirus NS3-4A cuts after Leu at 97–99% across four
   junctions; signalase after a small residue; Npro after its own Cys at 98%.
2. **Find the published site.** RefSeq `mat_peptide` coordinates on a reference
   genome count; so does a protein-sequencing paper. A prediction does not.
3. **Check the site transfers.** *Rodent hepacivirus* has published coordinates
   in `NC_021153` and they do **not** transfer: that reference reaches all 58
   failing genomes for C, NS3 and NS5B at ~40% identity and gives **no hit at
   all** for E2 or p7. Homology does not cross 40% identity for a hypervariable
   envelope or a 54-residue hydrophobic peptide.
4. **Validate the derived product**, not just its length: terminal residue
   against the site, internal stops, length distribution against the reference.
5. If it loses to the PSSM, **say so in the feature's own comment** so it is
   not tried again blind.

## The mechanism

The annotator already supports it. A PSSM feature declares
`non_pssm_partner: ["FEAT"]`, and the derived feature declares:

```json
"begin": {"begin_pssm": "P7",  "begin_pssm_loc": "STOP",  "begin_offset": 1},
"end":   {"end_pssm":   "NS3", "end_pssm_loc":   "START", "end_offset":   1}
```

Both sides must declare the partnership or it never fires; `install_module.py`
checks this. `min_len`/`max_len` are **nucleotides**, and `end_offset` is
subtracted. See silent-failures.md.

**A derived feature cannot anchor on another derived feature** — the machinery
anchors on PSSM matches only. That is why Pestiviridae p7 could not be derived:
it would need E2 and NS2 as anchors, and NS2 is itself derived.

## Where this does NOT rescue anything

When the source annotates only the polyprotein *and* no cleavage site is
published for that lineage, there is nothing to do. Ten divergent
Hepaciviridae lineages — *Hepacivirus glareoli*, koala, possum, shrew — carry a
polyprotein annotation and no mature peptides; `NC_038429` has one gene feature
and **zero** mat_peptides, and koala and possum hepacivirus have no RefSeq
record at all. Curator's ruling: *"if you can't find coordinates in the
literature for those mat peptides, and thus cannot build de novo PSSMs from
their sequences, then we are done."*

Record it as a documented ceiling and move on. Do not place a boundary by
assuming a neighbour's length.
