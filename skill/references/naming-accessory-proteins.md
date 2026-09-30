# Naming an accessory protein that the source called by its mass

BV-BRC annotated 26 closterovirid, bromovirid and betaflexivirid accessory
features as `p6 protein`, `p22A protein`, `14 kDa protein` and the like. Those
are masses, not descriptions, and `check_annotations.py` flags every one: *the
annotation is just the symbol; give a functional description and keep the symbol
in `gene_symbol`.*

Resolved 30 September 2026. The procedure below is the part worth reusing.

## The procedure

**1. Find out what the feature actually contains before looking anything up.**
A mass bin is not a homology group. Two checks, both cheap:

```bash
# which species sit in which alignment
grep -h '^>' PSSM-Alignments/<FEAT>/corrected_alis/*.fa

# is the bin ONE homology group, or several? all-vs-all within the feature,
# union-find over alignment ids
blastp -query <FEAT>.faa -db <FEAT>.faa -evalue 1e-3 -outfmt '6 qseqid sseqid'
```

Every one of the 26 was species-partitioned at the alignment level, and **6 of
26 turned out to be two to five non-homologous groups sharing a mass**. Naming
those after one member's function mislabels the rest.

**2. Only then go to the literature, and attribute it to the group you have.**
The trap is real: the literature says CTV `p20` is the functional homolog of BYV
`p21`, *not* of BYV `p20`. `Closterovirus_P20` holds CTV p20 (54 sequences, a
silencing suppressor) and BYV p20 (9, an Hsp70h-interacting vascular-transport
protein) in two separate alignments that share no HSP at 1e-3. They are two
different proteins. Had the bin been named from the CTV literature alone, 9
sequences would have been given a function they do not have.

**3. Apply the bar.** A functional string needs both:

- the feature is one homology group, **or** one group holds **≥ 90%** of its
  sequences, and
- a function is documented for that group — not a localization, not a
  prediction, not a phenotype of a different species' positional homolog.

Otherwise it is `Uncharacterized lineage-specific <SYM> protein`. That is an
existing S1 Table form, and `PER_FEATURE` in `check_annotations.py` exempts it
from the one-symbol-per-string rule, because the symbol is the whole point.

**4. Prefer an existing string to a new one, and measure membership.**
`Small hydrophobic protein` already exists for the paramyxovirus and
pneumovirus SH. The closterovirid hallmark array carries several small
hydrophobic genes per genome, so the string fits — but which of the Crinivirus
small features belong to that class is a measurement, not a guess:

```
feat   seqs  med len   GRAVY  max TM window (19 aa)   verdict
P6B       7       51    0.91                   2.37   hydrophobic
P7A      42       57    0.71                   2.38   hydrophobic
P7B       8       65    0.76                   3.12   hydrophobic
P8B      16       67    0.37                   1.61   hydrophobic
P8A      32       78   -0.05                   0.86   SOLUBLE -- excluded
```

`P8A` is the soluble one and stays uncharacterized. That is not a judgement
call dressed up as a number: LIYV `P5` is described as an ER-localized integral
transmembrane protein and LIYV `P9` as *soluble*, and the measurement picks out
which of these features is which. Kyte-Doolittle GRAVY plus the maximum
19-residue window; none of the five is homologous to any other at 1e-3.

Reusing the string forced one change to the checker. `want` is built from every
taxon *except* the module's own, so a per-taxon row can never satisfy it, and
four non-homologous Crinivirus genes cannot all be `SH` without losing the
distinction. `Small hydrophobic protein` joined `PER_FEATURE` for that reason,
with the reason in the source.

## Outcome

12 features took a functional string, 14 became uncharacterized
lineage-specific. Across all fourteen Alsuviricetes modules the checker now
reports zero VARIANTs, zero symbol mismatches and zero style flags.

| verdict | features |
|---|---|
| `RNA silencing suppressor <SYM> protein` | Closterovirus P23; Ampelovirus P20B; Crinivirus P22A, P22B, P22C, P25; Bromoviridae 2b |
| `Host range determinant <SYM> protein` | Closterovirus P13, P18 |
| `Membrane-associated viroporin P33 protein` | Closterovirus P33 |
| `Movement protein` / `Mov` | Closterovirus P6 |
| `Small hydrophobic protein` | Crinivirus P6B, P7A, P7B, P8B |
| `Uncharacterized lineage-specific <SYM> protein` | Closterovirus P19, P20; Ampelovirus P7, P20A, P21, P24; Crinivirus P8A, P27; Quinvirinae ORF2a, ORF5a, P14 |

The five that a single-species reading would have got wrong:

- **Closterovirus P20** — two proteins, as above.
- **Ampelovirus P24** — GLRaV-1 p24 is a well-documented silencing suppressor,
  but it is one of three groups and holds only 42 of 62 (67.7%). Below the bar.
- **Ampelovirus P7** — five groups, largest 55%. GLRaV-1 p7 is a small
  transmembrane protein, GLRaV-3 p7 a unique ORF of unknown function.
- **Ampelovirus P21** — two groups; GLRaV-1 p21 is explicitly reported as
  unknown function and GLRaV-3 p21 only as cytosolic.
- **Crinivirus P27** — one group across eight species, but the suppressor
  evidence is for TICV p27 and TICV is not in the collection. Extrapolating a
  function across species on positional grounds is the U-number trap wearing a
  different hat.

**Ampelovirus P20A** is the case that shows the bar doing work in the easy
direction: one clean homology group, GLRaV-3 only, and the only published
result is that it associates with microtubules. A localization is not a
function, so it stays uncharacterized.

## What this does not settle

Six features remain mass bins holding non-homologous groups. Naming them
honestly is done; **splitting them is not**, and until they are split those
features will keep calling two or more different proteins under one key. The
minority groups are small — 2 to 9 sequences in the Closterovirus cases, larger
in Ampelovirus P7 and P24 — and splitting is a change to shipped modules that
needs its own measurement and its own decision.
