# Bornaviridae — what the literature says, with verified PMIDs

Every PMID below was resolved against PubMed (esummary/efetch) in this build and
the title printed; none is written from memory. All of them are
`PMID_claude_generated` — no curator has read them for this module.

## The family

| | |
|---|---|
| ICTV profile | **34227935** — *ICTV Virus Taxonomy Profile: Bornaviridae* (J Gen Virol 2021). Enveloped, linear negative-sense non-segmented RNA, ~9 kb; hosts mammals, birds, reptiles, fish. |
| nuclear replication + splicing | **8035500** — Schneider, Schneemann & Lipkin, *RNA splicing in Borna disease virus, a nonsegmented, negative-strand RNA virus* (J Virol 1994): two introns identified in the 7.1- and 2.8-kb transcripts. |
| the same year, independently | **7831966** — *RNA splicing contributes to the generation of mature mRNAs of Borna disease virus* (Virus Res 1994). |
| intron III | **11172106** — Cubitt, Ly & de la Torre, *Identification and characterization of a new intron in Borna disease virus* (J Gen Virol 2001): a third intron by alternative 3′ splice-site choice, coding for two predicted proteins of 8.4 kDa and 165 kDa (p165, a deleted form of L). **mRNA detected; the proteins are predicted, not demonstrated.** |
| alternative splicing / SA3 | **11070091** — Tomonaga et al. (PNAS 2000): alternative 3′ splice site SA3, negatively regulated by an exon splicing suppressor inside the polymerase exon. |

### Gene order is NOT the same in every genus

**42358456** — *Continuous in silico screening of fish, amphibians, reptiles and
birds expands the diversity of bornavirids and reveals unexpected genomic
architectures* (Virus Evol 2026) states it explicitly:

    3'-N-X/P-M-G-L-5'          Orthobornavirus
    3'-N-X/P-G-M-L-5'          Carbovirus, Cultervirus
    3'-N-vp1-vp2-X/P-G-M-L-5'  Cartilovirus

**29462190** — Hyndman et al. (PLoS Pathog 2018) reports the carbovirus order
`3'-N-X-P-G-M-L` from Australian carpet pythons, "representing a transposition
of the G and M genes compared to other bornaviruses".

M and G swap places; **no genus lacks any of the six genes** in the material
this module was built from. That is why one module covers the family: the
protein set is identical and only the order differs, and gene order does not
enter a PSSM. Measured here: 16 reference contigs route 100% of the 218
near-complete genomes of all five genus groups (see BUILD_NOTES.md).

## The six proteins

Reference: BoDV-1 RefSeq NC_001607, fetched in this build —
N 54..1166 (370 aa), X 1223..1486 (87 aa), P 1272..1877 (201 aa),
M 1893..2321 (142 aa), G 2236..3747 (503 aa), L gene 3687..8822 with the CDS
annotated 3996..8822 (1608 aa). *Corroboration only; the CDS choice is
discussed under L.*

| key | protein | BoDV-1 | evidence |
|---|---|---|---|
| N | nucleoprotein (p40/p38) | 370 aa | **9188580** |
| X | X protein (p10) | 87 aa | **11796712**, **23100370**, **19397879** |
| P | phosphoprotein (p24) | 201 aa | **9268312**, **15992626** |
| M | matrix protein | 142 aa | **19237566** |
| G | glycoprotein precursor (p57/gp94) | 503 aa | **9557754**, **9420287** |
| L | RNA-dependent RNA polymerase | 1711 aa spliced | **10756058** |

- **N**: **9188580** — Pyper & Gartner, *Molecular basis for the differential
  subcellular localization of the 38- and 39-kilodalton structural proteins of
  Borna disease virus* (J Virol 1997). The abundant 38/39-kDa doublet is the
  nucleocapsid; **the two forms differ by 13 amino acids at the N-terminus of
  the 39-kDa form**, which carries a Pro-Lys-Arg-Arg nuclear localization
  signal; the 39-kDa form is nuclear, the 38-kDa form nuclear + cytoplasmic.
  Corroborated here: BoDV-1 N begins `MPPKRRLVDDTDAM...` — PKRR at residues
  3–6 and the second in-frame Met at residue 14, i.e. 13 residues in.
  See **"what is deliberately not modelled"** below.
- **15221540** — Schneider, Naegele & Staeheli (Arch Virol 2004): p39 and p38
  arise from alternative translation initiation; p38 is almost inactive in the
  polymerase complex and buffers the N:P ratio.
- **30541858** — Kojima et al. (J Virol 2019): N transcripts are **also**
  spliced, and with alternative initiation the N gene can yield at least six
  isoforms with different localisation, one of them ER-targeted by a signal
  peptide exposed by intron removal; the spliced isoforms inhibit the
  polymerase. Demonstrated by expression constructs and recombinant virus.
- **X**: **11796712** (importin-α binding motif directing nuclear import),
  **23100370** (X proteins of bornaviruses interfere with type I interferon
  signalling), **19397879** (Hsc70 controls replication via X).
- **P**: **9268312** (phosphorylated by PKCε and CKII), **15992626** (P is the
  central regulatory element assembling the polymerase complex; X is a
  negative regulator of polymerase activity).
- **M**: **19237566** — *Crystal structure of the Borna disease virus matrix protein (BDV-M) reveals ssRNA binding properties* (2009).
- **G**: see the cleavage section.
- **L**: **10756058** — Walker et al. (J Virol 2000). The 190-kDa L protein,
  detected by infected-rat serum, nuclear, phosphorylated, binds P. Critically:
  *"Expression of this protein from the third transcription unit of the viral
  genome is dependent on a splicing event that fuses a small upstream open
  reading frame in frame with the larger downstream continuous open reading
  frame."*

## G is cleaved by furin — the one proven cleavage in this family

**9557754** — Richt et al., *Processing of the Borna disease virus glycoprotein
gp94 by the subtilisin-like endoprotease furin* (J Virol 1998). ORF-IV encodes
p57, N-glycosylated to gp94; the protein is cleaved by furin, which *"most
likely recognizes one of three potential cleavage sites, namely, an arginine at
position 249 of the ORF-IV gene product"*. The furin inhibitor decRVKRcmk
reduces infectious virus, so the cleavage is required for full biological
activity. **9420287** (entry) and **26332529** (cell-to-cell spread) and
**14581559** (dissemination in neurons) establish the glycoprotein's role;
**33268525** shows that cleavage efficiency, governed by the signal peptide's
interaction with furin, determines mature particle production, and that the
property is conserved across mammalian and avian orthobornaviruses.

**Measured in this build, over the 154-member G collection** (search for
`R-X-[K/R]-R` in the 40–62 % window, extended to the end of the basic run):

    exactly one site        150 / 154   (4 ambiguous, excluded)
    P1  = Arg              150 / 150
    P1' = Asp              149 / 150    (Glu in 1)
    P4  = Arg              150 / 150
    BoDV-1 503 aa form      cut after residue 249 -> 249 aa + 254 aa

The BoDV-1 context is `...RLKRRRR | DTQQI...`, cutting after **R249**, exactly
the residue Richt 1998 names, and the modal products are 249 aa and 254 aa
(consistent with gp51/gp43). The site is enzyme-determined, so it re-locates
itself in every genus rather than being a column number.

**The signal peptide is NOT modelled.** No mature N-terminus of any bornavirus
G has been determined by sequencing, this build has no SignalP licence, and
RefSeq NC_001607 carries no `sig_peptide` or `mat_peptide` feature (checked).
`G_N` therefore runs from the precursor's own initiator Met and *includes* the
~20-residue signal peptide; it is the N-terminal furin fragment, not the fully
matured ectodomain. Stated again in §07 of the coverage audit.

## Cartilovirus carries two extra glycoprotein-derived ORFs

**38028148** — Eshak et al., *Diving deep into fish bornaviruses* (Virus Evol
2023). Little skate bornavirus (*Leucoraja erinacea*) has the architecture
`3'-N-Vp1-Vp2-X/P-G-M-L-5'`; Vp1 and Vp2 are *"ancient duplications of the gene
encoding the viral glycoprotein (G)"* — G/Vp1 28 % identity, G/Vp2 21 %,
Vp1/Vp2 44 %; both lack G's transmembrane domain and both carry a predicted
furin cleavage site. Corroborated here on the dump's single Cartilovirus
genome (3055759.3, BK063518): ORFs in exactly that order, Vp1 = 416 aa at
1172..2419 and Vp2 = 220 aa at 2451..3110, and tblastn of Vp1 against its own
genome hits the G locus at 31 % and Vp2 hits the Vp1 locus at 42 %.

**Vp3 is not modelled.** The same paper proposes a spliced Vp1/Vp2 hybrid, but
calls it *"potential"* and *"tentatively named"*: no protein was observed and no
deposit exists. A prediction is not a product.

## What is deliberately NOT modelled, and why

1. **N p38 (the 38-kDa isoform).** Real and well characterised (9188580,
   15221540), but it is an alternative *translation-initiation* product that is
   an exact C-terminal substring of p40, it is **not separately deposited**
   anywhere in this dump (no 357 aa orthobornavirus N), and no fixed offset is
   derivable: measured over the N collection the second in-frame Met sits at
   residue 14 in 69 members, 16 in 23 and 8 in 12. A nested profile would
   collide with N on N's own coordinates — the failure the skill records for
   parvovirus VP1/VP2.
2. **The spliced N isoforms of 30541858** — same reasons, plus no deposits.
3. **Intron III products p8.4 and p165 (11172106)** — predicted from mRNA; no
   protein demonstrated and no deposits.
4. **Cartilovirus Vp3, and the P–G hybrids predicted for reptilian
   culterviruses (42358456)** — predicted only.
5. **Avian "delgebornaviruses" that lack G, X and/or M (42358456)** — those
   genomes are not in this BV-BRC dump; the module declares no feature it
   cannot see.

## Splicing and L — measured, not assumed

Walker 2000 (10756058) says the L product depends on splicing. Measured here,
by taking every deposited L in the collection and finding the longest exact
match to its own genome in one frame:

    Orthobornavirus  exon 1 contributes 5 aa (MSFHA)  97 deposits
                                        6 aa (MSFHAR) 24
                                        0 aa (collinear)        19
    Cultervirus                         0 aa (collinear)         8
    Carbovirus                          0 aa (collinear)         2
    Cartilovirus                        0 aa (collinear)         1

So the orthobornavirus L is a two-exon protein whose **first exon encodes only
five or six residues**, and L is **collinear in the other three genera**. The
19 collinear orthobornavirus entries are the internal-AUG form: RefSeq
NC_001607 annotates exactly that (CDS 3996..8822 = 1608 aa) inside a gene
annotated 3687..8822, and 3996 − 3687 = 309 nt = **103 codons**, which is the
offset measured here between the 1608 aa and 1711 aa deposits.

**A trap worth recording.** A tblastn collinearity screen of all 745 collection
members against their own genomes reported *every one* as a single collinear
HSP, L included — because the genomic frame upstream of the splice acceptor
happens to read `LCFLTS` where the protein reads `MSFHAS`, a five-residue
difference that is invisible at a 95 %-identity threshold on a 1711-residue
protein. The junction was only found by requiring an **exact** match.
