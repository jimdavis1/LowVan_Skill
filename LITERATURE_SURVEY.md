# Nodaviridae — literature survey of proteins, cleavage products and non-canonical expression

Compiled 2026-11-21 for the LowVan Nodaviridae annotation module.

**Verification rule applied throughout.** Every PMID in this document was fetched in this
session from `eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi` (or `efetch.fcgi`) and the
returned title is reproduced verbatim next to the claim it supports. No PMID is written from
memory. Where a claim rests only on an abstract it is marked **(abstract only)**. Where a
length or feature rests only on a GenBank/RefSeq annotation with no supporting paper it is
marked **GenBank-only** and is *not* presented as established.

Evidence levels used in the tables:

| level | meaning |
|---|---|
| **EXP** | experimentally demonstrated (protein detected, mutagenesis, structure, MS, Edman) |
| **SEQ** | inferred from sequence in a published paper (ORF described, not protein-level verified) |
| **GB** | GenBank/RefSeq annotation only, no paper found |

---

## 1. Alphanodavirus (insect nodaviruses)

Members covered: Flock House virus (FHV), Nodamura virus (NoV), Black beetle virus (BBV),
Boolarra virus (BoV), Pariacoto virus (PaV), Newington virus.

### Genome architecture (established)

* **RNA1** encodes **protein A**, the RdRp, from the genomic message.
* During replication protein A makes a 3'-coterminal **subgenomic RNA3** (~370–414 nt, capped,
  not packaged) from the 3' end of RNA1.
* RNA3 carries **two overlapping ORFs**: **B1**, in the *same* frame as protein A (it is
  literally the C-terminal fragment of protein A, translated from an internal AUG), and **B2**,
  in the **+1 frame** relative to protein A.
  * PMID 15163762 — "Nodamura virus nonstructural protein B2 can enhance viral RNA accumulation
    in both mammalian and insect cells." (J Virol 2004). Full text, verbatim: *"Protein B1 is in
    the same reading frame as protein A and thus represents its C-terminal fragment"*, *"whereas
    protein B2 is in the +1 reading frame relative to protein A"*, and *"FHV, 106 amino acids;
    NoV, 134 and 137 amino acids"*. Table 1 footnote (for NoV): *"Residue 1 of protein B1
    corresponds to residue 908 of protein A; protein B1 is identical to the C terminus of
    protein A."*
  * PMID 6495657 — "Sequence of the black beetle virus subgenomic RNA and its location in the
    viral genome." (Virology 1984): RNA3 is 389 b, *"two overlapping reading frames for putative
    proteins of size 10,768 and 11,633 Da"*.
  * PMID 3839022 — "Structure of the black beetle virus genome and its functional implications."
    (J Mol Biol 1985): BBV RNA1 3105 b, protein A 101,873 Da, RNA3 389 b starting 16 bases inside
    the protein A ORF.
* **RNA2** encodes the capsid precursor **protein alpha**, which after particle assembly
  **autocatalytically cleaves between Asn363 and Ala364** to give mature **beta** (363 aa) and
  **gamma** (44 aa). The question posed in the brief is **confirmed exactly as stated**:
  * PMID 1404613 — "Maturation cleavage required for infectivity of a nodavirus." (J Virol 1992),
    verbatim: *"provirions, which mature by autocatalytic cleavage of the 407-residue coat
    precursor protein between asparagine residue 363 and alanine residue 364"*; N363D/T/A
    abolish cleavage and infectivity, all nine revertants restored Asn.
  * PMID 9765417 — "Specific encapsidation of nodavirus RNAs is mediated through the C terminus
    of capsid precursor protein alpha." (J Virol 1998), verbatim: *"the alpha subunits
    autocatalytically cleave between Asn363 and Ala364. This cleavage generates mature coat
    proteins beta (363 residues) and gamma (44 residues)"*.
  * PMID 2116525 — "Structural homology among four nodaviruses as deduced by sequencing and
    X-ray crystallography." (J Mol Biol 1990): cleavage-site region conserved across BBV, FHV,
    Boolarra and Nodamura.
  * **Pariacoto is the documented exception to the Asn/Ala rule**: PMID 10799587 —
    "Characterization and construction of functional cDNA clones of Pariacoto virus, the first
    Alphanodavirus isolated outside Australasia." (J Virol 2000), verbatim: *"Although the PaV
    capsid precursor is cleaved into the two mature capsid proteins beta and gamma, the amino
    acid sequence at the cleavage site, which is Asn/Ala in all other alphanodaviruses, is
    Asn/Ser in PaV."*

### Table — Alphanodavirus proteins

| Protein | Symbol | Segment / mRNA | Length aa (strain, accession) | Mechanism / frame | PMIDs (verified titles) | Evidence |
|---|---|---|---|---|---|---|
| RNA-dependent RNA polymerase (replicase) | **A** | RNA1 (genomic mRNA) | FHV 998 (NC_004146.1 CDS 40..3036); BBV 998 (NC_001411.2); Boolarra 998 (NC_004142.1); NoV 1043 (NC_002690.1 CDS 22..3153); PaV **973** (NC_003691.1) | single ORF, no cleavage products known; N-proximal transmembrane + mitochondrial targeting signal | PaV 973 aa from 10799587 "Characterization and construction of functional cDNA clones of Pariacoto virus…"; BBV protein A from 3839022 "Structure of the black beetle virus genome…"; targeting from 12208963 "Flock house virus RNA polymerase is a transmembrane protein with amino-terminal sequences sufficient for mitochondrial localization and membrane insertion."; 11689648 "Flock house virus RNA replicates on outer mitochondrial membranes in Drosophila cells."; 7815535 "Requirements for the self-directed replication of flock house virus RNA 1." | **EXP** (function, targeting); FHV/BBV/BoV/NoV lengths **GB** (literature gives Da, not aa) |
| FHV protein A N-terminal mitochondrial localization signal | — | within A | aa **1–46** sufficient | membrane insertion, N-terminus in IMS/matrix, C-terminus cytoplasmic | 12208963 (title above), verbatim: *"the N-terminal 46 amino acids of protein A were sufficient for mitochondrial localization and membrane insertion"* | **EXP** |
| Non-structural protein B1 | **B1** | sgRNA3, **in frame with A**, internal AUG; = C-terminal fragment of A | FHV **102** (NC_004146.1 CDS 2728..3036); BBV 102 (NC_001411.2); NoV ~**136** (A residues 908–1043) | internal initiation on sgRNA3; *not* a cleavage product | 15163762 "Nodamura virus nonstructural protein B2 can enhance viral RNA accumulation in both mammalian and insect cells."; 12573589 "Recovery of infectivity from cDNA clones of nodamura virus and identification of small nonstructural proteins." | NoV B1 protein detected **EXP** (SDS-PAGE + MALDI-TOF, 12573589); FHV/BBV lengths **GB**; frame relationship **EXP/SEQ** |
| Non-structural protein B2 (RNAi suppressor) | **B2** | sgRNA3, **+1 frame** relative to A | FHV **106**; BBV 106; Boolarra 106; NoV **137 and 134** (two forms); PaV **90** | +1-frame overlapping ORF; NoV's two forms arise from the **first and second AUG of RNA3** (alternative initiation), differing only at the N terminus | 12016316 "Induction and suppression of RNA silencing by an animal virus." (Science 2002); 15163762 (title above); 12573589 (title above) — *"showed synthesis of B1 and two forms of B2 (B2-134 and B2-137)"*; 16228003 "Dual modes of RNA-silencing suppression by Flock House virus protein B2."; 16270100 "The structure of the flock house virus B2 protein, a viral suppressor of RNA interference, shows a novel mode of double-stranded RNA recognition."; 19249868 "Structure of the RNA-binding domain of Nodamura virus protein B2, a suppressor of RNA interference."; 24115437 "RNA interference functions as an antiviral immunity mechanism in mammals."; PaV 90 aa from 10799587 | **EXP** (NoV B2-134/B2-137 by MS; FHV B2 106 aa stated in 15163762; function by structure + genetics). PaV B2 **SEQ** (ORF described, *"could encode a protein of 90 amino acids"*) |
| Capsid precursor | **alpha** (α) | RNA2 | FHV **407** (NC_004144.1); BBV **407** (PMID 6548308); Boolarra 403 (NC_004145.1); NoV 399 (NC_002691.1); PaV **401** (NC_003692.1) | precursor; cleaved post-assembly | 6548308 "Primary and secondary structure of black beetle virus RNA2, the genomic messenger for BBV coat protein precursor." — *"an open reading frame for a protein 407 amino acids long"*; 1404613; 9765417; 10799587; 2798110 "Nucleotide sequences of three Nodavirus RNA2's: the messengers for their coat protein precursors." (no abstract — sequence source only) | **EXP** |
| Mature capsid protein | **beta** (β) | RNA2, from alpha | FHV **363** (alpha 1–363) | autocatalytic cleavage at **Asn363 ↓ Ala364** (Asn363 ↓ **Ser** in PaV) | 1404613; 9765417; 10799587 | **EXP** |
| Mature gamma peptide | **gamma** (γ) | RNA2, from alpha | FHV **44** (alpha 364–407) | same cleavage; pentameric helical bundle at 5-fold axes; Phe402/405/407 required for specific RNA packaging | 9765417; 18077727 "Rescue of maturation-defective flock house virus infectivity with noninfectious, mature, viruslike particles."; 16378964 "Morphological changes in the T=3 capsid of Flock House virus during cell entry." | **EXP** |
| Minor capsid protein (N-truncated alpha) | **alpha'** (no standard symbol) | RNA2 | PaV **377** (= alpha minus 24 N-terminal residues) | **alternative initiation at the second AUG** of the capsid ORF; assembles and undergoes maturational cleavage; dispensable for infectivity | 13679619 "Virions of Pariacoto virus contain a minor protein translated from the second AUG codon of the capsid protein open reading frame." | **EXP** |
| Second ORF on RNA2 ("URF") | — | RNA2 | BBV **72** (NC_002037.1 CDS 1110..1328) | predicted ORF only, no product reported | 6548308 — *"A second open reading frame for a putative protein 72 amino acid residues long begins at base 1110."* | **computational prediction only** — do **not** build a feature |

### Lineage notes

* **Newington virus**: PubMed search returns **no** papers (`esearch` for "Newington virus
  nodavirus" → no hits) and NCBI holds **no RefSeq** for it. Any annotation for Newington rests
  on GenBank deposits alone — flag as **GB**.
* **Boolarra virus** description: PMID 4044199 — "Boolarra virus: a member of the Nodaviridae
  isolated from Oncopera intricoides (Lepidoptera: Hepialidae)." No B1 CDS is annotated on
  NC_004142.1 (only protein A and B2).
* **Pariacoto** structure: PMID 11135676 — "The structure of pariacoto virus reveals a
  dodecahedral cage of duplex RNA."
* **Nodamura virus** is the only alphanodavirus that also infects mammals; B2 is the RNAi
  suppressor in mice (PMID 24115437, title above).
* **Protein B is dispensable for RNA1 replication**: PMID 7815535 — *"protein B could be
  eliminated without inhibiting replication by mutations that disrupted the second of the two
  overlapping open reading frames on RNA 3"*. (Note PMID 7884931 is the *erratum* to this paper,
  not the paper.)

---

## 2. Betanodavirus (fish nervous necrosis viruses)

Genotypes: SJNNV, RGNNV, BFNNV, TPNNV, plus AHNV (Atlantic halibut, BFNNV-type) and GGNNV/OSGNNV
(RGNNV-type).

**Premise correction:** GenBank **D30814 is SJNNV RNA2** (coat protein gene, strain SJ-G91,
1410 bp), *not* RNA1. SJNNV RNA1 is AB025018 / NC_003448.

### Table — Betanodavirus proteins

| Protein | Symbol | Segment / mRNA | Length aa (strain, accession) | Mechanism / frame | PMIDs (verified titles) | Evidence |
|---|---|---|---|---|---|---|
| RNA-dependent RNA polymerase | **A** (protein A) | RNA1 | SJNNV **983** (SJNag93; NC_003448.1 CDS 79..3030); RGNNV **982** (SGWak97 NC_008040.1; also *E. akaara* isolate, 982 aa); GGNNV **982** (Singapore strain); AHNV **981**; BFNNV 981 (NC_013458.1); TPNNV 982 (NC_013460.1) | single ORF; **no cleavage products, no alternative forms documented** | 10580064 "Sequence of the non-structural protein gene encoded by RNA1 of striped jack nervous necrosis virus."; 11172107 "Determination of the complete nucleotide sequences of RNA1 and RNA2 from greasy grouper (Epinephelus tauvina) nervous necrosis virus, Singapore strain."; 15109133 "Complete sequence of RNA1 and subgenomic RNA3 of Atlantic halibut nodavirus (AHNV)." — *"a main open reading frame encoding a polypeptide of 981 amino acids"*; 22270757 "Complete sequence of a viral nervous necrosis virus (NNV) isolated from red-spotted grouper (Epinephelus akaara) in China." — *"a 982-amino-acid open reading frame … encoding protein A"*; 18717064 "Comparisons among the complete genomes of four betanodavirus genotypes."; 11602777 "Establishment of an infectious RNA transcription system for Striped jack nervous necrosis virus, the type species of the betanodaviruses." | **EXP/SEQ** (lengths literature-stated for SJNNV, RGNNV, GGNNV, AHNV; BFNNV/TPNNV **GB**) |
| Protein A mitochondrial targeting signals | — | within A | GGNNV MTS **aa 215–255**; AHNV **aa 1–40** *and* **aa 225–246** | membrane association, mitochondrial replication site | 15163743 "Membrane association of greasy grouper nervous necrosis virus protein A and characterization of its mitochondrial localization targeting signal."; 17602779 "RNA-dependent RNA polymerase from Atlantic halibut nodavirus contains two signals for localization to the mitochondria." (both **abstract only** for the exact residue ranges) | **EXP** |
| Non-structural protein B1 | **B1** | **sgRNA3**, **same frame as protein A**, internal AUG, shares the protein A stop | RGNNV **111** (gene 336 nt); SJNNV **112** (computed from NC_003448.1 in-frame AUG at nt 2692 — *not* a published coordinate) | internal initiation on sgRNA3; **in-frame C-terminal fragment of protein A**; not a cleavage product | 19136133 "Betanodavirus non-structural protein B1: A novel anti-necrotic death factor that modulates cell death in early replication cycle in fish cells." — *"The B1 gene (336 nt) was cloned from the redspotted grouper nervous necrosis virus (RGNNV) genome"*; 36560822 "Characterization of Nervous Necrosis Virus (NNV) Nonstructural Protein B2 and Its Enhancement on Virus Proliferation." — *"The B1 protein consists of 336 nucleotides and its open reading frame (ORF) matches the 3′ termini of the RdRp reading frame"*; 29449573 "Nuclear targeting of the betanodavirus B1 protein via two arginine-rich domains induces G1/S cell cycle arrest mediated by upregulation of p53/p21."; 36740084 "Red-grouper nervous necrosis virus B1 protein inhibits fish IFN response by targeting Ser5-phosphorylated RNA polymerase II to promote viral replication."; 42829137 "Construction of red-spotted grouper nervous necrosis virus (RGNNV) nonstructural protein B1/B2 mutant strains and dynamics during infection."; 32050492 "Betanodavirus and VER Disease: A 30-year Research Review." (notes *"there has been some controversy about the existence of B1"*) | **EXP** (protein expressed, knocked down, functional); **start coordinate is not published** — no RefSeq annotates a B1 CDS |
| Non-structural protein B2 (RNAi suppressor; necrosis factor) | **B2** | **sgRNA3**, **+1 frame** relative to protein A | **75** in SJNNV (NC_003448.1 CDS 2756..2983, annotated "protein B"), RGNNV (NC_008040.1 2753..2980), AHNV, BFNNV, TPNNV | +1-frame overlapping ORF; sgRNA3 is capped and begins at **RNA1 nt 2730**, made by **premature termination**, not splicing | 16186236 "Characterization of Striped jack nervous necrosis virus subgenomic RNA3 and biological activities of its encoded protein B2."; 15109133 (AHNV: *"an open reading frame encoding a peptide of 75 amino acids corresponding to a hypothetical B2 protein"*); 16352533 "Betanodavirus B2 is an RNA interference antagonist that facilitates intracellular viral RNA accumulation."; 16809288 "Sequestration and protection of double-stranded RNA by the betanodavirus b2 protein."; 17669362 "Identification of critical residues in nervous necrosis virus B2 for dsRNA-binding and RNAi-inhibiting activity through by bioinformatic analysis and mutagenesis." (abstract only); 19116179 "Betanodavirus non-structural protein B2: A novel necrotic death factor that induces mitochondria-mediated cell death in fish cells."; 20870718 "Betanodavirus B2 causes ATP depletion-induced cell death via mitochondrial targeting and complex II inhibition in vitro and in vivo."; 25008790 "Beta-nodavirus B2 protein induces hydrogen peroxide production, leading to Drp1-recruited mitochondrial fragmentation and cell death via mitochondrial targeting."; 36560822 (*"its ORF is in a + 1 reading frame relative to RdRp"*, *"consists of 75 amino acids"*); 41146288 "Giant grouper nervous necrosis virus subgenomic RNA3 is transcribed by a premature termination mechanism via intragenomic RNA-RNA pairing."; 42218873 "Nervous necrosis virus B2 non-structural protein shut off host translation via sequestration of polyadenylate-binding protein." (abstract only) | **EXP** |
| Capsid / coat protein — **NOT cleaved** | **CP** (sometimes written "protein alpha" by analogy only) | RNA2 | SJNNV **340** (SJ-G91, **D30814**; NC_003449.1); RGNNV **338** (NC_008041.1; also *E. akaara* isolate); GGNNV **338**; TPNNV 340 (NC_013461.1); BFNNV 338 (NC_013459.1); OSGNNV 338 (KT071606, PDB 4WIZ) | single uncleaved chain. N-ARM 1–33 (motif 23RRRANNRRRSN33), N-arm 34–51, S domain 52–213, linker 214–220, P domain 221–338 — **structural regions, not mat_peptides** | 9049363 "Comparison of the coat protein genes of five fish nodaviruses, the causative agents of viral nervous necrosis in marine fish." — *"a single ORF of 1023 bases encoding a protein of 340 amino acids"*; 22270757 (*"a 338-amino-acid major ORF encoding coat protein"*); 11172107; 26491970 "Crystal Structures of a Piscine Betanodavirus: Mechanisms of Capsid Assembly and Viral Infection."; 12642096 "Identification of nucleolus localization signal of betanodavirus GGNNV protein alpha."; 12185275 "Intramolecular disulfide bonding is essential for betanodavirus coat protein conformation."; 39087906 "Molecular Mechanism of pH-Induced Protrusion Configuration Switching in Piscine Betanodavirus Implies a Novel Antiviral Strategy."; 34133901 "Nervous Necrosis Virus Coat Protein Mediates Host Translation Shutoff through Nuclear Translocalization and Degradation of Polyadenylate Binding Protein." | **EXP** |

### The capsid-cleavage question, settled

**Betanodavirus CP is not proteolytically cleaved; there is no alpha → beta + gamma maturation.**
The primary experimental paper is PMID 9672600 — "A fish encephalitis virus that differs from
other nodaviruses by its capsid protein processing." (Arch Virol 1997). Abstract, verbatim: the
RNA2 ORF *"carries the catalytic D-75 residue but lacks the site for autocatalytic proteolysis,
the process yielding the two capsid proteins of insect nodaviruses"*; the observed 43–45 kDa
doublet *"did not result from the assembly-dependent autocatalytic cleavage of a protein
precursor"* and was also shown **not** to come from alternative initiation, a −1 frameshift, or
degradation — the two forms are co-synthesised and attributed to co-translational modification.
(abstract only; the abstract carries the full result.)

Corroboration: PMID 32050492 (review, full text) — *"In Alphanodavirus … protein α, is
auto-catalytically cleaved into proteins β and γ … This mechanism was not observed in
betanodavirus"*; and PMID 26491970, whose T=3 structure is built from the **intact 338-residue
CP** and whose T=3 VLPs assemble from a single uncleaved polypeptide expressed in *E. coli*.

**Module consequence:** do not define beta/gamma mat_peptides for betanodaviruses. Do not let
the alphanodavirus Asn/Ala rule propagate across the family.

### Other betanodavirus ORFs

None beyond protein A, B1, B2 and CP has experimental support. RNA2 carries a single ORF
(PMIDs 9672600, 9049363, 11172107). RefSeqs NC_003448/NC_003449/NC_008040/NC_008041 annotate
only protein A, "protein B" (= B2, SJNNV only) and coat protein — **no B1 CDS is annotated
anywhere**, so a B1 feature must be added by hand.

---

## 3. Crustacean nodaviruses (the proposed but unratified "Gammanodavirus")

MrNV, PvNV, CMNV. (XSV is a satellite virus in a *different family* — see §3.4.)

### Table — crustacean nodavirus proteins

| Protein | Symbol | Segment | Length aa (strain, accession) | Mechanism / frame | PMIDs (verified titles) | Evidence |
|---|---|---|---|---|---|---|
| RdRp / protein A | **A** / RdRp | RNA1 | MrNV **1045** (NC_005094.1 CDS 24..3161) — GenBank-only; literature gives only a **partial** 1034 aa; PvNV **1033** (NC_014978.1) — GenBank-only; CMNV **1043** (MT270124.1) — literature-backed | single ORF; methyltransferase + RdRp separated by a predicted TM region (model only) | 23318596 "Genetic analysis of RNA1 and RNA2 of Macrobrachium rosenbergii nodavirus (MrNV) isolated from India." (partial RdRp 1034 aa); 32659308 "Genomic characterization of covert mortality nodavirus from farming shrimp: Evidence for a new species within the family Nodaviridae." — *"the RdRp of length 1043 amino acids"*; 21790065 "Ultrastructural and sequence characterization of Penaeus vannamei nodavirus (PvNV) from Belize." (ORFs described as *"may encode"*); 36757137 "Structural and functional characterization of RNA dependent RNA polymerase of Macrobrachium rosenbergii nodavirus (MnRdRp)." (AlphaFold model, abstract only) | MrNV/PvNV **GB**; CMNV **SEQ** |
| B2 protein | **B2** | RNA1, overlapping 3' end of protein A in a shifted frame | MrNV **133** (NC_005094.1 CDS 2725..3126); PvNV **131** (NC_014978.1, "B2-like protein"); CMNV **132** (MT270124.1) | overlapping-frame ORF. **No subgenomic RNA3 has been demonstrated for any of these three viruses** — do not claim one | 23318596 — *"another a B2-like protein with a length 133 amino acids"*; 32659308 — *"the protein B2 with a length of 132 amino acids"*; 21790065 (PvNV, predicted); 25982399 "Interaction study of a novel Macrobrachium rosenbergii effector caspase with B2 and capsid proteins of M. rosenbergii nodavirus reveals their roles in apoptosis." (abstract only) | **SEQ** for the ORFs. **RNAi-suppressor activity has NOT been demonstrated for MrNV, PvNV or CMNV B2** — it is generic-nodavirus inference. Word the product as "B2 protein", not "B2 RNAi suppressor" |
| B1 protein | — | — | **not documented** | — | no paper found; not annotated in NC_005094.1, NC_014978.1 or MT270124.1 | **absent / undocumented** |
| Capsid protein — **not cleaved** | **CP** | RNA2 | MrNV **371** (NC_005095.1; literature-backed twice); PvNV **368** (NC_014977.1) GenBank-only; CMNV **437** (MT270123.1) literature-backed | single uncleaved polypeptide; **no gamma peptide**. Structure: NTA/N-arm 1–~30 disordered (RNA-binding region residues 20–29), S domain 62–242, linker 241–258, **dimeric** P-domain spikes (unlike the trimeric spikes of alpha-/betanodaviruses); PDB 6H2B | 23318596 — *"a single ORF encoding a capsid protein of 371 amino acids"*; 30346944 "Structure of the Macrobrachium rosenbergii nodavirus: A new genus within the Nodaviridae?" — *"The full-length Mr_NV CP is a polypeptide of 371 amino acids"*, *"Betanodaviruses and Mr_NV both have their C-termini at the capsid exterior and do not encode a γ-peptide"*; 28522842 "Cryo-Electron Microscopy Structure of the Macrobrachium rosenbergii Nodavirus Capsid at 7 Angstroms Resolution."; 15660790 "White tail disease of the giant freshwater prawn, Macrobrachium rosenbergii: separation of the associated virions and characterization of MrNV as a new type of nodavirus." — virion capsid is *"a single polypeptide of 43 kDa"*; 24878641 "RNA-binding region of Macrobrachium rosenbergii nodavirus capsid protein."; 32659308 — *"RNA2 is 1448 bp and encodes a capsid protein of 437 amino acids"* | **EXP** (MrNV); CMNV **SEQ**; PvNV **GB** |

### The one reported MrNV capsid "cleavage" is an expression artefact

PMID 34445426 — "Regulation of Proteolytic Activity to Improve the Recovery of Macrobrachium
rosenbergii Nodavirus Capsid Protein." (Int J Mol Sci 2021): Edman sequencing located a cut at
**Arg26** of recombinant MrNV CP, caused by an *E. coli* **cysteine protease** (blocked by
E-64). This is a heterologous-expression artefact, **not** viral maturation. Do not model it as
a mat_peptide. Likewise the N-terminal and C-terminal truncations in PMIDs 24878641, 27693291
("C-terminal domain on the outer surface of the Macrobrachium rosenbergii nodavirus capsid is
required for Sf9 cell binding and internalization.") and 30707739 ("Expression, purification and
characterization of the dimeric protruding domain of Macrobrachium rosenbergii nodavirus capsid
protein expressed in Escherichia coli.") are designed constructs.

### 3.4 XSV — satellite virus, *Sarthroviridae*, not Nodaviridae

| Protein | Symbol | Genome | Length aa | Mechanism | PMIDs | Evidence |
|---|---|---|---|---|---|---|
| Capsid protein 17 | **CP17** | monocistronic ssRNA, 796 nt (NC_043197.1) / 872 nt (NC_043492.1) | **174** (NC_043492.1 CDS 89..613) | — | 14993649 "Characteristics of the monocistronic genome of extra small virus, a virus-like particle associated with Macrobrachium rosenbergii nodavirus: possible candidate for a new species of satellite virus." — ssRNA 796 nt encoding *"a single structural protein, capsid CP-17"*; 17928085 "Cloning and sequencing of capsid protein of Indian isolate of extra small virus…" | **EXP** for CP17; the 174 aa figure is **GB** |
| Capsid protein 16 | **CP16** | same ORF | **163** (NC_043492.1 CDS 122..613, same stop) | **alternative in-frame downstream initiation** (not cleavage) | 15660790 — virion *"capsid exhibits two polypeptides of 16 and 17 kDa"*; the two-CDS layout itself is **GB** (NC_043492.1) plus the ICTV Sarthroviridae chapter; the primary sequence paper is in a non-PubMed Chinese-language journal | two polypeptides in virions **EXP**; the dual-initiation mechanism **GB / ICTV chapter** |

Taxonomy: PMID 30507372 — "ICTV virus taxonomy profile: Sarthroviridae." — *"a single genus,
Macronovirus … a single species, Macrobrachium satellite virus 1 … satellite viruses of
Macrobrachium rosenbergii nodavirus, an unclassified virus related to members of the family
Nodaviridae."* XSV therefore belongs in a separate module, not the Nodaviridae one.

Also useful: PMID 31734381 — "In experimental challenge with infectious clones of Macrobrachium
rosenbergii nodavirus (MrNV) and extra small virus (XSV), MrNV alone can cause mortality in
freshwater prawn (Macrobrachium rosenbergii)."

---

## 4. Caenorhabditis nodaviruses — Orsay, Santeuil, Le Blanc, Mělník

These are the only nodavirus-like viruses with a **documented programmed ribosomal frameshift**
and **non-AUG initiation**, and the only ones that **lack B1/B2 entirely**.

### Table — Caenorhabditis nodavirus proteins

| Protein | Symbol | Segment | Length aa (strain, accession) | Mechanism | PMIDs (verified titles) | Evidence |
|---|---|---|---|---|---|---|
| RdRp / protein A | **A** / RdRp | RNA1 | Orsay **937** (JU1580, NC_028097.1 = HM030970, CDS 527..3340); Santeuil **982** (NC_015069.1); Le Blanc **942** (NC_028133.1); Mělník **≥745 partial** (MK774657.1/MK774658.1) | ordinary AUG start, single ORF. **RNA1 has only this one CDS — no sgRNA3, no B1, no B2** | 21283608 "Natural and experimental infection of Caenorhabditis nematodes by novel viruses related to nodaviruses."; 23043172 "Complete genome sequence of Le Blanc virus, a third Caenorhabditis nematode-infecting virus."; 24503084 "Orsay virus utilizes ribosomal frameshifting to express a novel protein that is incorporated into virions."; 31434736 "Noda-Like RNA Viruses Infecting Caenorhabditis Nematodes: Sympatry, Diversity, and Reassortment." | **SEQ** (lengths **GB**/RefSeq; genome papers describe the ORF) |
| Capsid protein (alpha) | **CP** / α | RNA2 | Orsay **391** (NC_028098.1 CDS 89..1264); Santeuil **389** (NC_015070.1, N-term NKNNKK); Le Blanc **388** (NC_028134.1, N-term NKHRTH); Mělník ≥375 partial | **non-AUG initiation**: mature capsid N-terminus is **NKNNT**, corresponding to an **AAC codon at nt 89–91**, with **no AUG 5′ of it** in RNA2. RefSeq encodes this as `/transl_except=(pos:89..91,aa:Asn)` | 24503084, verbatim: *"the N-terminal sequence was NKNNT… The asparagine (Asn) corresponds to nucleotides 89–91, an AAC codon, which is not traditionally thought to be a suitable translation initiation site. In the complete RNA2 segment, there was no AUG codon 5′ to the codon corresponding to the N-terminal Asn amino acid."*; 25136116 "Crystal structure of a nematode-infecting virus." | **EXP** (Edman) for Orsay; Santeuil/Le Blanc **SEQ** — the Le Blanc RefSeq says the product is *"predicted based on analogy with Orsay virus"* |
| CP-δ fusion protein | **CP-δ** | RNA2 | Orsay **765** (NC_028098.1 `join(89..1261,1261..2385)`); Santeuil **765**; Le Blanc **827**; Mělník ≥751 partial | **−1 programmed ribosomal frameshift**. Slippery heptamer **G_UCA_AAC = nt 1255–1261** of Orsay RNA2, immediately 5′ of the alpha-ORF stop (UAG, nt 1262–1264); a stimulatory **stem-loop begins ~8 nt downstream**; **no pseudoknot** detected | 24503084, verbatim: *"the presence of a putative heptameric slippery sequence, G_UCA_AAC … immediately before the alpha ORF stop codon"*, *"a putative RNA stem loop structure beginning eight nucleotides downstream of the slippery sequence, consistent with a canonical eukaryotic −1 frameshift cassette"*, *"no evidence for a pseudoknot"*, and the transframe peptide was MS-sequenced; 28241071 "Structure of a pentameric virion-associated fiber with a potential role in Orsay virus entry to host cells."; 32817218 "Orsay Virus CP-δ Adopts a Novel β-Bracelet Structural Fold and Incorporates into Virions as a Head Fiber." | **EXP** (MS/MS transframe peptide, gradient co-purification) for Orsay; Le Blanc/Santeuil/Mělník **SEQ** (Le Blanc shares the slippery sequence; Santeuil's differs) |
| Protein delta (free) | **δ** | RNA2, own AUG | Orsay **346** (CDS 1345..2385, ADW54430.1, **HM030971.2 only — dropped from NC_028098.1**); Santeuil **332** (HM030973); Le Blanc 398 (JQ943580.1, superseded coords) | **separate AUG initiation** inside the δ frame; δ = CP-δ residues **420–765** | 21283608 — Santeuil RNA2 *"appeared to encode a second ORF of 332 amino acids at the 3′ end… named ORF δ"*; 29743360 "Orsay δ Protein Is Required for Nonlytic Viral Egress." — an **ATG→CTG mutant** that makes normal CP-δ but no free δ has a viral-exit defect (**abstract only**); 32817218; 28241071 | **EXP** (reverse genetics; the 2014 paper had been cautious, the 2018 paper settled it) |

### Caenorhabditis-specific cautions

* **No B2.** PMID 21283608, verbatim: *"All known nodavirus B2 proteins overlap with the
  C-terminus of the RNA-dependent RNA polymerase and are encoded in the +1 frame relative to the
  polymerase. No open reading frame with these properties was predicted in the 3′ end of the RNA1
  segment."* PMID 24503084: *"they appear to lack B1/B2 proteins."* Functional corroboration:
  PMID 23468484 "Characterization of virus-encoded RNA interference suppressors in Caenorhabditis
  elegans." — *"No VSR activity was detected for either of the two Orsay viral proteins proposed
  previously as VSRs."*
* **The non-AUG mechanism is unresolved in the literature.** The 2014 paper floats an upstream
  **UUG at nt 56–58** followed by proteolysis, then argues against it (*"the UUG codon in Orsay
  virus lacks a strong initiation context (ggc_UUG_c)"*; Santeuil has an in-frame stop five
  codons upstream), and prefers initiation directly at the NKNNT site *"via some unknown
  mechanism"*. Annotate the **fact** (Edman-verified N-terminus, non-AUG start), not a mechanism.
* **Nesting:** Orsay CP = CP-δ residues 1–391; free δ = CP-δ residues 420–765. Three overlapping
  products from one 2.3 kb region.
* **Annotation traps:** Santeuil RefSeq NC_015070.1 labels the *fusion* `/product="protein delta"`
  and gives the *capsid* `/note="protein A"`. Orsay RefSeq NC_028098.1 **drops the free-δ CDS**
  that HM030971.2 carries. Mělník has **no complete CDS and no RefSeq** (MK774657–MK774660,
  2 isolates, both ends open).
* Tropism and review: PMID 24314656 "Orsay, Santeuil and Le Blanc viruses primarily infect
  intestinal cells in Caenorhabditis nematodes."; PMID 31424970 "Natural Viruses of Caenorhabditis
  Nematodes."

---

## 5. Other noda-like viruses (brief)

| Virus | Notable features | Lengths / accessions | PMIDs (verified titles) | Evidence |
|---|---|---|---|---|
| **Mosinovirus (MoNV)** | The most unusual member: **two subgenomic RNAs** (RNA3 580 nt, RNA4 292 nt), with **B2 encoded on RNA4**; capsid gene acquired from outside the family (16% identity to Lake Sinai virus 2); virions ~50 nm; mature capsid ~56 kDa and *"do not seem to be cleaved from a longer precursor"* | RNA1 2972 nt, RNA2 1801 nt | 25210176 "A unique nodavirus with novel features: mosinovirus expresses two subgenomic RNAs, a capsid gene of unknown origin, and a suppressor of the antiviral RNA interference pathway." | **EXP** |
| **Wuhan nodavirus (WhNV)** | RNA1 3149 nt, three predicted ORFs (A, B1, B2), but **B1 protein is not produced**: *"expression of viral proteins, ProA, B2 and Proalpha, with the absence of B1 production"*. sgRNA3 370 nt, 5′ capped, 3′ co-terminal, **single B2 ORF with the first AUG as the authentic start**. B2 is an RNAi suppressor (aa 1–20 essential) | RNA2 1562 nt, one ORF (protein alpha) | 20441781 "Characterization of Wuhan Nodavirus subgenomic RNA3 and the RNAi inhibition property of its encoded protein B2."; 16780981 "Isolation and RNA1 nucleotide sequence determination of a new insect nodavirus from Pieris rapae larvae in Wuhan city, China."; 24100569 "Crystallization and preliminary crystallographic analysis of a viral RNA-silencing suppressor encoded by Wuhan nodavirus." | **EXP** |
| **Culannivirus / CulV** (*Culex annulirostris* virus, Australia) | RdRp phylogeny resolves **six clades**, only two matching existing genera; the paper argues more genera are needed but **proposes no names** | RNA1 MH794142.1 (3030 nt, RdRp 946 aa); RNA2 MH794143.1 (1385 nt, capsid precursor alpha 406 aa) — both **GB** | 30596795 "The taxonomy of an Australian nodavirus isolated from mosquitoes." | ORF lengths **GB**; phylogeny **EXP/SEQ** |
| **Culex mosquito virus 6** | **No verifiable PubMed paper found.** Do not cite one. | — | — | — |

---

## Settled questions (a)–(g)

**(a) Is the betanodavirus capsid proteolytically cleaved? — NO.**
There is no alpha→beta+gamma maturation cleavage in betanodaviruses. The RNA2 ORF *"carries the
catalytic D-75 residue but lacks the site for autocatalytic proteolysis"* (PMID 9672600, "A fish
encephalitis virus that differs from other nodaviruses by its capsid protein processing.",
abstract only but the abstract carries the full result, which also excludes alternative
initiation and −1 frameshifting as explanations for the observed 43/45 kDa doublet). Confirmed by
the review PMID 32050492 and by the T=3 crystal structure built from the intact 338-residue CP
(PMID 26491970). **Annotate CP as one chain, 338 aa (340 aa in SJNNV and TPNNV). No mat_peptides.**

**(b) Do betanodaviruses and MrNV encode B1 as well as B2?**
*Betanodaviruses:* **yes, both.** B2 = **75 aa**, +1 frame relative to protein A, on sgRNA3
(PMIDs 16186236, 15109133, 36560822). B1 = **111 aa** in RGNNV (gene 336 nt, PMID 19136133),
**112 aa** in SJNNV (computed, not published), in the **same frame** as protein A — it is the
C-terminal fragment of protein A, translated from an internal AUG on sgRNA3, and shares the
protein A stop. B1 is experimentally real (expressed, knocked down, functionally characterised:
PMIDs 19136133, 29449573, 36740084) but **its start coordinate has never been published and no
RefSeq annotates a B1 CDS** — a B1 feature must be hand-placed. The 2020 review notes
*"there has been some controversy about the existence of B1"*, and the AHNV RNA3 paper
(15109133) and the 2025 GGNNV sgRNA3 paper (41146288) describe RNA3 as carrying only B2.
*MrNV:* **B2 only, 133 aa** (PMID 23318596). **No B1 is documented or annotated.** Also: no
subgenomic RNA has been demonstrated for MrNV — do not assert an sgRNA3. And MrNV B2's
RNAi-suppressor activity has **not** been shown; that is generic-nodavirus inference.

**(c) Orsay/Santeuil/Le Blanc: δ and CP-δ.**
CP-δ is made by a **−1 programmed ribosomal frameshift**. Orsay slippery heptamer
**G_UCA_AAC at RNA2 nt 1255–1261**, immediately 5′ of the alpha-ORF stop codon (UAG, nt
1262–1264); stimulatory stem-loop begins ~8 nt downstream; **no pseudoknot**. Demonstrated by
MS/MS of the transframe peptide (PMID 24503084). RefSeq encodes it as
`join(89..1261,1261..2385)`, nt 1261 used twice — this is a frameshift, **not** a splice.
**A free δ does exist**: translated from its own AUG at RNA2 nt 1345–1347 (346 aa in Orsay;
332 aa in Santeuil), and an ATG→CTG mutant that still makes CP-δ but no free δ has a non-lytic
**egress** defect (PMID 29743360, abstract only). Note the free-δ CDS is present in HM030971.2
but **absent from RefSeq NC_028098.1**.
**The capsid (and hence CP-δ) initiates at a non-AUG codon**: Edman sequencing gives N-terminus
**NKNNT** at the **AAC codon nt 89–91**, with no upstream AUG in RNA2 (PMID 24503084). The
*mechanism* is explicitly unresolved by the authors and no later paper settles it. **Protein A
initiates at a normal AUG.**

**(d) Any documented spliced or RNA-edited nodavirus product? — NO, none, in any lineage.**
Basis, stated positively rather than as absence alone: (i) these are cytoplasmic positive-sense
RNA viruses with no nuclear phase; (ii) the one subgenomic message in the family, RNA3, is made
by **premature termination** of negative-strand synthesis / internal initiation, not splicing —
demonstrated for SJNNV (PMID 16186236), GGNNV (PMID 41146288, via an 8-bp long-distance
intragenomic RNA–RNA helix), and Wuhan nodavirus (PMID 20441781, and PMID 21325414 "Internal
initiation is responsible for synthesis of Wuhan nodavirus subgenomic RNA"); (iii) the only
`join()` in any nodavirus RefSeq is the Orsay/Santeuil/Le Blanc **frameshift** join, which reuses
one nucleotide and must not be treated as an intron; (iv) the only "editing" in the literature is
**host** ADAR-type hyperediting of viral dsRNA replication intermediates, antagonised by B2
(PMID 16809288) — a host process, not a programmed edited viral ORF.
PubMed searches for "nodavirus ribosomal frameshifting", "nodavirus RNA editing splicing" return
no family-wide hits other than the Orsay frameshift paper.

**(e) Does protein A have cleavage products or alternative forms? — NO cleavage, anywhere.**
No nodavirus protein A is proteolytically processed in any lineage. The only "products" sharing
protein A sequence are **B1** (alphanodaviruses, betanodaviruses) — a *separate translation
product* of sgRNA3 from an internal AUG in the same frame, not a cleavage product. Protein A's
documented modifications are localisation signals, not processing: FHV aa **1–46** sufficient for
mitochondrial targeting and membrane insertion (PMID 12208963); GGNNV MTS **aa 215–255**
(PMID 15163743); AHNV **two** signals, aa **1–40** and aa **225–246** (PMID 17602779, abstract
only). Related: PMIDs 11689648, 19244330 ("Mitochondrion-enriched anionic phospholipids
facilitate flock house virus RNA polymerase membrane association."), 17301137 ("Nodavirus RNA
replication protein a induces membrane association of genomic RNA.").

**(f) Additional ORFs with experimental support.**
* **Nodamura B2 exists in two lengths, 137 and 134 aa**, from the **first and second AUG of
  RNA3**, differing only at the N terminus; both detected by SDS-PAGE + MALDI-TOF
  (PMIDs 12573589, 15163762). FHV B2 is 106 aa; Pariacoto B2 is 90 aa (sequence-predicted).
* **Pariacoto makes a minor capsid protein from the second AUG of the capsid ORF**, 24 residues
  shorter than alpha (377 aa), present in wild-type virions, assembly- and cleavage-competent,
  dispensable for infectivity (PMID 13679619). This is the clearest alternative-start product in
  the family after NoV B2.
* **BBV RNA2 carries a predicted 72-aa second ORF** at base 1110 (PMID 6548308; RefSeq
  NC_002037.1 "URF"). **Computational prediction only — no product reported. Do not build it.**
* **XSV CP16/CP17**: two polypeptides in virions (PMID 15660790); the dual-initiation
  explanation is GenBank/ICTV-chapter, not PubMed literature.
* **Mosinovirus has two sgRNAs (RNA3 580 nt, RNA4 292 nt)**, with B2 on RNA4 (PMID 25210176) —
  the only nodavirus with more than one subgenomic message.
* **Wuhan nodavirus B1 is predicted but not produced** (PMID 20441781).
* **CMNV ORFs are all "putative"** — RdRp 1043 aa, B2 132 aa, capsid 437 aa, sequence-derived
  only (PMID 32659308). CMNV has **no RefSeq**; use MT270124.1 / MT270123.1.
* **PvNV ORFs are predicted only** (PMID 21790065: *"2 overlapped open reading frames (ORF A and
  B), which may encode a RdRp and a B2 protein"*).

**(g) ICTV genera in Nodaviridae.**
The ICTV Report chapter is **PMID 30431412 — "ICTV Virus Taxonomy Profile: Nodaviridae."**
(J Gen Virol 2019;100(1):3–4), verbatim: *"The family Nodaviridae includes two genera,
Alphanodavirus and Betanodavirus."* There is **no newer Nodaviridae profile in PubMed**.
The current Master Species List (**MSL41, 2025**) still lists **exactly two genera and nine
species**: *Alphanodavirus* (A. boolarraense, A. flockense, A. heteronychi, A. nodamuraense,
A. pariacotoense) and *Betanodavirus* (B. epinepheli, B. pseudocarangis, B. takifugui,
B. verasperi) — binomial names adopted at MSL39.
**There is no ratified Gammanodavirus, Deltanodavirus or Culannivirus genus.**
* "Gammanodavirus" is an author proposal only: PMID 23318596 suggests the name; echoed by
  PMID 30270534 "Molecular biology of Macrobrachium rosenbergii nodavirus infection in giant
  freshwater prawn."; supported structurally by PMIDs 28522842 and 30346944 (whose title is
  itself a question: *"…A new genus within the Nodaviridae?"*).
* "Culannivirus" appears in NCBI taxonomy as a genus-rank name **under "unclassified
  Nodaviridae"** (taxid 2501412) — i.e. not ICTV-ratified. In the source paper (PMID 30596795)
  it is a *virus* name (culannivirus, CulV), and the paper proposes no genus names.
* **MrNV, PvNV, CMNV and the Caenorhabditis viruses are all unclassified**, listed under
  "Related, unclassified viruses" in the ICTV chapter; PMID 30507372 calls MrNV *"an unclassified
  virus related to members of the family Nodaviridae."*
* **XSV is outside Nodaviridae**: family *Sarthroviridae*, genus *Macronovirus* (PMID 30507372).
* The only other family in order Nodamuvirales is *Sinhaliviridae*.

---

## Machine-readable literature length ranges

Ranges are literature/annotation-supported aa lengths across the named strains of each lineage.
`pmid` is the primary verified citation for the range (null where the range is GenBank-only).

```json
{
  "Alphanodavirus": {
    "protein_A":        {"min": 973,  "max": 1043, "pmid": "10799587", "note": "PaV 973 aa literature-stated; FHV/BBV/Boolarra 998 and NoV 1043 are RefSeq CDS lengths (GenBank-only; the 1985 BBV paper gives 101,873 Da, not aa)"},
    "B1":               {"min": 102,  "max": 136,  "pmid": "15163762", "note": "in-frame C-terminal fragment of protein A from an internal AUG on sgRNA3; FHV 102 (GenBank); NoV B1 starts at protein A residue 908 (paper), giving ~136 aa; not a cleavage product; Pariacoto and Boolarra have no annotated B1"},
    "B2":               {"min": 90,   "max": 137,  "pmid": "15163762", "note": "+1 frame vs protein A on sgRNA3. FHV 106 and NoV 134+137 literature-stated (two NoV forms from the 1st and 2nd AUG of RNA3, PMID 12573589). PaV 90 aa sequence-predicted (PMID 10799587)"},
    "capsid_alpha":     {"min": 399,  "max": 407,  "pmid": "6548308",  "note": "BBV/FHV 407 literature-stated; PaV 401 (PMID 10799587); Boolarra 403 and NoV 399 GenBank-only"},
    "capsid_alpha_short": {"min": 377, "max": 377, "pmid": "13679619", "note": "Pariacoto only: minor virion protein from the SECOND AUG of the capsid ORF, 24 residues shorter than alpha; assembly- and cleavage-competent; dispensable"},
    "capsid_beta":      {"min": 355,  "max": 363,  "pmid": "9765417",  "note": "FHV beta = alpha 1-363, 363 aa, literature-stated; other members inferred by the conserved Asn/Ala (Asn/Ser in PaV) site, so the lower bound is approximate"},
    "capsid_gamma":     {"min": 44,   "max": 48,   "pmid": "9765417",  "note": "FHV gamma = alpha 364-407, 44 aa, literature-stated; PaV gamma is the Asn/Ser variant (PMID 10799587); other members not individually measured"},
    "RNA2_second_ORF":  {"min": 72,   "max": 72,   "pmid": "6548308",  "note": "BBV only, base 1110, COMPUTATIONAL PREDICTION ONLY, no product reported - do not build"}
  },
  "Betanodavirus": {
    "protein_A":        {"min": 981,  "max": 983,  "pmid": "15109133", "note": "SJNNV 983 (PMID 10580064), RGNNV/GGNNV 982 (PMIDs 22270757, 11172107), AHNV 981 (PMID 15109133); BFNNV/TPNNV GenBank-only"},
    "B1":               {"min": 111,  "max": 112,  "pmid": "19136133", "note": "RGNNV gene 336 nt = 111 aa (literature); SJNNV 112 aa computed from the in-frame AUG at RNA1 nt 2692 - NOT a published coordinate. Same frame as protein A, shares its stop. No RefSeq annotates a B1 CDS"},
    "B2":               {"min": 75,   "max": 75,   "pmid": "16186236", "note": "+1 frame vs protein A; sgRNA3 begins at RNA1 nt 2730 by premature termination; 75 aa exactly in SJNNV, RGNNV, AHNV, BFNNV, TPNNV"},
    "capsid_CP":        {"min": 338,  "max": 340,  "pmid": "9049363",  "note": "SJNNV 340 (PMID 9049363) and TPNNV 340; RGNNV/GGNNV/BFNNV 338 (PMIDs 22270757, 11172107). NOT CLEAVED - no beta/gamma, no mat_peptides (PMID 9672600)"}
  },
  "Crustacean_nodavirus_MrNV_PvNV_CMNV": {
    "protein_A":        {"min": 1033, "max": 1045, "pmid": "32659308", "note": "CMNV 1043 literature-stated; MrNV 1045 and PvNV 1033 are GenBank-only (MrNV literature gives a 1034 aa PARTIAL, PMID 23318596)"},
    "B2":               {"min": 131,  "max": 133,  "pmid": "23318596", "note": "MrNV 133 (literature), CMNV 132 (literature, PMID 32659308), PvNV 131 (GenBank-only). Overlapping frame at the 3' end of protein A. NO subgenomic RNA demonstrated; RNAi-suppressor activity NOT demonstrated for any of the three"},
    "B1":               {"min": null, "max": null, "pmid": null,       "note": "not documented and not annotated in any of the three - do not build"},
    "capsid_CP":        {"min": 368,  "max": 437,  "pmid": "30346944", "note": "MrNV 371 (literature x2), PvNV 368 (GenBank-only), CMNV 437 (literature, PMID 32659308). NOT cleaved; no gamma peptide. The Arg26 cut (PMID 34445426) is an E. coli cysteine-protease artefact, not a mat_peptide"}
  },
  "Caenorhabditis_nodavirus_Orsay_Santeuil_LeBlanc_Melnik": {
    "protein_A":        {"min": 937,  "max": 982,  "pmid": "24503084", "note": "Orsay 937, Le Blanc 942, Santeuil 982 (RefSeq CDS lengths); Melnik >=745 partial. Ordinary AUG start. RNA1 has ONE CDS: no sgRNA3, no B1, no B2"},
    "capsid_CP":        {"min": 388,  "max": 391,  "pmid": "24503084", "note": "Orsay 391, Santeuil 389, Le Blanc 388, Melnik >=375 partial. NON-AUG initiation: N-terminus NKNNT at the AAC codon, RNA2 nt 89-91 in Orsay (Edman-verified); first residue is N, not M. Santeuil/Le Blanc predicted by analogy"},
    "CP_delta_fusion":  {"min": 765,  "max": 827,  "pmid": "24503084", "note": "-1 programmed ribosomal frameshift; Orsay slippery heptamer G_UCA_AAC at RNA2 nt 1255-1261, stem-loop ~8 nt downstream, NO pseudoknot. Orsay 765, Santeuil 765, Le Blanc 827, Melnik >=751 partial. RefSeq join(89..1261,1261..2385) is a frameshift, NOT a splice"},
    "delta_free":       {"min": 332,  "max": 398,  "pmid": "29743360", "note": "separate AUG (Orsay RNA2 nt 1345); Orsay 346, Santeuil 332, Le Blanc 398 (superseded coords). delta = CP-delta residues 420-765. Required for non-lytic egress. ABSENT from RefSeq NC_028098.1 - present only in HM030971.2"},
    "B2":               {"min": null, "max": null, "pmid": "21283608", "note": "ABSENT. No +1 ORF at the 3' end of RNA1 (PMID 21283608); no VSR activity for either Orsay protein tested (PMID 23468484). Do not build a B2 feature"}
  },
  "Other_noda_like": {
    "MoNV_protein_A":   {"min": null, "max": null, "pmid": "25210176", "note": "RNA1 2972 nt, RNA2 1801 nt; aa lengths not stated in the abstract"},
    "MoNV_capsid":      {"min": null, "max": null, "pmid": "25210176", "note": "~56 kDa, acquired from outside the family; 'do not seem to be cleaved from a longer precursor'"},
    "MoNV_B2":          {"min": null, "max": null, "pmid": "25210176", "note": "encoded on sgRNA4 (292 nt) - the only nodavirus with TWO subgenomic RNAs (RNA3 580 nt, RNA4 292 nt); binds long dsRNA, blocks Dicer-2"},
    "WhNV_B2":          {"min": null, "max": null, "pmid": "20441781", "note": "sgRNA3 370 nt, single B2 ORF, first AUG authentic, aa 1-20 essential for RNAi suppression; WhNV B1 is PREDICTED but NOT produced"},
    "CulV_protein_A":   {"min": 946,  "max": 946,  "pmid": null,       "note": "MH794142.1 GenBank-only; the paper (PMID 30596795) gives no aa length"},
    "CulV_capsid":      {"min": 406,  "max": 406,  "pmid": null,       "note": "MH794143.1 'putative capsid protein precursor alpha', GenBank-only"},
    "XSV_CP17":         {"min": 174,  "max": 174,  "pmid": "14993649", "note": "SATELLITE VIRUS, family Sarthroviridae - belongs in a separate module. 796 nt monocistronic genome encoding CP-17; the 174 aa figure is GenBank-only (NC_043492.1)"},
    "XSV_CP16":         {"min": 163,  "max": 163,  "pmid": "15660790", "note": "two polypeptides of 16 and 17 kDa in virions is literature-backed; the dual-initiation (same stop, downstream AUG) layout is GenBank/ICTV-chapter only"}
  }
}
```
