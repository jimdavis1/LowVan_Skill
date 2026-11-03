# Mature peptides — what gets built, per module, and why

Curator's rule, 3 Nov 2026: **"Don't build them unless you can find
literature, or base them on reasonable homology."**

This file records the decision for every candidate cleavage, with the evidence
grade that justifies it. Nothing here is built on the human adenain consensus
alone, because `Barthadenovirus/MATURE_PEPTIDES.md` shows that consensus does
not travel between genera.

Evidence grades used below:

| Grade | Meaning |
|---|---|
| **MEASURED** | the peptide was sequenced from purified virions |
| **HOMOLOGY** | the site is measured in another genus and sits at an unambiguously aligned position here, with the same core motif |
| **NOT BUILT** | neither applies, or the literature says the site is unused |

---

## 1. pVII — BUILD IN ALL FOUR MODULES

The best-evidenced cleavage in the family, and the only one that reaches every
module.

**Measured in two genera, independently, in the same paper.** Table 1 of
PMID_claude_generated 8659111 reports N-terminal sequencing of pVII peptide II
from purified virions of *both* OAV287 (Barthadenovirus) and human Ad5
(Mastadenovirus):

| Virus | Measured peptide II | Immediately follows |
|---|---|---|
| OAV287 | `GIGTSSMRA` | `…SPSDN **TGW** ↓ G` |
| Human Ad5 | `GLRFPSKMF` | `…SPSNN **TGW** ↓ G` |

The paper also states directly that pVII's **N-terminal 16 residues are highly
conserved between viruses of human, animal, and avian origin**, which is the
literature's own warrant for transferring the site.

**Measured in our own collections, all four modules:**

| Module | Consensus N-terminus | Cut after residue | Deposits carrying the site |
|---|---|---|---|
| Mastadenovirus | `MSILISPSNNTGWG` | 13 | **281 / 302** |
| Aviadenovirus | `MSILISPNDNRGWG` | 13 | **38 / 38** |
| Barthadenovirus | `MSILVSPSDNTGWG` | 13 | **15 / 16** |
| Siadenovirus | `MHSVVYSPGDSRGWG` | 14 | **11 / 11** |

The motif is `SP-[AGS]-[DN]-[NS]-[TRS]-G-W ↓ G`. The `GW↓G` core is invariant;
only the position shifts by one in Siadenovirus, and the alignment is
unambiguous.

**Decision:** build `PVII` (CDS, precursor) + `PVII_NTERM` (mat_peptide,
residues 1–13/14) + `PVII_MAT` (mat_peptide, mature protein VII) in all four
modules. **MEASURED** in Mastadenovirus and Barthadenovirus, **HOMOLOGY** in
Aviadenovirus and Siadenovirus.

## 2. pVII second cut (`MRAT↓G`) — BARTHADENOVIRUS ONLY

**MEASURED** in OAV287: mature VII begins `GLKFSKKQPVRVRPYYR` (Table 1, peak 9),
giving an 88 aa mature protein from a 111 aa precursor, with peptide II
occupying residues 14–23.

**NOT transferred.** The paper states this site "does not conform to the
specificity previously described for human Ad proteinases". `MRAT` is present
in 9 of 16 Barthadenovirus deposits and essentially absent elsewhere —
Siadenovirus shows only one partial match (`MRGSG`). Transferring it would be
exactly the mistake this file exists to prevent.

## 3. p32K → mature p32K — BARTHADENOVIRUS ONLY

**MEASURED.** Type I site `LAGG↓AYR`; mature N-terminus `AYRKRKKKFQRPKPR`
(Table 1, peak 7). Confirmed in genome `114430.16` at residue 12.
286 aa precursor → **274 aa** mature.

p32K is **genus-specific**; there is nothing to transfer it to.

**Collection split applied:** genome `114430` deposits both forms and both had
landed in `P32K`. The shorter is now moved to `P32K_MAT` on two signals that
must both hold — it is an exact C-terminal substring of the longer deposit
*from the same genome*, and it does **not** begin with Met, because a cleaved
N-terminus cannot.

## 4. pVIII — BARTHADENOVIRUS ONLY

**MEASURED.** The peptide `GSNRLTSYPDIPPILK` (Table 1, footnoted as equivalent
to Ad12 pVIII peptide III) begins at residue 118 of the 217 aa OAdV-7
precursor, preceded by `SELNGW` — so the site is `LNGW↓G`. Fragments of 117
and 100 aa, against the paper's predicted 13.2 and 11.4 kDa.

**Not transferred to other modules** without checking the aligned position
there first; Mastadenovirus pVIII cleavage is separately documented and will be
taken from its own literature rather than from this site.

## 5. pVI C-terminal proteinase cofactor — BARTHADENOVIRUS

**MEASURED/established.** A type II site near the C-terminus releases the
**11-residue peptide cofactor** that activates the adenovirus proteinase —
`GVNFNTRRYCY` in OAV, `GVQSLKRRRCF` in Ad2. The penultimate Cys required for
dimerisation and protease activation is conserved, as are 2 of 3 Arg. Active
consensus: `G-XXXXX-R-R-[yr]-C-[fy]`.

⚠️ The truncated Ad2 peptide `KRRRCF` does **not** activate recombinant
proteinase, so the boundary matters and the full 11 residues must be taken.

## 6. NOT BUILT — and the reason in each case

| Candidate | Module | Why not |
|---|---|---|
| **pVI N-terminus** | Barthadenovirus | The type I site is conserved ~30 residues in, but **Phe rather than Ala occupies P1′** at residue 29 and the paper says the site **may be masked**. No N-terminal sequence could be obtained. |
| **pIIIa** | Barthadenovirus | The type II site `LGGS↓G` used elsewhere is not strictly conserved; a type I site `LIGH↓G` sits 20 residues from the C-terminus, but **no peptide with the predicted sequence was isolated from purified virus**, so the paper says plainly it is "unclear whether or not the site is used". Unclear is not a boundary. |
| **pX** | Barthadenovirus | **Explicitly not cleaved.** The Ad2 res-32 and Ad12 res-27 sites are absent; the Ad12 res-42 site is **not used** (P1′ is Phe) and no C-terminal fragment was recovered. Table 2 lists no mature size for protein X. |
| **p32K, pVIII, pVI, MRAT cut** | Siadenovirus | No virion-proteomics or N-terminal sequencing study exists for TAdV-3, and these sites are either genus-specific or measurably non-conforming. Only the pVII `GW↓G` cut transfers. |
| PIIIA 573→516, DBP 381→363, RH2 186→103, RH4 157→144, E4_3 225→217 | Barthadenovirus | Structural nesting observed in the data, **no literature**. A shorter deposit nested at the C-terminus is equally consistent with a submitter annotating a downstream start codon. Recorded as observations in `Barthadenovirus/MATURE_PEPTIDES.md`; none built. |

## 7. Mastadenovirus — still to do

Human adenovirus cleavage is far better documented than any of the above, and
the original family-level notes (adenain `(M/I/L)XGX↓G` / `(M/I/L)XGG↓X`,
PMID_claude_generated 12367749, 11683632) **do** apply to Mastadenovirus, which
is the genus they were derived from. The six precursors — pIIIa, pVI, pVII,
pVIII, pX, pTP — are to be built from that literature for this module only,
with each site checked against the collection the same way pVII was above.
The one signal peptide in the family, E3 gp19K, is also Mastadenovirus-only.
