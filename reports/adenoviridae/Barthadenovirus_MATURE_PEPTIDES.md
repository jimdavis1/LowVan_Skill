# Barthadenovirus mature peptides — from the measured cleavage sites

**Source: PMID_claude_generated 8659111** (Vrati *et al.*, *Virology* 220:186–199,
1996), read in full. This is **experimental** evidence: proteins and peptides
were purified from OAV287 virions, separated, and **sequenced from the
N-terminus**. It is not a prediction, and it supersedes anything derived from
the human adenovirus consensus.

Every boundary below was re-measured against the module's own collections
before it was written here. Where the paper's measured peptide and the
deposited sequence agree character-for-character, that is stated; where the
paper says a site is *not* used, no peptide is modelled.

---

## 0. The consensus does not travel, and here is the proof

The human adenovirus proteinase recognises (Webster 1989; Anderson 1990, as
cited by the paper):

```
type I    Hy-X-G-X ↓ G
type II   Hy-X-G-G ↓ X          Hy = Leu, Ile or Met
```

Applying that to this genus would have been wrong in **four separate ways**,
each one measured:

| What the human consensus would do | What OAV287 actually does |
|---|---|
| cleave **pX** at the conserved site | **pX is not cleaved at all.** The sites used in Ad2 (res 32) and Ad12 (res 27) are absent from OAV, and the Ad12 res-42 site is **not used** — no C-terminal peptide was recovered from virions, probably because P1′ is Phe. The protein runs at ~10 kDa, uncleaved. |
| cleave **pVI** at its N-terminal type I site | The site is conserved ~30 residues in, but **Phe rather than Ala occupies P1′** at residue 29, so the site **may be masked**. No N-terminal sequence could be obtained. |
| miss **pVII**'s real sites entirely | pVII's ORF *does not predict* a precursor, and there is no site at residue 25 where Ad12 uses `MYGG↓A`. Two **non-conforming** sites are used instead — see §2. |
| never find **p32K** | p32K has no human counterpart at all. |

The paper states the extension directly: in OAV,
**`(N/M/L/I)-X-G-X ↓ G`** and **`(M/L/I)-X-A-X ↓ G`** are *also* recognised.
The `N` at P4 and the `A` at P2 are outside the human specificity.

---

## 1. p32K → mature p32K  ✅ MODEL THIS

**Measured.** N-terminal sequencing of the mature virion protein gave
`AYRKRKKKFQRPKPR` (Table 1, HPLC peak 7), produced by a **type I** cleavage at

```
LAGG ↓ AYR
```

**Verified against the collection.** Genome `114430.16` (OAdV-7) carries
`LAGG↓AYR` at residue 12, and the residues following it are
`AYRKRKKKFQRPKPRK` — the measured sequence exactly.

| | |
|---|---|
| precursor | 286 aa |
| cleavage | after residue 12 |
| mature | **274 aa**, starting `AYRKRKKKFQRPKPRK` |

**Both forms are already in the dump.** Genome `114430` deposits the 286 aa
precursor *and* the 274 aa mature protein, and the shorter is an **exact
C-terminal substring** of the longer. They are currently in one collection,
which would train a single profile on both forms — the VP1/VP2 defect from
Parvoviridae in a new costume.

**The site is not universally present**: only 6 of 20 P32K deposits carry a
type I site in the first 60 residues (`LAGG`, `LMGG`, `LVGG`, `ITGG`, `LMGG`
×2), all at residue 12–13 except one at 34. Per-species treatment is required;
a single motif scan across the genus will not do it.

> **Note on precursor length.** The 1996 paper models this ORF at 177 residues
> (21.5 kDa) with a 129-residue mature product. The **same group's** 1998
> transcript-mapping paper (PMID_claude_generated 9614874, Table 1) gives
> `p32  32.1 kDa` — ~286 aa — which agrees with the deposits (272–355 aa,
> median 304) and with the database. The 1998 figure supersedes the 1996 one.
> What survives both, because it was *measured* rather than modelled, is the
> mature N-terminus `AYRKRK…` and the `LAGG↓AYR` site.

## 2. pVII → peptide 1 + peptide 2 + mature VII  ✅ MODEL THIS

The best-evidenced cleavage in the genus: **three measured peptides, all three
matching the deposited sequence exactly.**

Two cleavages, neither conforming to the human consensus:

```
NTGW ↓ G      (N at P4 — outside human specificity)
MRAT ↓ G      (paper: "does not conform to the specificity previously
               described for human Ad proteinases")
```

OAdV-7 `114430.16`, precursor **111 aa**:

```
MSILVSPSDNTGW GIGTSSMRAT GLKFSKKQPVRVRPYYRAQWGQLNARTSLEKLKTKLKYYEK…
└─ 1–13 ────┘ └ 14–23 ─┘ └─ 24–111, mature VII ───────────────────────
```

| Fragment | Residues | Length | Paper's measured peptide (Table 1) | Match |
|---|---|---|---|---|
| N-terminal peptide | 1–13 | 13 aa | `SILVSPSDN` | contained ✅ |
| peptide II | 14–23 | 10 aa | `GIGTSSMRAT` | **exact** ✅ |
| **mature VII** | 24–111 | **88 aa** | `GLKFSKKQPVRVRPYYR` | **exact** ✅ |

**Conservation across the genus:** `NTGW` sits at position 9 in **16 of 16**
deposits. `MRAT` is at 19–20 in 9 of 16. The N-terminal peptide is the region
the paper calls highly conserved between human, animal and avian adenoviruses,
and the module's own data agrees — every deposit begins `M-[STPA]-[IV]-L-…-SPSDNTG`.

## 3. pVIII → two fragments  ✅ MODEL THIS

**Measured.** A peptide equivalent to Ad12 pVIII peptide III was recovered from
OAV virions: `GSNRLTSYPDIPPILK` (Table 1), confirming that **one of the two
possible sites is used**.

**Located in the collection.** In OAdV-7 (`114430.14`, `114430.16`, 217 aa) that
peptide starts at residue 118, preceded by `SELNGW` — so the site is

```
LNGW ↓ G        (another Hy-X-G-X↓G with W at P1, as in pVII's NTGW↓G)
```

| Fragment | Residues | Length | Paper's predicted size | Check |
|---|---|---|---|---|
| N-terminal | 1–117 | 117 aa | 13.2 kDa | ~12.9 kDa ✅ |
| C-terminal | 118–217 | 100 aa | 11.4 kDa | ~11.1 kDa ✅ |

The paper notes neither product was *detected* on the gel, probably because
they co-migrate with proteins VII and X — but the peptide itself was recovered
and sequenced, so the cleavage is established.

Only the two OAdV-7 deposits contain `GSNRLT` literally; the site must be found
per species by motif, not by string.

## 4. pVI → C-terminal proteinase cofactor only  ⚠️ PARTIAL

- **C-terminus:** a **type II** site generating the **11-residue peptide
  cofactor** that activates the adenovirus proteinase is conserved in OAV (and
  in mouse adenovirus). The OAV peptide is `GVNFNTRRYCY` (cf. Ad2
  `GVQSLKRRRCF`). The penultimate Cys, required for dimerisation and protease
  activation, is conserved, as are 2 of 3 Arg. The paper deduces the active
  consensus as `G-XXXXX-R-R-[yr]-C-[fy]`.
  ⚠️ A truncated Ad2 peptide `KRRRCF` did **not** activate recombinant
  proteinase, so the full 11 residues matter.
- **N-terminus: DO NOT MODEL.** The type I site is conserved ~30 residues in,
  but **Phe occupies P1′ at residue 29** and the paper says the site **may be
  masked**. No N-terminal sequence could be obtained — the N-terminus may be
  blocked.

## 5. Explicitly NOT cleaved — do not model  🚫

| Protein | Why |
|---|---|
| **pX** | Ad2 res-32 and Ad12 res-27 sites are **absent** from OAV. The Ad12 res-42 site, conserved in other adenoviruses, is **not used** — no C-terminal peptide was recovered, probably because P1′ is Phe. Observed at ~10 kDa, consistent with no cleavage. Table 2 lists **no mature size** for protein X. |
| **pIIIa** | In other adenoviruses a type II site `LGGS↓G` cuts ~19 residues from the C-terminus. That site is **not strictly conserved** in OAV. A type I site `LIGH↓G` is present 20 residues from the C-terminus, **but no peptide with the predicted sequence was isolated**, so the paper says plainly it is "unclear whether or not the site is used". Unclear is not a boundary. |
| **protein V** | **Absent from the genome.** |
| **protein IX** | **Absent from the genome.** |

## 6. Other precursor/mature pairs found in the data, not yet in the literature

Scanning every Barthadenovirus collection for a shorter deposit that is an
**exact C-terminal substring** of a longer one *from the same genome* — the
structural signature of a cleaved precursor — turns up six pairs:

| Feature | Genome | Precursor → mature | Removed | Propeptide ends | Status |
|---|---|---|---|---|---|
| **P32K** | 114430 | 286 → 274 | 12 | …NNTA **LAGG** | ✅ confirmed by §1 |
| PIIIA | 130499 | 573 → 516 | 57 | …ARKNETAD | ⚠️ **N-terminal**, so NOT the C-terminal pIIIa site of §5 — a different event, unexplained |
| DBP | 10511 | 381 → 363 | 18 | …KSKSATEE | ⚠️ DBP is not a known cleaved precursor |
| RH2 | 130499 | 186 → 103 | 83 | …VANFNPIS | ⚠️ 83 residues is large; more likely a truncated deposit than a mature form |
| RH4 | 114430 | 157 → 144 | 13 | …VELPTEIW | ⚠️ unexplained |
| E4_3 | 130499 | 225 → 217 | 8 | …MNLVGICK | ⚠️ unexplained |

Only the P32K pair has literature behind it. **The other five are recorded as
observations, not as mature peptides**, and none is being modelled. A shorter
deposit nested at the C-terminus is equally consistent with a submitter
annotating a downstream start codon, and nothing here distinguishes the two.
The RH2 pair is the same 7-sequence length spread already flagged in
`LITERATURE_INVENTORY.md`.

---

## What to build

| Feature | Model | Evidence grade |
|---|---|---|
| `P32K` precursor + `P32K_MAT` (274 aa) | ✅ | measured N-terminal sequence |
| `PVII` precursor + `VII` mature (88 aa) + 2 peptides | ✅ | three exact measured peptides |
| `PVIII` precursor + 2 fragments (117 / 100 aa) | ✅ | measured peptide + size agreement |
| `PVI` C-terminal 11-residue proteinase cofactor | ✅ | conserved, function established |
| pVI N-terminus, pIIIa, pX | 🚫 do not model | paper says masked / unclear / not used |

## Still missing

Siadenovirus has **no equivalent virion-proteomics paper**. Its cleaved
products therefore rest on the human consensus alone, which this file has just
shown does not travel across genera. Until a TAdV-3 proteomics or N-terminal
sequencing study is found, Siadenovirus mature peptides must be labelled as
**inference from the human consensus**, or not built.
