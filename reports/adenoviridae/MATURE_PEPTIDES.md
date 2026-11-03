# Mature peptides and signal peptides — Adenoviridae

The curator's scope is every known protein **and mature peptide, including
signal peptides**. This file records why they cannot come from the dump and
what will be done instead.

## They are not in the data. Measured, not assumed.

- **Zero `mat_peptide` features** in all four reference genomes
  (NC_001405 Human adenovirus C, NC_001720 Fowl adenovirus A,
  NC_000899 Fowl adenovirus D, NC_002513 Bovine adenovirus 2).
- **Zero mature-peptide deposits** among the 21,947 distinct sequences.
  Every deposit named "protein VIII" is 227 aa — the same length as every
  deposit named "pVIII". The database calls the precursor by both names.

This is the `ABSENT` verdict of the skill's step 2b: not a defect in the
triage rules, and not a shortfall to record and move past. It means the
sequences have to be sourced elsewhere.

## Where they come from instead

The adenovirus protease (adenain, L3 23K) cleaves six precursors at
**(M/I/L)XGX↓G** and **(M/I/L)XGG↓X** (PMID_claude_generated 12367749,
"Substrate specificity of adenovirus protease"; DNA-binding and activation,
11683632).

| precursor | products | note |
|---|---|---|
| pIIIa | IIIa | N-terminal cut |
| pVI | pVIn, **VI**, pVIc | cut at **both** ends; the 11-residue pVIc is the protease's own cofactor |
| pVII | VII | |
| pVIII | VIII | |
| pX | **mu** | |
| pTP | iTP, TP | |

So `build_mature.py` takes each precursor collection, finds the consensus
cleavage motif in the aligned column where the literature places it, and cuts
there — rather than inventing a coordinate per sequence.

**The cut must be verified against the motif in each member, not applied by
offset.** A fixed offset from the N-terminus is wrong the moment one species'
precursor is a few residues longer, and nothing downstream would show it.

## The one signal peptide

**E3 gp19K**, Mastadenovirus only. A cleaved N-terminal signal sequence
directing the protein into the ER, with a C-terminal ER-retention signal
(PMID_claude_generated 8350398, 8249282). It is the only signal peptide in the
family's proteome.

SignalP 6.0 is **not installed** on this machine and is academic-licensed, so
it cannot be installed unattended. That does not block anything: the skill's
own decision table says a documented cleavage site beats a prediction, and
gp19K has one. If the curator wants SignalP for corroboration they have to
request the download themselves.

---

## CORRECTION, 3 Nov 2026 — the human consensus does not travel

Everything above derives mature peptides from the published **human**
adenovirus adenain specificity, `(M/I/L)XGX↓G` / `(M/I/L)XGG↓X`
(PMID_claude_generated 12367749, 11683632). Applying that to the whole family
is an assumption, and for Barthadenovirus it is contradicted by the primary
literature.

Vrati *et al.* 1996 (PMID_claude_generated 8659111) analysed proteins and
peptides from **purified OAV287 virions** and identified the cleavage sites
the adenovirus proteinase actually uses in that virus. Some matched sites
known from human adenovirus proteins. Others were new, and the paper states
explicitly that some of them **did not conform to the known specificity of the
human Ad proteinase**.

So for Barthadenovirus:

1. The human consensus **cannot** be used to predict cleavage sites, and any
   mature peptide generated that way would be a fabricated boundary wearing a
   citation that does not support it.
2. The authoritative sites are the measured ones in 8659111, which is **not
   open access** — only the abstract has been read. The sites themselves are
   not yet in hand.
3. The same paper reports a **processed 28 kDa virion polypeptide** from the
   left-hand end, i.e. p32K is itself a precursor (32K → 28K). That is a
   mature peptide this module does not currently model.

**Mature peptides for Barthadenovirus are therefore BLOCKED on obtaining
PMID 8659111**, not merely unbuilt. Siadenovirus has no equivalent
virion-proteomics paper found so far, so its cleaved products rest on the
human consensus alone and must be labelled as the inference they are.
