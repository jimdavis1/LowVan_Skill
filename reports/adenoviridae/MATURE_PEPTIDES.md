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
