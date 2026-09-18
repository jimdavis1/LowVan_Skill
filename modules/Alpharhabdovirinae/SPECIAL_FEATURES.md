# MERHA_L_SPLICED — real, buildable, not yet built

**Status: deferred, deliberately. Do not delete it.**

`Alpharhabdovirinae` declares one `"special": "splice"` feature,
`MERHA_L_SPLICED`. It has no reference set, so it never fires.

## The biology is solid

**Kuwata et al. 2011, J Virol 85:6185 (PMID 21507977)**, "RNA splicing in a new
rhabdovirus from Culex mosquitoes". Culex tritaeniorhynchus rhabdovirus (CTRV)
carries a **76-nt intron** in the L coding region with canonical GU-AG
spliceosomal donor/acceptor sites, a predicted branch point and a
polypyrimidine tract. The boundaries were mapped by strand-specific RT-PCR.
CTRV replicates in the **nucleus**, which is what makes splicing possible at
all — unique among Mononegavirales outside Bornaviridae.

The generator's coordinates are antigenome **8648-8723**, which is exactly 76 nt
and matches the paper. Excision shifts +2 -> +3 and yields one 2123 aa ORF.

CTRV is classified in genus **Merhavirus**, so the `MERHA_` prefix is correct.

## Why it matters more than a missing feature

Without the splice, L is **emitted truncated at the in-frame stop at 8681** by
the module's ordinary `L` profile. A CTRV genome does not merely lose a call —
it gets a **wrong** one. The current state is the worst of both: declared,
unable to fire, and the mis-call happens anyway.

## Scope

BV-BRC holds **13 CTRV genomes, all near-complete**, against 95 other
Merhavirus records (46 near-complete *M. merida*) that have an uninterrupted L
and are called correctly today. So the feature affects 13 genomes, and all 13
are viable reference candidates — the skill's rule for splice and
transcript-edit sets is one reference per genome, not one per clade.

## The failure mode to test when building it

The generator flags it already: the ordinary `L` profile **also** hits CTRV and
stops at 8681, so `get_splice_variant_features.pl` must **replace** that call,
not add to it. If it duplicates, a truncated L has been traded for a duplicate
call, which `evaluate_module.py` treats as always a defect. Check this on CTRV
contigs specifically, before and after.

## One correction to the Rhabdoviridae build notes

Those notes list `MERHA_L_SPLICED` among "32 features with NO COLLECTION" and
recommend removing all of them. That is a category error: a `special` feature
never has a collection or a PSSM, by design. What it lacks is a nucleotide
reference set under `Splice-Variants/<Module>/<FEAT>.fasta`. Deleting it on
those grounds would discard a documented, buildable feature for the wrong
reason.
