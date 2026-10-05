# Cardiovirus build notes

## Module boundary

One module for the genus: Ldr (zinc-finger host antagonist, 67-76 aa),
VP4-VP2-VP3-VP1, 2A, 2B, 2C, 3A, one VPg, 3C, 3D. 15 rep contigs, 100%
routing. Built 2 October 2026 at the documented defaults for every feature
(60 PSSMs); leader key renamed L -> Ldr on 3 October.

Separate from Aphthovirus because the leader slot holds a different protein:
67-76 aa zinc-finger here against a 201 aa protease there.

## Measured 4 October 2026

- Full-taxon coverage: 227 genomes >= 5 kb, 106 exemplars at 95% nt identity,
  106 routed, 93 with every feature. Lowest features VPg 97/106, 2A and 3A
  98/106, Ldr 101/106. No feature trails the best by 20 points.
- Held-out panel: 40 genomes (rep-contig genomes excluded), 40 clean, no
  genome, contig or feature flag. Cleanest panel of the ten modules.
- Self-recall 1083/1086.

## L* and 2B* -- BUILT 5 October 2026

Both are out of frame with the polyprotein, so no amount of tuning the
thirteen declared features reaches them. Both are now declared and built;
the module carries 15 features and 72 PSSMs. `special_cardio.py` extracts
both from the annotated genomes.

**L*** - alternative AUG in the leader region, different frame. Chen 1995
Nat Med 1:927 (PMID 7585219) shows mutating its start cripples demyelination;
Kong & Roos 1991 J Virol 65:3395 (PMID 2033677) mapped the initiation site.
Searching the 106 exemplars for an out-of-frame ORF >= 100 codons starting at
an AUG within 150 nt of the polyprotein start finds 41 (23 C. theileri,
6 C. ranori, 4 rodent cardiovirus, singletons elsewhere), median 156 aa.
L* is collinear with itself, so it needs only a collection and an ordinary
PSSM feature - tblastn finds an overlapping ORF in any frame. Built that way:
key LSTAR, anno "L* protein", symbol "L*", CDS, kmers 0 (it overlaps POLY),
upstream_ext 1 / downstream_ext 0. 47 sequences over 209 annotated genomes,
median 156 aa, window 125-187, bit_cutoff 70, 7 profiles (2 first-pass +
5 leftover).

**2B*** - -1 programmed ribosomal frameshift at a conserved G_GUU_UUY motif
just after the 2A/2B junction. Loughran 2011 PNAS 108:E1111 (PMID 22025686),
EMCV, 128-129 aa transframe product, small-plaque on knockout; Finch 2015
J Virol 89:8580 (PMID 26063423), TMEV, same site at 74-82% efficiency.
Reading the heptamer at its correct codon phase (the lone purine is the third
base of a codon) finds it intact in 21 exemplars, all C. rueckerti (= EMCV),
and the construct translates to 128-130 aa in every one - an independent match
to the published length. Theilovirus genomes carry the motif 16 nt before the
called 2B start. 2B* is NOT collinear, so it is declared
`special: transcript_edit` with a per-genome nucleotide reference set and no
PSSM: key 2BSTAR, anno "TransFrame protein" (reusing the Togaviridae TF
vocabulary string), symbol "2B*", mat_peptide, window 127-133.

Over all 209 annotated genomes the motif is in phase and intact in 80, every
one C. rueckerti, slip codon at position 8 in 77 of them. Reference set:
80 constructs, all ending in a stop (the product terminates there, so the
stop codon is kept and the protein carries a trailing *), 0 internal stops,
0 ambiguous bases, 129 aa x76 / 130 x3 / 131 x1. The N-terminus is
PFMAKPKKQVF..., which is 2B's own start - the shared 11-12 residues the
frameshift papers describe, reproduced without being imposed.

SELF ROUND-TRIP: 12 donor genomes through the full chain (make_gto ->
annotate_by_viral_pssm-GTO -> get_transcript_edited_features) recover their
own protein BYTE-IDENTICALLY, 12/12, 129 aa, none missing, none differing.
Note the called feature reaches the GTO but not the flat feature table -
the table has no representation for a special feature, so evaluate_coverage
cannot see 2B* and reports it as absent rather than as low.

Five distinct source strings name L* and one names 2B*, so BV-BRC records
both proteins the module does not call.

All six PMIDs above are model-proposed and belong in PMID_claude_generated
when the features are built.
