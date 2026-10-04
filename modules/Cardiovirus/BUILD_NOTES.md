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

## Known gaps: two proven proteins not declared

Both are out of frame with the polyprotein, so no amount of tuning the
thirteen declared features reaches them. Sequences extracted and counted;
neither feature is built.

**L*** - alternative AUG in the leader region, different frame. Chen 1995
Nat Med 1:927 (PMID 7585219) shows mutating its start cripples demyelination;
Kong & Roos 1991 J Virol 65:3395 (PMID 2033677) mapped the initiation site.
Searching the 106 exemplars for an out-of-frame ORF >= 100 codons starting at
an AUG within 150 nt of the polyprotein start finds 41 (23 C. theileri,
6 C. ranori, 4 rodent cardiovirus, singletons elsewhere), median 156 aa.
L* is collinear with itself, so it needs only a collection and an ordinary
PSSM feature - tblastn finds an overlapping ORF in any frame.

**2B*** - -1 programmed ribosomal frameshift at a conserved G_GUU_UUY motif
just after the 2A/2B junction. Loughran 2011 PNAS 108:E1111 (PMID 22025686),
EMCV, 128-129 aa transframe product, small-plaque on knockout; Finch 2015
J Virol 89:8580 (PMID 26063423), TMEV, same site at 74-82% efficiency.
Reading the heptamer at its correct codon phase (the lone purine is the third
base of a codon) finds it intact in 21 exemplars, all C. rueckerti (= EMCV),
and the construct translates to 128-130 aa in every one - an independent match
to the published length. Theilovirus genomes carry the motif 16 nt before the
called 2B start. 2B* is NOT collinear, so it needs `special: transcript_edit`
and a per-genome nucleotide reference set, not a PSSM.

Five distinct source strings name L* and one names 2B*, so BV-BRC records
both proteins the module does not call.

All six PMIDs above are model-proposed and belong in PMID_claude_generated
when the features are built.
