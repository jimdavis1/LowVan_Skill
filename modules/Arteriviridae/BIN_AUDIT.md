# Bin audit -- Arteriviridae

Does every collection hold ONE protein, and its OWN protein?

**Unresolved: 0**  (defects 11, signed off 11, paralogues 1, collections 24)

## FRAGMENT -- 11 bin(s), 1999 sequence(s)

| bin | n | median aa | finding | signed off |
|---|---|---|---|---|
| `NSP1` | 6 | 260 | partial of ORF1A (99.4% id over 15% of ORF1A) -- identity only, no co-deposit to confirm it | mature peptide of pp1a, not a fragment: nsp1 is released by its own PCPalpha/beta autoprotease (PMID 19553313). Partial-of-ORF1A at high identity is the DEFINING property of every nsp, so the FRAGMENT test cannot discriminate here. Discriminated instead by (a) junction chemistry -- every member begins M and ends G at the measured autoproteolytic site -- and (b) the literature length window in length_refs.json, both applied in build_collections.py before this audit ran. |
| `NSP1A` | 137 | 180 | partial of ORF1A (100.0% id over 8% of ORF1A) -- identity only, no co-deposit to confirm it | mature peptide of pp1a (PRRSV/LDV/SHFV nsp1alpha, PCPalpha autocleavage at M...MH, PMID 19553313). Same argument as NSP1; junction M/MH and literature window enforced at bin time. |
| `NSP1B` | 140 | 205 | partial of ORF1A (100.0% id over 9% of ORF1A) -- identity only, no co-deposit to confirm it | mature peptide of pp1a (nsp1beta, PCPbeta autocleavage, PMID 19553313). Junction AS/G and literature window enforced at bin time. |
| `NSP2` | 964 | 1059 | partial of ORF1A (100.0% id over 45% of ORF1A) -- identity only, no co-deposit to confirm it | mature peptide of pp1a, the largest one. Its C-terminus is released by PLP2, nsp2's own papain-like cysteine protease, which cleaves G|G (PMID 11602729). Partial-of-ORF1A at 100% identity is the DEFINING property of every nsp, so the FRAGMENT test cannot discriminate here; discriminated instead by junction chemistry and the literature window 700-1400 aa, applied at bin time. n=964: 409 deposits retained, plus 555 members re-derived from the genomes on the nsp1beta/nsp3 adjacency and the PLP2 G|G site (nsp2fix/extract_nsp2_nsp3.py, 21 Nov 2026). 229 deposits ending in Isoleucine were dropped to nsp2fix/NSP2.dropped_Iterm.faa: PLP2 cleaves after Gly, so an Ile terminus is a submitter ORF convention, not a cleavage product, and those members were clipping the profile 216 codons short at an internal LPG|GARY site. All 555 extracted members end in Gly. |
| `NSP3` | 543 | 230 | partial of ORF1A (100.0% id over 10% of ORF1A) -- identity only, no co-deposit to confirm it | mature peptide of pp1a. Its N-terminus is the PLP2 G|G site and its C-terminus the nsp4 3CLSP cut after Glu (PMID 11602729). n=543: 8 deposits plus 535 re-derived on the same adjacency as NSP2. All 535 extracted members start Gly and end Glu, matching all 8 deposits exactly -- independent confirmation of the NSP2 repair from the opposite end of the junction. Literature window 200-350 aa enforced at bin time. |
| `NSP4` | 7 | 204 | partial of ORF1A (100.0% id over 8% of ORF1A) -- identity only, no co-deposit to confirm it | mature peptide of pp1a, the 3C-like serine protease itself, G/E termini (PMID 11602729). Junction G/E and literature window enforced at bin time. |
| `NSP5` | 10 | 170 | partial of ORF1A (100.0% id over 7% of ORF1A) -- identity only, no co-deposit to confirm it | mature peptide of pp1a, 3CLSP cleavage GS/E (PMID 11602729). Junction and literature window enforced at bin time. |
| `NSP6` | 5 | 16 | partial of ORF1A (100.0% id over 1% of ORF1A) -- identity only, no co-deposit to confirm it | mature peptide of pp1a. Median 16 aa is CORRECT, not a truncation: nsp6 is a 16-residue peptide (PMID 11602729). The 1%-of-ORF1A coverage this audit reports is the expected consequence of its true size. Junction G/E and literature window 16-16 enforced at bin time. |
| `NSP7A` | 14 | 123 | partial of ORF1A (100.0% id over 7% of ORF1A) -- identity only, no co-deposit to confirm it | mature peptide of pp1a, nsp7alpha, S/E termini (PMID 11602729). Junction and literature window enforced at bin time. |
| `NSP7B` | 164 | 120 | partial of ORF1A (100.0% id over 5% of ORF1A) -- identity only, no co-deposit to confirm it | mature peptide of pp1a, nsp7beta, AN/E termini (PMID 11602729). Junction and literature window enforced at bin time. |
| `NSP8` | 9 | 45 | partial of ORF1A (100.0% id over 2% of ORF1A) -- identity only, no co-deposit to confirm it | mature peptide of pp1a, the C-terminal product of pp1a -- no downstream cleavage, so its C-terminus is the pp1a stop (PMID 30353315). Median 45 aa is the literature size. Junction AG/- and literature window enforced at bin time. |

## SPLIT -- 1 bin(s), 39 sequence(s)  (reported, not retired)

| bin | n | median aa | finding | signed off |
|---|---|---|---|---|
| `NSP12` | 39 | 143 | 9 groups, largest holds 44% of the bin -- this is more than one protein and resembles no named feature |  |

