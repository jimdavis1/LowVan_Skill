# Alphaflexiviridae triage

5,420 genomes, 14,619 features, 6,774 unique sequences, 388 distinct
(string, genus) pairs. One module, per `work/PARTITIONING.md`: Potexvirus,
Allexivirus and Lolavirus share the replicase + triple gene block + CP layout,
with Allexivirus and Lolavirus adding two accessories.

Rep contigs: **84.9% of 1,348 genomes at the 25 budget** (clusters at 70%
identity: 102 total, largest 407). Ships. Splitting into Potexvirus (93.4%)
and Allexivirus+Lolavirus (97.4%) would reach ~95% but spend 50 references
instead of 25; the single-module decision stands unless the curator prefers
to spend them.

## Six features

| key | protein | length | genera |
|---|---|---|---|
| REP | replicase, RdRp domain | ~1,417 Potex / ~1,543 Allexi | all |
| TGB1 | triple gene block protein 1 | 224-265 | all |
| TGB2 | triple gene block protein 2 | 90-125 | all |
| TGB3 | triple gene block protein 3 | 63-110 | Potexvirus |
| CP | coat protein | 215-260 | all |
| ALLEXI_40K | 40K protein | 357-364 | Allexivirus |
| NABP | nucleic acid-binding protein | 100-130 | Allexivirus |

## Three strings that must not be trusted on their own

Every one of these is the U-number trap: a positional label covering
non-homologous proteins. All were caught by histogramming length against
(string, genus), which is why that step is not optional.

**`tgb3` in Allexivirus is bimodal and mostly not a TGB3.** n=233, with
p5=99 but p25 through p95 at 357-364 aa. A potexvirus TGB3 is 63-89 aa.
Allexiviruses carry a 40K protein at the position where potexviruses carry
TGB3, and submitters have labelled it by position. Binning on the string
would have put 233 sequences averaging 364 aa into a collection whose real
members are under 90, and the resulting profile would have called
"triple gene block protein 3" on the 40K locus across the whole genus.

**`rna dependent rna polymerase` without the hyphen is bimodal.** n=84,
p25-p50 = 160 aa, p75-p95 = 1,417-1,456 aa. The hyphenated form is clean at
1,417. Roughly half these records are not a polymerase at all.

**`polymerase` is mostly not the polymerase.** n=43, p50 = 230 aa, which is
TGB1-sized, with real 1,439 aa polymerases only in the tail.

`hypothetical protein` (n=45) spreads from 52 to 325 aa with no mode. It
stays unassigned; there is nothing to bin it on.

## Two pairs length cannot separate

Length gates did the work for the tobamovirus replicases and they will not
work here:

- **TGB2 (~103) against NABP (~104-127).** Both Allexivirus, both small.
- **CP (215-260) against TGB1 (224-265).** Overlapping across the range.

Both pairs are genuinely different proteins, so profiles built from correctly
labelled sequences separate them by homology. The discriminator is
`bit_cutoff` plus `qc_cross_feature.py`, not `min_len`/`max_len`, and the
length windows must be set wide enough not to lie about what the feature is.

## Rule shape

Bin on (string, genus, length window) together. A string-only rule is wrong
for at least three of the six features, and a length-only rule is wrong for
two of them.
