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

## The coat protein has two real forms, and one of them is Mandarivirus

`coat protein` in Potexvirus is not one length distribution. 2,029 records sit
at 200-239 aa and a separate mode of **118 sits at 320-359**, with a third
group of 290 Allexivirus records at 120-159. All three say only "coat protein"
or "capsid protein", so the string decides nothing and length alone would
throw away real sequences. Homology settles all three.

**Allexivirus 120-159 aa: partial CP records.** 40 of 40 sampled hit the
confident Allexivirus CP core (240-279 aa) at a median **99% identity**, and
none hits NABP or the 40K at all. The alignment geometry says fragment: the
query aligns end to end (q1-135, q1-147) against the *back half* of the
subject (s125-259, s113-259), covering 52-60% of it. These are truncated
deposits of the real CP, so they stay out of the training set. A profile built
from them would be a truncated copy of the CP profile at 99% identity, which
is the one thing the curator's checklist forbids outright.

**Potexvirus 320-359 aa: Mandarivirus.** 40 of 40 hit the CP core at 86-87%,
but the geometry is the opposite of a fragment and not a mis-called start
either: every one aligns **q91-323 against s1-233**, so the first 90 residues
are an N-terminal extension the core form does not have. It is not an
over-called upstream Met -- all 40 begin with Met and **none has a Met at
residue 91**, so there is no shorter reading to prefer. The extension is real.

The genomes are **Indian citrus ringspot virus** and **Citrus yellow vein
clearing virus**: Mandarivirus, whose CP is ~34 kDa against a potexvirus
22-27 kDa. `PARTITIONING.md` puts Mandarivirus in this module, and an earlier
query for it found nothing because **BV-BRC labels these genomes Potexvirus**.
The module is recovering its own Mandarivirus members in spite of the genus
column -- the same hazard that put a genuine NABP inside a "Trichovirus" in
Betaflexiviridae_MP.

So the CP window is wide (180-400) and covers both real forms. They differ by
90 residues at 87% identity, so mmseqs separates them at `-mi 0.8` and each
gets a profile with its own clean N-terminus; `qc_truncation_symmetric.py` will
report the containment, and at 87% it is a divergent form to keep rather than a
truncation to retire -- the same call made for CP cluster 21 in
Betaflexiviridae_MP.

## Windows, after all of the above

| key | window | note |
|---|---|---|
| REP | 1,300-1,800 | excludes the 160-aa "polymerase" records |
| TGB1 | 200-280 | |
| TGB2 | 85-135 | overlaps NABP; homology separates them |
| TGB3 | 50-120 | Potexvirus; excludes the Allexivirus 40K |
| CP | 180-400 | both full-length forms; excludes the 120-159 partials |
| ALLEXI_40K | 320-400 | 289 records: 210 "tgb3", 29 "40 kda protein", 9 "serine-rich protein", 8 "p42" |
| NABP | 90-145 | Allexivirus/Lolavirus |
