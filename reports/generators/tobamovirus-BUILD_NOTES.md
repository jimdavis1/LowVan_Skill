
## 2026-09-17 — annotation string collision, and why the first quality run read 0% good

The first GTO quality run scored **0 of 89 genomes good**, with
`Genome has too many HSPs for: Nonstructural polyprotein; Count = 2` on 87
genomes and a scatter of 52 "Feature is too short" / 38 "Feature is too long".

Nothing was wrong with the profiles, the cutoffs, the rep contigs or the
genomes. `REP126` and `REP183` had both been declared with the annotation
string `Nonstructural polyprotein`, and `viral_genome_quality.pl` keys its
`%essential` hash by the annotation string rather than by the feature key:

```perl
$essential->{$anno}->{max_len} = ...;   # second feature read wins
```

So one feature's length window silently replaced the other's, and `%anno_count`
summed both features' calls against the single surviving `copy_num`. The script
runs once per genome, so Perl's hash randomisation picked a different winner
each time — 52 genomes where REP183's 1369-1805 window applied and the 126K
call (median 1115 aa) read too short, 38 where REP126's 955-1259 window applied
and the 183K call (median 1616 aa) read too long. 52 + 38 = 90, one flagged
feature per genome: the coin flip, not biology.

**The fix is two distinct strings**, and fixing it exposed a second error. The
tobamovirus replicase is **not** proteolytically processed — there is no
protease in the tobamovirus genome, and the 126K and 183K proteins are separate
translation products of the same ORF, not cleavage products of a precursor. So
`polyprotein` was factually wrong for both features regardless of the
collision. The source data says what they are:

| feature | top BV-BRC strings | now |
|---|---|---|
| REP126 | small replicase subunit (76), 129K protein (68), replication-associated protein (58) | `Methyltransferase and helicase replication protein`, symbol `126K` |
| REP183 | RNA-dependent RNA polymerase (167), 186K protein (72), 183 kDa protein (39) | `RNA-dependent RNA polymerase`, symbol `183K` |

REP183 reuses a string already in the controlled vocabulary (Togaviridae nsP4);
REP126 adds one new string, in the compound-functional style the coronavirus
entries already use. Four `Tobamovirus` rows were added to
`annotation-vocabulary.tsv` with the PubMed column left **blank**, because every
citation in this module is model-proposed and lives in `PMID_claude_generated`.

**Tooling.** `check_annotations.py` now fails with a non-zero exit on two
features in one module sharing an `anno` string, and the rule is written up in
`references/json-schema.md` under "The annotation string is a key, not a label".
A sweep of all 38 installed modules found no other instance.
