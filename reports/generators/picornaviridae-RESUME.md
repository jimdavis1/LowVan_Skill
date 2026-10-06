# Picornaviridae — 14 modules complete and QC'd (6 October 2026)

Repo: https://github.com/jimdavis1/LowVan_Skill — head `f3ca5fed`.
Runtime of record: `Alsuviricetes/Viral_Annotation` (70 modules).
Work trees `work/<Module>/`, pipeline `work/_pipeline/`, QC `work/_qc/`.

| module | feat | PSSMs | refs | routed | panel clean |
|---|---|---|---|---|---|
| Enterovirus | 12 | 701 | 25 | 2856/2856 | 36/40 |
| Aphthovirus | 13 | 132 | 13 | 431/431 | 39/40 |
| Kobuvirus | 12 | 115 | 23 | 313/313 | 32/40 |
| PICO_LDR_VP4 | 15 | 100 | 25 | 182/201 | 36/39 |
| Hepatovirus | 12 | 86 | 25 | 113/116 | 36/39 |
| PICO_NOLDR_VP4 | 13 | 72 | 25 | 83/83 | 26/40 |
| PICO_LDR_VP0 | 13 | 70 | 25 | 119/136 | 29/39 |
| Cardiovirus | 15 | 67 | 15 | 106/106 | 40/40 |
| Parechovirus | 11 | 66 | 25 | 244/245 | 39/39 |
| Sapelovirus | 13 | 63 | 25 | 208/209 | 34/40 |
| Teschovirus | 13 | 49 | 7 | 160/160 | 39/40 |
| PICO_NOLDR_VP0 | 11 | 28 | 14 | 53/53 | 38/40 |
| Senecavirus | 13 | 22 | 3 | 10/10 | 39/40 |
| Avihepatovirus | 13 | 22 | 9 | 22/22 | 39/40 |

**Family coverage: 29 of 66 genera, 94.4% of 39,757 genomes.**
44 artifacts published; `ARTIFACTS.md` has the links.

## QC sweep, 5 Oct — all 14 modules (`work/_qc/`)

install --check 14/14 clean · json_canon 14/14 canonical ·
rebuild_pssms **1,591 profiles, 0 stale** · check_dlits clean ·
qc_cross_feature **0 real mislabels**.

`qc_cross_feature`'s 10 rows are all a mature peptide matching POLY, which is
definitional. **The tool needs the parent CDS excluded before the next
polyprotein family** or it reports one false positive per module.

## New features built

- **Cardiovirus L\*** — out-of-frame AUG in the leader region. Ordinary PSSM
  feature; tblastn reads all six frames. Called 43/106, and the distribution
  is the proof: 23/53 *C. theileri*, 8/8 *ranori*, **1/26 EMCV** — EMCV does
  not encode it.
- **Cardiovirus 2B\*** — −1 frameshift transframe product. `transcript_edit`,
  80 per-genome nucleotide references, no PSSM. 129 aa ×76 against a published
  128–129. **Self round-trip 12/12 byte-identical.**
- **Avihepatovirus 2A3** — RefSeq lumps 2A2+2A3 as one 285 aa record; split at
  the HLPR motif the DHAV-3 reference annotates separately.
- **Aphthovirus VPg** — every species now called (reedi 32/32, bogeli 23/23,
  burrowsi 5/5, FMDV 360/368, **422/431** overall, zero duplicates).

## The rule the grouped modules produced

**A genus contributing 0–1 sequences to a feature gets 0 calls for it.**
Exact, no exceptions, and genome count does not predict it:

- Shanbavirus: 44 exemplars, **1** sequence per mature peptide except 3D (10).
  Called on POLY and 3D only.
- Hunnivirus: **122** polyproteins, **0** annotated VP1. VP1 called 0/76.
- PICO_NOLDR_VP0: four genera at 1–3 sequences; a leftover-pass rebuild moved
  the result by one genome.

**Count sequences per genus per feature BEFORE building.** Nothing downstream
reports it — the profiles that do get built recover their own training data
perfectly.

## Outstanding

1. **Projection rebuild of the four grouped modules.** The technique that took
   the Sapelovirus leader from 6/209 to 170/208 was applied to five genus
   modules and not to these. Hunnivirus (122 polyproteins, 0 annotated VP1) is
   the clearest case in the family.
2. **PICO_NOLDR_VP4 VP4 should be split.** Per-genus medians are Cosavirus 68,
   Dicipivirus 44, Tremovirus 20. Flagged as borderline before building and
   left alone; Tremovirus gets 0/10. The 2A split at the same spread (31 vs
   141–164) worked.
3. **353 of 526 no-genus genome types do not route.** 173 do, and get 10–12
   features each. (An earlier figure of "94% unrouted" was wrong — it came
   from a 10 kb completeness floor against a ~7.5 kb genome.)
4. **Theilovirus short-ORF 2B\*** undeclared; **Coypu aphthovirus** VPg has no
   source record.
5. **Every citation is model-proposed** and sits in `PMID_claude_generated`.
   None confirmed by a curator.

## Rules learned the hard way

- **Never lower `-mi`.** It redefines "same protein". Lower the MEMBER floor
  and let the leftover pass (which never lowers identity) reach divergent
  members. Two-member clusters are fine; a floor of 1 is not — the one
  exception is Senecavirus 2A, approved and recorded in its BUILD_PARAMS.
- **A declared feature with no profile is silent at runtime.**
  `install_module.py` now refuses it (`--force-known-bad` to override). This
  reached a shipped module twice in one day before the gate existed.
- **Key ≠ gene symbol.** The feature table reports the SYMBOL. Counting by key
  returned 0 for Cardiovirus L\* and Enterovirus VPg and put a wrong number on
  a published page.
- **One annotation string = one symbol.** `VPg` was `VPG` in 8 modules and
  `VPg` in 6. Fixed; 21 Picornaviridae strings now registered in the
  vocabulary.
- **Never edit a script while it is running** — bash reads incrementally.
- **Never `taskpolicy -b` the pipelines**: measured 5 Oct, coverage went from
  4 genomes/min to 0 and Box got no faster. The Box fix that mattered was not
  rsyncing per-genome `ann/` trees to Box.

## Environment

Scratch `~/lowvan_scratch` (local; `/private/tmp` is wiped on reboot and cost
a day on 3 Oct). Long jobs in `screen -r lowvan`. `nice_daemon.sh` keeps
compute at nice +5.
