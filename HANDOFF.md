# Handoff — LowVan Alsuviricetes, as of 30 September 2026

Point a fresh session at this file. It is the current state of play, not a
history; `work/BUILD_ORDER.md` holds the chronological record and the reasoning
behind each decision.

## Where things are

| thing | path |
|---|---|
| the kit (the repo that matters) | `LowVan_Skill/`, remote `https://github.com/jimdavis1/LowVan_Skill` |
| per-module JSON, PSSMs, rep contigs | `LowVan_Skill/modules/<Module>/` |
| the installed runtime this machine annotates with | `Viral_Annotation/Viral_PSSM.json` + `Viral_Annotation/Viral-PSSMs/` |
| build scratch, dumps, collections, measurements | `work/<Module>/` |
| published pages and their generators | `LowVan_Skill/reports/`, index in `LowVan_Skill/ARTIFACTS.md` |
| the skill itself | `LowVan_Skill/skill/` |

**Push only to `jimdavis1/LowVan_Skill`.** CEPI-dxkb is downstream, is behind
the kit, and syncing it is Jim's call.

## Status

**All 14 Alsuviricetes modules are built, evaluated, installed, and shipped.**
Each has its four pages published (string collapse, PSSM registry, vocabulary
saturation, coverage audit) — 56 of 56, none missing.

Allexivirus · Alphaflexiviridae · Ampelovirus · Bromoviridae · Capillovirus ·
Carlavirus · Closterovirus · Crinivirus · Marafivirus · Quinvirinae ·
Tobamovirus · Trivirinae · Tymovirus · Velarivirus

**The homology rescue is done.** Backfilled 1 October 2026 on the six modules
that never had it; Allexivirus turned out never to have needed it. Measured
against a control rather than against the old baseline, it is worth **+20
genomes across 770**, 70.3% -> 72.9% genuinely clean. See "What the rescue was
worth" below.

Invariants that held as of the last commit, and are worth re-checking after any
change:

```bash
# every module JSON canonical, and identical to the installed runtime block
python3 LowVan_Skill/skill/scripts/json_canon.py --check \
        LowVan_Skill/modules/*/*_Viral_PSSM.json Viral_Annotation/Viral_PSSM.json

# zero VARIANTs, zero symbol mismatches, zero style flags across all 14
for m in Allexivirus Alphaflexiviridae Ampelovirus Bromoviridae Capillovirus \
         Carlavirus Closterovirus Crinivirus Marafivirus Quinvirinae \
         Tobamovirus Trivirinae Tymovirus Velarivirus; do
  python3 LowVan_Skill/skill/scripts/check_annotations.py \
    --json LowVan_Skill/modules/$m/${m}_Viral_PSSM.json \
    --vocab LowVan_Skill/skill/assets/annotation-vocabulary.tsv
done
```

Closed by curator decision and **not to be reopened or generalised**:
Carlavirus at 76.3% routing, Bromoviridae at 74.2% with a 25.8% no-call gap,
Quinvirinae `P14`-in-`CP` and `ORF2A`-in-`TGB1` as real out-of-frame
overlapping ORFs. Each is recorded in its module's `BUILD_NOTES.md`.

## Deferred, in priority order

### 1. ~~The homology rescue~~ — done, 1 October 2026

Ran on Ampelovirus, Bromoviridae, Carlavirus, Closterovirus, Quinvirinae and
Velarivirus. **Allexivirus never needed it**: its profiles are selected from
Alphaflexiviridae, whose rescue did run, and re-running it over the parent
adopts 3 more sequences, all *Potexvirus*, none *Allexivirus*. Six modules, not
seven.

What it adopted, and what that was worth:

| module | adopted | training seqs | profiles | control -> after | +/- |
|---|---|---|---|---|---|
| Carlavirus | 425 of 524 | 3139 -> 3463 | 316 -> 338 | 81.5% -> **85.7%** | +12 / -0 |
| Quinvirinae | 664 of 836 | 4026 -> 4374 | 116 -> 115 | 79.6% -> **83.3%** | +4 / -0 |
| Closterovirus | 420 of 640 | 1857 -> 2091 | 80 -> 97 | 27.0% -> **36.5%** | +13 / -6 |
| Bromoviridae | 816 of 1253 | 2904 -> 3568 | 124 -> 135 | 84.9% -> 84.9% | +0 / -0 |
| Ampelovirus | 679 of 1185 | 2274 -> 2819 | 139 -> 147 | 52.5% -> 50.5% | +3 / -5 |
| Velarivirus | 110 of 683 | 627 -> 705 | 67 -> 69 | 50.0% -> 48.4% | +0 / -1 |
| **total** | | | 842 -> 901 | **70.3% -> 72.9%** | **+32 / -12** |

**Read the control column, not the old baseline.** `measurements/quality.tsv`
predates the 30 September `copy_num` and vocabulary edits, and on three modules
those edits dominate the rescue completely — dropping `copy_num` from six
Closterovirus accessories retires 184 "missing essential feature" flags on its
own and makes the rescue look like 0% -> 36.5% when it is 27.0% -> 36.5%. The
control is the shipped module scored on the same panels with the current JSON.

**11 of the 12 regressions are one pre-existing defect**: a profile trained on
the ORF1a+ORF1b frameshift fusion matches both halves and calls
`Replicase polyprotein` twice. POLY was already 27% fusion-length before the
rescue. Net across the two closterovirid modules the rescue *helps* —
duplicate POLY calls go 51 -> 40 — but it cannot fix the cause, which is the
`transcript_edit` item under "Smaller open items". The twelfth regression is
Velarivirus HSP90, below.

Panels are rebuilt from each module's own `measurements/quality.tsv` and live
in `work/<Module>/gto/`; `eval_control.tsv` and `eval_after.tsv` sit beside
them.

### 2. Seven features are mass bins holding more than one homology group

**Jim's rule: a bin must be homologous; it cannot hold different features.**
These are therefore defects, not open questions. No pair of groups within a bin
shares a blastp HSP at 1e-3.

| bin | groups |
|---|---|
| Closterovirus P18 | CTV 75 \| Thesium 2 |
| Closterovirus P19 | GLRaV-2 51 \| wheat + Triticum 2 |
| Closterovirus P20 | CTV 54 \| BYV 9 |
| Closterovirus P23 | CTV 94 \| Thesium 2 \| Raspberry 2 |
| Ampelovirus P7 | GLRaV-3 47 \| GLRaV-1 29 \| 3 × GLRaV-13 group of 3 |
| Ampelovirus P21 | GLRaV-3 + PMWaV-6 40 \| GLRaV-1 38 |
| Ampelovirus P24 | GLRaV-1 42 \| SCMMV 11 \| PMWaV 9 |

`P20` is the clearest: CTV p20 is a silencing suppressor, BYV p20 is an
Hsp70h-interacting vascular-transport protein, and the literature puts CTV
p20's real homolog at BYV *p21*.

**Undecided:** several groups are 2–3 sequences (Thesium, Raspberry, wheat, two
GLRaV-13 groups). Proposed default, not yet approved: send groups below 5
sequences to `NOT_MODELLED.tsv` rather than give each its own feature, making it
7 bins → 14 features instead of 19. Ask before implementing.

### 3. The order, corrected

Item 1 is done, so the splits in item 2 now stand alone — but doing them still
re-runs `build_collections.py`, which **truncates the files the rescue appended
to**. Re-run the rescue immediately afterwards with the same windows, then
diff every collection count against the previous build.

The order in the previous handoff was also **missing a step**, and it cost
7 genomes before it was caught:

```
split the triage rules -> build_collections -> rescue -> rebuild PSSMs
  -> REFRESH min_len/max_len FROM THE COLLECTIONS   <-- this one
  -> curate -> re-evaluate -> re-measure copy_num -> regenerate artifacts
```

`min_len`/`max_len` in the module JSON are **derived from the collection** —
verified across 36 features, where every untouched bound matches its
collection's min and max exactly. A rescue that adds sequences therefore
invalidates them, and the stale bound makes `viral_genome_quality.pl` flag
"Feature is too short" on a protein the module now models correctly.
13 features across all six modules needed refreshing; `work/fix_bounds.py`
does it and touches nothing else.

Measured costs, revised — evaluation is **far cheaper than the previous
estimate of 88.6 s per genome**:

| step | cost |
|---|---|
| `rescue_unassigned.py` | 5.5 s per module |
| PSSM rebuild, whole module | 1-3 min |
| re-evaluation | **~16 s per genome**, 17-96 min per module |
| the whole six-module pass, rescue to shipped | ~9 h wall clock, most of it evaluation |

Evaluation **must still be serial** — `run_gto_eval.py` refuses to start if
another `annotate_by_viral_pssm.pl` is running.

### 4. Smaller open items

- **Velarivirus's rescue is a net loss of one genome — keep or revert?** Default
  is keep: +12% training data against one genome whose `HSP90` is called at
  497-505 aa where the collection minimum is 506 and the source says 513
  (`p60 protein`, genome 1654603). The cause is HSP90's clustering changing
  (8 alignments to 7), not shorter training data.
- **Velarivirus's real gap is declared features, not training data.** Of the 573
  sequences the rescue refused, 185 are `hypothetical protein` and 85 are named
  by mass — `p25`, `p21`, `p60`, `p27`, `p4`, `P22` across *Velarivirus* and
  *Olivavirus*. The module declares seven features and no accessory, so they
  have nowhere to go. Accessory-discovery clustering, as was done for
  Closterovirus P6-P33, is the route; it has not been scoped.
- **`-m` in `features.json` never reached the command line.** The file writes
  `"-m": "3"`; `DEFAULTS` is keyed `m`, so `params.update()` adds a key the
  command never reads. The run keeps `-m 5` while `BUILD_PARAMS` records
  `departures: {"-m": "3"}` and the build log prints `[-m=3]`. **28 features in
  seven modules** — Ampelovirus, Carlavirus, Closterovirus, Pegivirus,
  Pestiviridae, Quinvirinae, Velarivirus — claim a member floor they were not
  built with. Fix written at `work/patch_build_features.py`, **deliberately not
  applied**: it changes clustering on those features, and applying it during
  the rescue would have made neither effect measurable. Decide separately.
- **The `work/<Module>/*_Viral_PSSM.json` files are stale.** They predate the
  30 September vocabulary sweep — `Coat protein`/`CP` for
  `Nucleocapsid protein`/`N`, no Closterovirus CP/CPM paralogy comment, a
  missing Velarivirus `copy_num`. Installing from `work/` silently reverts the
  sweep. The kit is the source of truth for the JSON as well as the scripts;
  copies kept as `*.json.workdir-stale`.

- **Quinvirinae P14's call rate does not reconcile**: §05 of its audit says
  60.7% of 150 routed genomes (91 calls); the only run on disk says 22.2% of 108
  (24 calls), cross-checked against that run's feature-count distribution. They
  disagree on ORF5a too (8.7% vs 25.9%), so it is not rounding, and the
  150-genome output is not on disk. Both figures are named on the page.
  Reproduce one before quoting either.
- **`POLY` collides with Orthoflavivirus's `POLY`** (`Genome polyprotein`) across
  6 modules, and Ampelovirus `P7` with Fimoviridae's `P7`. Reported by
  `check_annotations.py` as collisions; no decision taken.
- **Endornaviridae** tabled at ~30% coverage at budget 25, and the 18 small taxa
  Jim scoped out.
- **The four rhabdovirus modules** (Alpharhabdovirinae, Betarhabdovirinae,
  Dichorhavirus, Novirhabdovirus) are in the kit but **not installed in this
  runtime**, and their PSSM directories are not here either. Almost certainly
  just this being the Alsuviricetes working runtime, but nothing here exercises
  them.

## Environment traps, all hit at least once on this project

- **`mmseqs` lives only in the `lowvan` conda env.** Outside it,
  `fasta-cluster-pssm-2.pl` "succeeds" in 0.06 s, produces zero PSSMs and leaves
  a hash-named directory. Run builds under `conda run -n lowvan`.
  `FCP_Main_Utils.pm` needs `PERL5LIB=/Users/jdavis/Viral_Annotation/Other_Scripts`.
- **`git` invoked from inside Python** resolves to a broken xcrun shim and
  returns empty output with rc=1, so verification scripts silently pass. Shell
  out to `/usr/bin/git`, or dump `git show` output to a file first from bash.
- **Never run two evaluations at once.** This machine hit load 138–148 three
  times. `pkill -f run_eval.sh` leaves reparented children alive.
- **Box sync flattens timestamps.** Do not estimate durations from file mtimes —
  a build.log "span" read 3617 min for 219 s of work. Use the per-feature
  timings the build logs record, and work in a temp dir, syncing Box at the end.
- **Re-running `build_collections.py` discards the rescue.** Documented, and
  tripped anyway. Re-run the rescue immediately after, and diff every
  collection's count against the previous build.
- BV-BRC writes `uniq.seq` and `uniq.id_ann` with a **trailing tab**. That is
  what crashed `rescue_unassigned.py` on every Alsuviricetes dump; the kit's
  copy is fixed, but `work/rescue_unassigned.py` was a stale pre-fix duplicate
  and silently reintroduced the crash. It has been synced.
- **A bulk-cleaned dump is recoverable from its `dump.log`.** The log holds the
  whole skeleton verbatim in three blocks — genome ids, representative feature
  ids, md5s — verified line-for-line against Ampelovirus, whose real files
  still exist. Only metadata, products and sequences have to be re-fetched.
  Carlavirus and Quinvirinae came back with every count matching the original
  log and zero md5 drift against the 24 September snapshot.
- **Trust the collection over the build log.** Bromoviridae's `REP2.fasta` held
  587 sequences against a build log saying 590, and the three missing ones were
  restored on that basis. They were removed deliberately: blastp puts all three
  at 97-99.6% identity to `REP1` across all 993 residues, so they are `1a`
  proteins mislabelled `RNA-dependent RNA polymerase` at source. They raised
  `REP2`'s `max_len` to 993 and cost 8 genomes. **A collection whose maximum
  length jumps after a rescue is a binning change, not a coverage gain** —
  check the extremes, not just the count.
- **Cluster numbers are not stable across rebuilds.** `POLY.11` before a
  rebuild is not `POLY.11` after it; Ampelovirus went from 24 to 29 POLY
  profiles with 6 and 9 retired. Compare profiles by membership, never by name.

## A working copy lives outside the project

`/private/tmp/lowvan_scratch/repo` is a full checkout of the upstream
`Viral_Annotation` repo with the pipeline patches applied and every module's
built PSSMs and rep contigs in place (~1.8 GB), and `/private/tmp/lowvan_runtime`
is a second copy. **Nothing unique should live there.** As of 30 September 2026
three pipeline edits did, and were captured into `LowVan_Skill/patches/` with
`patches/README.md` describing each — `viral_genome_quality.pl`,
`annotate_by_viral_pssm-GTO.pl` and `get_splice_variant_features.pl`. Before
clearing either directory, run `git status` inside it and diff anything modified
against the patches here.

## How Jim wants to work

- Don't interrupt for anything that isn't a real decision. Publish modules when
  ready.
- Ask **one** thing at a time, with a default, and show the measurement behind it.
- Answer his questions as answers. Facts go to the literature, not to him.
- A decision on one taxon is **not** a rule. Do not generalise it into the skill.
- `copy_num` marks a feature essential — measure the call rate after the module
  runs, never assume it.
