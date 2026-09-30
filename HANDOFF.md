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

### 1. The homology rescue has never run on 7 of 14 modules

**This is the item that most affects what the modules actually call.** Shelved
30 September 2026 at Jim's instruction, and noted in all seven coverage audits
so it cannot be lost.

`rescue_unassigned.py` takes every protein that matched no annotation-string
rule, blastp's it against that module's own collections, and adopts each into
the feature it matches at ≥80% identity over ≥60% of the query. It recovers
what a string rule cannot reach — `hypothetical protein`, masses (`64kDa
protein`), positional labels (`ORF3`), typos (`Polypretein`, `Coat prtein`) —
and refuses everything below threshold rather than sweeping it into a grab-bag.

Measured by dry run, 30 September 2026:

| module | unassigned | would adopt | now binned | gain | dump |
|---|---|---|---|---|---|
| Bromoviridae | 1256 | 819 | 2901 | +28% | local |
| Ampelovirus | 1185 | 679 | 2274 | +30% | local |
| Closterovirus | 640 | 420 | 1857 | +23% | local |
| Velarivirus | 683 | 110 | 627 | +18% | local |
| Allexivirus | — | — | — | — | needs re-fetch |
| Carlavirus | — | — | — | — | needs re-fetch |
| Quinvirinae | — | — | — | — | needs re-fetch |

Why it matters beyond the counts: the uncovered population *is* the uncalled
population. On Alpharhabdovirinae, of 115 genomes with a called N and M but no
P, all 115 had a full-length recognisable P, and 85 resembled a leftover
sequence while **zero** resembled a clustered one. Neither self-recall nor a
reference panel shows this.

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

### 3. Do 1 and 2 in one pass — the order is forced

Splitting requires re-running `build_collections.py`, **which truncates the
files the rescue appends to**. A rescue done first is discarded. Correct order:

```
split the triage rules → build_collections → rescue → rebuild PSSMs
  → curate → re-evaluate → re-measure copy_num → regenerate artifacts
```

Doing them separately pays the cost twice. Measured costs:

| step | cost |
|---|---|
| `rescue_unassigned.py` | 5.5 s per module |
| PSSM rebuild, whole module | 2.2–3.6 min (from the build logs' own timings) |
| **re-evaluation** | 88.6 s per genome × 64–142 genomes = **95–210 min per module** |
| copy_num + artifacts + republish | ~15 min per module |

≈ 10 h of compute for the four local-dump modules, and it **must be serial** —
`run_gto_eval.py` refuses to start if another `annotate_by_viral_pssm.pl` is
running. Curation is manual, in Jalview, and is Jim's time, not estimable here.

### 4. Smaller open items

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
- BV-BRC writes `uniq.seq` and `uniq.id_ann` with a **trailing tab**.

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
