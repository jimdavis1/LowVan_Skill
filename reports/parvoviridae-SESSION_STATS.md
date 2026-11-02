# Session statistics — Parvoviridae, 1–2 November 2026

Session id `6863ff3d-4d0e-4b62-92a1-55d96589b743`. Written 2 November 2026 at
the close of the Parvoviridae work, before Adenoviridae was started.

**Read the caveat in "Tokens" before quoting any number from it.** The
measured figures are taken from the session transcript and from `git`; the
derived ones are marked as derived.

## Tokens

| | |
|---|---|
| context window | **15,000,000 tokens** |
| remaining at close | ~14,997,000 |
| **output tokens** (measured) | **1,294,006** |
| **cache-creation input tokens** (measured) | **3,621,045** |
| cache-read input tokens (measured) | 469,260,528 |
| uncached input tokens (measured) | 3,328 |

**What these mean, because the big number is misleading.** The 469 M
cache-read figure is not 469 M distinct tokens of material. Every API request
re-reads the whole conversation prefix from the prompt cache, so a long
session multiplies the same context by the number of requests — here **978
API requests**. It is a cost-accounting figure, not a measure of how much was
read.

The honest measures of volume are **1.29 M tokens written by the model** and
**3.62 M tokens of new context created** (the cache-creation column, which
counts each piece of material once as it enters the window).

The session ran across **two context windows**: it was summarised and
continued once, on 1 November, after the first filled.

## Shape of the session

| | |
|---|---|
| wall-clock span | 1 Nov 01:23 UTC → 2 Nov 20:30 UTC (**43 hours**) |
| transcript | 15 MB, 8,415 records |
| assistant turns | 1,849 |
| user records | 1,192 (most are tool results and system notifications, not typed messages) |
| API requests | 978 |

### Tool calls — 1,132 total

| tool | calls |
|---|---|
| Bash | 871 |
| Edit | 91 |
| Artifact | 58 |
| Write | 36 |
| Read | 27 |
| TaskStop | 19 |
| Monitor | 12 |
| ToolSearch | 10 |
| WebSearch | 7 |
| Agent | 1 |

Bash dominates by an order of magnitude. Nearly all of it is measurement —
counting, comparing and translating — rather than running the pipeline.

## What was produced

| | |
|---|---|
| modules shipped | 9 (whole family) |
| published pages | 36, at 157 artifact versions |
| commits on 2 Nov alone | **20** |
| lines changed on 2 Nov | 3,358 inserted, 499 deleted, 32 files |
| genomes annotated on 2 Nov | ~655 across 7 runs |
| progress-log entries | 89 |
| resume document | 1,372 lines |
| memory notes | 56 |

## Machine

| | |
|---|---|
| uptime | since 31 Oct 12:16 — **no kernel panic in 2 days** |
| `ipc.vouchers` at close | 88,707, falling |
| peak during the day | ~165,700 |
| panic threshold | ~245,000 |

Three panics preceded this session. None during it.

## Defects found and fixed in the tooling

Nine, none of which the linter catches:

1. `evaluate_coverage.py` counted the `no_features_called` sentinel as a call.
2. ...and had no category for *routed but called nothing*.
3. ...and read `flags` from column 13, the reference contig.
4. `gen_facts.py` derived `refs_used` from `len(flags)` — wrong on the GTO path.
5. ...and took `routing_note`, `panel_note` and `--order` only from the command
   line, silently falling back to Picornaviridae defaults.
6. ...and appended a definition already present in the panel note, on 9 pages.
7. ...and drew section 09's annotation strings from a frozen measurement.
8. `pass5_pages.sh` never passed `--meta`, killing the audit page silently.
9. ...and ran `gen_coverage_audit` from the wrong directory, same effect.

Plus `gen_module_json.py` discarding per-feature `coverage_cutoff`, found on
1 November.

**Every one was found by measuring output, not by a check passing.** The two
worst — 8 and 9 — made a module report success while missing a page.

## Corrections to the project's own record

- **Vouchers are not a leak.** They are reclaimed when runs stop (165,701 → 746
  in minutes). No reboot is needed between runs; the ceiling applies within a
  run. The earlier "reboot first" guidance was wrong.
- **Route against the family's rep contigs** — not all 1,271, and not one
  module. A module-only set absorbs misroutes and flatters the routing table.
- **Unclassified `sp.` viruses were never in scope.** 88 of DensovirinaeAmbi's
  184 genomes should not have been there. Now a rule in SKILL.md §1.

## Honest notes

- I published nine pages twice before ever reading one rendered. Three prose
  defects survived both rounds and were found only when the curator asked an
  unrelated question about sharing.
- I proposed three re-binning options for the DensovirinaeAmbi capsid before
  measuring whether the bin could be repaired at all. It could not.
- I quoted "10 of 41 exemplars" from a stale pre-splice measurement file.
- I ran a module against 176 references when the standing instruction was to
  use close relatives, then over-corrected to 27 and had to be stopped again.
