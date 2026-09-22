# Built modules

Every module built with this kit. Fourteen of them, 156 features, 2,803
alignments, 215 rep contigs, two transcript-edited feature sets.

| module | features | alignments | rep contigs | special |
|---|---|---|---|---|
| `Alpharhabdovirinae` | 43 | 529 | 9 | `MERHA_L_SPLICED` — reference set never built, see SPECIAL_FEATURES.md |
| `Betarhabdovirinae` | 14 | 340 | 8 | |
| `Novirhabdovirus` | 8 | 21 | 2 | |
| `Dichorhavirus` | 8 | 26 | 1 | |
| `Merhavirus` | 8 | 16 | 8 | `L_SPLICED` splice, 12 references |
| `Togaviridae` | 14 | 228 | 9 | `TF` transcript_edit, 234 references included |
| `Matonaviridae` | 7 | 16 | 4 | |
| `Hepeviridae` | 4 | 63 | 24 | |
| `Tobamovirus` | 4 | 98 | 25 | `REP183` readthrough (`internal_stop`) |
| `Trivirinae` | 5 | 157 | 25 | |
| `Alphaflexiviridae` | 6 | 259 | 25 | Potexvirus + Lolavirus |
| `Allexivirus` | 7 | 76 | 25 | TGB3 `upstream_ext: 0`, non-AUG start |
| `Hepaciviridae` | 12 | 637 | 25 | hepatitis C; `F` `transcript_edit` (+1 frameshift, 844 refs) — the first module needing the deletion direction |
| `Orthoflavivirus` | 16 | 337 | 25 | 15-product polyprotein; ends validated against cleavage chemistry. `NS1P` `transcript_edit`, 104 references — cannot distinguish the SA14-14-2 vaccine lineage, see the audit |

`Alphaflexiviridae` and `Allexivirus` are one taxon split in two. They are
listed separately because the split is the point: Allexivirus TGB3 initiates
at a CUG rather than an AUG (Lezzhov 2015, PMID 26296665), so it needs
`upstream_ext: 0` where the same feature in Potexvirus needs `1`. A single
module cannot hold both settings. Every feature's `upstream_ext` in both
modules is set from a measurement of how often the Met scan succeeds, not
from policy — see either module's build script.

Leftover profiles are numbered as plain integers continuing the main pass.
Older builds named them `lo1`, `lo2` …; that prefix leaked a build-pipeline
detail into `family_assignments` in every annotated genome and has been
removed. `skill/scripts/renumber_leftover_profiles.py` performs the migration
and must move the working directory's `pssms/` as well as the repository's,
or the next install puts the old names back.

The four Rhabdoviridae modules, Togaviridae and Matonaviridae were built before
the Alsuviricetes work and imported here afterwards: the kit had always carried
the **lessons** from those builds — Rhabdoviridae is cited 65 times across
SKILL.md and the references, Togaviridae 12 — but never the **products**. The
reasoning was captured and the modules were not, which meant a reader could
learn how to build one but could not obtain one already built.


One directory per finished module, carrying everything needed to install it
into a `Viral_Annotation` checkout and everything needed to audit how it was
made.

```
modules/<Module>/
  <Module>_Viral_PSSM.json     the module block: features, cutoffs, lengths, citations
  Rep-Contigs/                 routing references + close_genomes.json
  PSSM-Alignments/<FEAT>/corrected_alis/   the curated alignments every profile
                                           was built from
```

## Why the PSSMs are not here

They are **derived**, and they are large: Tobamovirus is 26 MB of profiles
against 2.1 MB of alignments, a tenfold difference, and across a whole class
that is a gigabyte of files that no human will ever read. The alignments are
the provenance — they are text, they diff, and they are what a curator edits.

Rebuild the profiles from them:

Point a working directory at them and generate:

```bash
mkdir -p mywork/Alignments
cp -R modules/<Module>/PSSM-Alignments mywork/Alignments/<Module>
cp modules/<Module>/<Module>_Viral_PSSM.json mywork/
cp -R modules/<Module>/Rep-Contigs mywork/

python3 skill/scripts/pssms_from_alignments.py --workdir mywork \
        --module <Module> --write
python3 skill/scripts/install_module.py --workdir mywork \
        --repo /path/to/Viral_Annotation --check
```

The alignments ship under `<FEAT>/corrected_alis/`, which is the layout the
tools already expect, so nothing has to be moved first. Verified end to end
from these files alone: 98 Tobamovirus profiles regenerated and
`install_module.py --check` reported 0 problems.

**Not `rebuild_pssms.py`** — that one walks the `pssms/` directory and refreshes
each profile from its alignment, so it detects drift but cannot create a
profile that is absent. Pointing it at an alignments-only tree regenerates
nothing and reports success, which is how this README came to claim something
untrue in its first draft.

`pssms_from_alignments.py` walks the alignments instead and writes one profile
per alignment in the same format, so the two tools agree. Verified two ways. Deleting all 27 Tobamovirus CP profiles and regenerating
them reproduces all 27 **byte-identically** — both were made by the same
psiblast path — and `rebuild_pssms.py` then reports `0 stale` across the
module. Against profiles the *pipeline* built, regeneration reproduces the
**score matrix exactly** (Novirhabdovirus: 21 of 21, integer for integer) but
not the bytes: `lambdaUngapped` differs in its last digit between
`BlastInterface::alignment_to_pssm` and the psiblast CLI. That is floating
point and has no effect on scoring. Compare score matrices, not checksums.

## Installing one

```bash
python3 skill/scripts/install_module.py --workdir <module workdir> \
        --repo /path/to/Viral_Annotation --check      # then without --check
```

## Before you build anything against a Viral_Annotation checkout

**Copy `build/*` and `annotate/*` from this kit over it first.** The kit's own
directories are the source of truth. `patches/` is a diff against pristine
upstream and drifts behind them as work accumulates, and
`CEPI-dxkb/Viral_Annotation` has not received the annotator work at all — it
still has no `internal_stop` support, so any module relying on a readthrough
will be silently cropped. That repo gets updated from this one when the
maintainer decides the work is ready, not incrementally. This cost one rediscovered bug and one module built
against an older pipeline before it was noticed.
