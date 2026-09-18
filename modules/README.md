# Built modules

Every module built with this kit. Eight of them, 102 features, 1,319
alignments, 83 rep contigs.

| module | features | alignments | rep contigs | special |
|---|---|---|---|---|
| `Alpharhabdovirinae` | 43 | 529 | 10 | `MERHA_L_SPLICED` — reference set never built, see SPECIAL_FEATURES.md |
| `Betarhabdovirinae` | 14 | 340 | 8 | |
| `Novirhabdovirus` | 8 | 21 | 2 | |
| `Dichorhavirus` | 8 | 26 | 1 | |
| `Togaviridae` | 14 | 228 | 9 | `TF` transcript_edit, 234 references included |
| `Matonaviridae` | 7 | 16 | 4 | |
| `Hepeviridae` | 4 | 61 | 24 | |
| `Tobamovirus` | 4 | 98 | 25 | `REP183` readthrough (`internal_stop`) |

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
