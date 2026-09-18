# Built modules

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
per alignment in the same format, so the two tools agree. Verified: deleting
all 27 Tobamovirus CP profiles and regenerating them from the alignments
reproduces all 27 **byte-identically**, and `rebuild_pssms.py` then reports
`0 stale` across the whole module.

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
