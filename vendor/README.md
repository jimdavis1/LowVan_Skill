# Vendored third-party code

Open-source code from other projects, redistributed here so a clone of this
repository works without hunting down separate installs. All of it is licensed
for redistribution. Each file keeps its own licence and copyright; the MIT
licence at the repository root covers our code, not these.

## `bv-brc/`

The minimum set of BV-BRC / SEED modules needed to build and score Genome Type
Objects, so that a clone of this repository can run the GTO path with **no
BV-BRC installation at all**. It is a dependency closure computed from three
entry points — `GenomeTypeObject.pm`, `IDclient.pm` and `rast-create-genome.pl`
— not a copy of the whole distribution.

```
lib/GenomeTypeObject.pm                          the GTO itself
lib/IDclient.pm                                  mints feature ids
lib/P3DataAPI.pm  lib/P3Utils.pm                 BV-BRC data API, for the dump
lib/BlastInterface.pm  lib/gjoparseblast.pm      BLAST wrappers used by the
lib/gjoalignment.pm  lib/AminoAcidMatrix.pm        clustering pipeline
lib/SeedUtils.pm  lib/gjoseqlib.pm               sequence and genetic-code utilities
lib/SeedAware.pm  lib/Sim.pm  lib/P3AuthToken.pm
lib/BasicLocation.pm  lib/BBasicLocation.pm
lib/FBasicLocation.pm                            feature coordinates
lib/Bio/KBase/GenomeAnnotation/Client.pm         service client
lib/Bio/KBase/GenomeAnnotation/CmdHelper.pm
plbin/rast-create-genome.pl                      creates a GTO, issuing its id
bin/rast-create-genome                           wrapper (ours, see below)
```

Four closures, from `GenomeTypeObject.pm`, `IDclient.pm`, `P3DataAPI.pm` and
`BlastInterface.pm`. **This covers the whole toolchain** — the annotator, the
clustering pipeline, the BV-BRC dump and the GTO path all run against
`vendor/bv-brc/lib` alone. No BV-BRC install and no `seed_gjo` checkout.

**Licence.** The SEED Toolkit Public License, which **permits redistribution** —
that is what makes bundling these possible. Copyright (c) 2003-2013 University
of Chicago and the Fellowship for Interpretation of Genomes. Full text at
<http://www.theseed.org/LICENSE.TXT>; upstream at
<https://github.com/BV-BRC/BV-BRC-CLI> and <https://github.com/TheSEED>. Keep
the per-file notices and this directory's attribution intact and you are within
its terms.

**The wrapper is ours, not upstream's.** The BV-BRC bundle ships
`plbin/rast-create-genome.pl` but no `bin/` wrapper for it — only
`rast-create-genome-from-RAST` has one. `bin/rast-create-genome` is that missing
wrapper, written to resolve relative to wherever this repository was cloned
rather than to a fixed install prefix. It honours `$LOWVAN_PERL` if you need to
point it at a particular interpreter.

## What is still needed from CPAN

The vendored tree resolves every BV-BRC module. Two CPAN modules are not
vendorable and must be present in whichever perl you use:

```
File::Slurp     hard dependency of GenomeTypeObject and CmdHelper
Data::UUID      (or UUID) -- GenomeTypeObject mints feature ids with it
```

```bash
cpanm File::Slurp Data::UUID
```

`make_gto.py` checks for both before doing any work and prints this command if
either is missing. Without the check, a missing UUID module is silent until
`create_uuid()` is reached and the run dies partway through with
`No UUID generator found`.

Everything else these modules load — `FIG_Config`, `Config::Simple`,
`Data::Structure::Util`, `FIGO` — is guarded by `eval` or sits in a code path
this workflow never enters, and none of it needs installing.

## Verifying it works standalone

```bash
./check-environment.sh          # checks binaries, CPAN modules, and compiles

env -u PERL5LIB vendor/bv-brc/bin/rast-create-genome \
    --scientific-name "Test virus" --domain Viruses --genetic-code 1 \
    --ncbi-taxonomy-id 11292 --contigs some.fna -o test.gto
```

`-u PERL5LIB` is the point: it proves the vendored tree alone is sufficient.
A genome id is issued by the ID server, so this needs network access.
