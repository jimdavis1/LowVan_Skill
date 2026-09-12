# Vendored file manifest

Source: BV-BRC desktop bundle version 1.046, `/Applications/BV-BRC.app/deployment/`.
Upstream: <https://github.com/BV-BRC/BV-BRC-CLI/releases> and <https://github.com/TheSEED>.

Licence: SEED Toolkit Public License. Copyright (c) 2003-2013 University of
Chicago and the Fellowship for Interpretation of Genomes. Six of these files
carry no per-file header upstream; the licence covers them all, and this
manifest is the provenance record rather than a header added by us -- the
files are byte-for-byte as shipped.

To refresh: copy the same relative paths out of a newer deployment and
re-run the closure check in `vendor/README.md`.

| file | sha256 (first 16) | bytes |
|---|---|---|
| `lib/BBasicLocation.pm` | `5938c7d6e30c5370` | 15438 |
| `lib/BasicLocation.pm` | `353eef96d9b33ab2` | 25160 |
| `lib/Bio/KBase/GenomeAnnotation/Client.pm` | `cce12506c42c3c1c` | 2387022 |
| `lib/Bio/KBase/GenomeAnnotation/CmdHelper.pm` | `17bfa8b40e477895` | 8570 |
| `lib/FBasicLocation.pm` | `265cb9e58bef9b0d` | 15392 |
| `lib/GenomeTypeObject.pm` | `e270cc82545d2da8` | 53443 |
| `lib/IDclient.pm` | `91838d6ca44ac8a0` | 2338 |
| `lib/P3AuthToken.pm` | `ea9b48ae4d1df563` | 3318 |
| `lib/SeedUtils.pm` | `30db2be41e7ca678` | 75374 |
| `lib/Sim.pm` | `465ab40fe913752c` | 11661 |
| `lib/gjoseqlib.pm` | `8579df9cdd997da1` | 91069 |
| `plbin/rast-create-genome.pl` | `1604a93d8f11d638` | 3162 |

`bin/rast-create-genome` is **not** upstream: it is the wrapper this
repository adds, MIT like the rest of our own code.
