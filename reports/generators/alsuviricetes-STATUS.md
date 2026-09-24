# Alsuviricetes module status

Last updated 24 September 2026. The repository of record is
`LowVan_Skill/` (github.com/jimdavis1/LowVan_Skill). `Viral_Annotation/` is a
runtime target only and has no git remote.

## Shipped in the kit

Sixteen modules, 178 features, 3,144 alignments, 265 rep contigs,
3 transcript-edited reference sets.

| module | features | profiles | rep contigs | routing | genuinely clean | artifacts |
|---|---|---|---|---|---|---|
| Alpharhabdovirinae | 43 | 529 | 9 | — | — | registry, audit |
| Betarhabdovirinae | 14 | 340 | 8 | — | — | — |
| Novirhabdovirus | 8 | 21 | 2 | — | — | — |
| Dichorhavirus | 8 | 26 | 1 | — | — | — |
| Merhavirus | 8 | 16 | 8 | — | — | — |
| Togaviridae | 14 | 228 | 9 | — | — | — |
| Matonaviridae | 7 | 16 | 4 | — | — | all four |
| Hepeviridae | 4 | 63 | 24 | 93.0% | **89.1%** | all four |
| Tobamovirus | 4 | 98 | 25 | 97.1% | **89.4%** | all four |
| Trivirinae | 5 | 157 | 25 | 99.5% | **88.5%** | all four |
| **Alphaflexiviridae** | 6 | 259 | 25 | 90.1% | **92.4%** | **all four** |
| **Allexivirus** | 7 | 76 | 25 | 99.0% | **89.3%** | **all four** |
| **Orthoflavivirus** | 16 | 337 | 25 | 98.3% | **62.0%** | **all four** |
| **Hepaciviridae** | 12 | 808 | 25 | 94.7% | **53.8%** | **all four** |
| **Pestiviridae** | 13 | 103 | 25 | 98.8% | **28.9%** | **all four** |
| **Pegivirus** | 9 | 67 | 25 | 92.8% | **20.8%** | **all four** |

Quality is **genuinely clean**: no genome, contig *or* feature flag.
`run_gto_eval.py` decides good-versus-poor on genome and contig flags alone,
so its headline runs a few points higher.

**Pestiviridae declares one feature that has no PSSM at all.** NS2 is cut at
the Arg1589-Gly1590 site fixed by protein sequencing (PMID 12827468) rather
than modelled, because BV-BRC holds 69 NS2 records of which 34 are the right
length and 24 of those end in the wrong residue. It calls on 111 genomes with
110 ending in the correct Arg. The same was tried on Npro and **rejected on
measurement** — 58% ending in the conserved Cys168 against the PSSM's 96% —
because a derived end inherits its anchor's boundary and the capsid profile
does not place that junction reliably.

**Hepaciviridae's 53.8% is a species-capped figure and means almost nothing on
its own.** The panel gives every species equal weight, so 44 genomes of rare
divergent hepaciviruses count for as much as 51,573 *Hepacivirus C*. On HCV
itself the eleven PSSM features call at 96.8&ndash;98.4% and F at 87.4%. The
features that read low &mdash; E1, E2, NS2, p7 &mdash; do so because ten
lineages have a polyprotein annotation and no mature peptides at all, which is
0.07% of the database. See §05 and §07 of the audit.

**The two newest rows have the hardest denominators, and the low numbers are
not what they look like.** Both declare many essential features over one
polyprotein, so a genuinely-clean verdict requires every one of them at once;
Alphaflexiviridae above declares six features in total. Orthoflavivirus's
per-feature rates are 91.8–99.8%. Neither row is comparable to the five above
it.

**Denominators are not comparable.** Hepeviridae (n=622), Alphaflexiviridae
(n=170) and Allexivirus (n=103) were scored on full routed sets. Tobamovirus
(n=94) and Trivirinae (n=96) rest on ~95-genome samples whose 95% confidence
intervals are about eleven points wide, so the apparent convergence at 88-89%
across those two was never measurable.

## Allexivirus was split out of Alphaflexiviridae

One feature needed opposite settings in two genera. Lezzhov et al. 2015
(J Gen Virol 96:3159-64, PMID 26296665) showed by site-directed mutagenesis
that shallot virus X initiates TGB3 at a **CUG**, translated by leaky scanning
from a bicistronic template shared with TGB2. The coding sequence is conserved
across the genus; only the initiator is non-canonical.

`upstream_ext` is now set per feature from measurement rather than policy.
`scan_to_met_start` was instrumented across all 273 annotated genomes and
every firing recorded as finding a Met or running to the previous in-frame
stop:

```
                Allexivirus      Potexvirus
  feature     found / ran-to    found / ran-to    ext
  CP             12 / 2            26 / 6          1 / 1
  TGB2            0 / 0            13 / 4          1 / 1
  ALLEXI_40K     19 / 4             - / -          1 / -
  NABP            0 / 0             1 / 0          1 / 1
  REP             1 / 4             3 / 6          0 / 0
  TGB1            0 / 18            7 / 18         0 / 0
  TGB3            2 / 190          38 / 10         0 / 1
```

TGB3 is the row that cannot be reconciled in one module. Two results would
have been missed by setting this from policy: **TGB1 and REP fail in both
genera** and are 0 in both modules.

Measured effect of the split, same genomes and tools before and after:

| | combined | after |
|---|---|---|
| genuinely clean | 227/271 (83.8%) | 249/273 (**91.2%**) |
| `Feature is too long` | 30 | **4** |
| TGB3 calls not at a Met | 100/266 (37.6%) | 5/164 (**3.0%**) in Alphaflexiviridae |

The Allexivirus TGB3 non-Met rate is 96% and is **expected**, not a defect:
there is no AUG. Its 5' coordinate is approximate because the CUG cannot be
located mechanically — 65% of genomes have some in-frame CTG upstream at 1-47
codons with no consistent offset, and 35% have none before a stop.

## Defects found and fixed this unit

- **`install_module.py` never pruned stale alignments.** It pruned stale
  profiles inside a rebuilt feature but had no equivalent for alignments, so
  every rebuild left the old ones behind. Alphaflexiviridae shipped 276
  profiles against 289 alignments. Fixed; backfilled Hepeviridae 63/61 -> 63/63
  and Tobamovirus 98/100 -> 98/98.
- **223 rhabdovirus alignments were ragged** — raw unaligned sequence sets with
  no gap characters, 152 in Alpharhabdovirinae and 71 in Betarhabdovirinae.
  The kit ships alignments and rebuilds profiles from them, so those profiles
  could not be rebuilt at all. Realigned with MAFFT; verified end to end by
  regenerating 529/529 Alpharhabdovirinae profiles from kit files alone. This
  closes the "223 of 1474" item.
- **`loN` profile names removed.** The prefix leaked a build-pipeline detail
  into `family_assignments` in every annotated genome. 684 profiles renumbered
  to plain integers across twelve modules; generator fixed.
- **`renumber_leftover_profiles.py` missed the working directory's `pssms/`.**
  `install_module.py` installs from there, so the first reinstall after
  renumbering put 124 `loN` profiles back against integer-named alignments.
  Fixed and repaired across six working directories (182 profiles).

## Deliberately not changed

**`scan_to_met_start` runaway extension.** When no upstream Met exists the
function returns the codon after the previous in-frame stop rather than
declining, so the CDS runs to the edge of the upstream ORF. A patch was
written and measured (+18 points genuinely-clean on a 148-genome subset, a
no-op where profiles are adequate) and then **reverted at the curator's
direction**: the original behaviour stands. The per-feature `upstream_ext`
settings above are the containment.

**`coverage_cutoff` does not measure coverage.** Computes the alignment
against its own query span, so it is 1.0 for any gapless alignment and cannot
reject a partial profile match. Reported, deliberately unpatched.

## Next

| module | measured routing | notes |
|---|---|---|
| Quinvirinae | 81.0% | **no working directory exists** — this is a from-scratch build, not a rename |
| Capillovirus | 100% | apple stem grooving; CP inside the ORF1 polyprotein, needs a mature peptide |
| Bromoviridae | unmeasured | cucumber mosaic virus; 3 segments, needs per-segment clustering |
| Botrexvirus_Platypuvirus_Sclerodarnavirus | unmeasured | 42 genomes, fungal; may not ship |

Unmeasured beyond that: Crinivirus, Virgaviridae minors, Kitaviridae,
Mayoviridae, Benyviridae.

## Tabled on the rep-contig budget

| module | coverage at 25 | why |
|---|---|---|
| Endornaviridae | **71.8%** | 202 genomes in 69 clusters. The curve does not bend — 41% at 5 refs, 72% at 25, 91% at 50 — because nearly every genome is its own cluster. **No split can help**: one ORF, no capsid, uniformly divergent |
| Closteroviridae | **77.0%** | six genera in one module; **not yet split**, and the split is the lever to try |
| Quinvirinae | **81.0%** | 1,694 genomes in 142 clusters; **not yet split** |

Re-measured 23 September with one fixed criterion — records at ≥6,000 nt for a
monopartite taxon, merged by `genome_name` for a multipartite one. It
reproduces the independent 22 September figures (Closteroviridae 77.0% against
76.6%, Quinvirinae 81.0% against 81.1%), which is why the numbers above
supersede the older ones.

**Tymoviridae is no longer tabled.** Lumped it gated 74.6%; gated per genus its
three members all clear — Marafivirus 100%, Tymovirus 91.2%, Maculavirus 100%.
Same genomes, same budget. The lumped figure was measuring the module boundary,
not the taxon, and it is the precedent for splitting the two above.

## Open questions for the curator

- Every citation in the Claude-made modules is model-proposed and lives in
  `PMID_claude_generated`. None has been read and confirmed, including the
  Lezzhov paper that justifies the Allexivirus split.
- **Alphaflexiviridae NABP is trained backwards** — 338:54 Allexivirus to
  Potexvirus, inherited from before the split, for genomes that now route
  elsewhere, and called on 5 of 170. Rebuilding it from Potexvirus sequences
  alone would be more honest.
- **Allexivirus TGB3 rests on three sequences** from one virus. Two more exist
  in the collection (Blackberry calico, Garlic virus C AS16) and never
  clustered. Recovering them is the cheapest improvement available.
- **Lolavirus is named but unmodelled** — 10 genomes, four sequences per
  feature, no rep contig, riding on Potexvirus homology.
- `copy_num` in `viral_genome_quality.pl` couples "expected copy number" to
  "required present", so a genus-restricted accessory can have neither.
- ~30 hand-made modules in the runtime directory are not in the kit.
