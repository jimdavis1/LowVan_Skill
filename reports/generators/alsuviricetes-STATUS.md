# Alsuviricetes module status

Last updated 17 September 2026. The repository of record is
`LowVan_Skill/` (github.com/jimdavis1/LowVan_Skill). `Viral_Annotation/` is a
runtime target only and has no git remote.

## Shipped in the kit

| module | features | profiles | rep contigs | routing | quality | reports |
|---|---|---|---|---|---|---|
| Alpharhabdovirinae | 43 | 529 | 10 | — | — | registry, audit |
| Betarhabdovirinae | 14 | 340 | 8 | — | — | — |
| Novirhabdovirus | 8 | 21 | 2 | — | — | — |
| Dichorhavirus | 8 | 26 | 1 | — | — | — |
| Togaviridae | 14 | 228 | 9 | — | — | — |
| Matonaviridae | 7 | 16 | 4 | — | — | all four |
| Hepeviridae | 4 | 63 | 24 | 622/669 | **89.1%** | all four |
| **Tobamovirus** | 4 | 98 | 25 | 99/102 | **89.4%** | **all four** |
| **Betaflexiviridae_MP** | 5 | 157 | 25 | 407/409 | **88.5%** | **all four** |

Every module's installed profile count now equals its shipped alignment count.

Quality is reported as **genuinely clean**: no genome, contig *or* feature
flag. `run_gto_eval.py` decides good-versus-poor on genome and contig flags
alone, so its headline is 1-4 points higher than these figures. Scored
consistently the three modules converge at 88-89%.

**Hepeviridae was checked against the denominator bug and is unaffected.**
All 47 of its unscored genomes are *unrouted*, with no `viral_family` field,
so they fail a different and correct check and are already counted in the
622/669 routing figure. Its denominator of 622 is exactly the routed set.
Tobamovirus was inflated (95.5% reported, 90.4% true) because its five had
routed and then cleared no profile, which is a different situation.

## Next, in order

Chosen on the two rules the user set: easiest first, and weight crop/human
impact. A module is only started once its rep contigs are known to fit the
25-reference budget.

| module | measured routing | notes |
|---|---|---|
| Alphaflexiviridae | 84.9% at budget | Potexvirus + Allexivirus + Lolavirus as one module, per PARTITIONING.md. one module, 7 features and 123 profiles built; collections triaged and Mandarivirus recovered by homology three times under a Potexvirus label. Splitting would reach ~95% for 50 references |
| Betaflexiviridae_TGB | 81.0% | pome and stone fruit; the TGB half of the family already split |
| Capillovirus | 100% | apple stem grooving; needs a mature peptide (CP inside the polyprotein) |
| Bromoviridae | unmeasured | **cucumber mosaic virus**; 3 segments, needs per-segment clustering |
| Alphaflexiviridae_noTGB | unmeasured | 42 genomes, fungal; may not ship |

Unmeasured beyond that: Crinivirus, Virgaviridae minors, Kitaviridae,
Mayoviridae, Benyviridae.

## Tabled on the rep-contig budget

| module | coverage at 25 | why |
|---|---|---|
| Endornaviridae | ~30% | 141 species, 141 clusters. One ORF, so trivially easy to build and impossible to route |
| Closteroviridae (mono) | 74.6% | try splitting before retrying |
| Tymoviridae (lumped) | 74.6% | try splitting before retrying |

## Deferred by the user

**MERHA_L_SPLICED.** Real: Kuwata 2011 (PMID 21507977) describes a 76-nt
GU-AG intron in the citrus tristeza-related CTRV L gene at 8648-8723. All 13
CTRV genomes are near-complete. Before it ships, the splice call must be shown
to *replace* the truncated L call rather than duplicate it.

## Open questions for the curator

- `coverage_cutoff` does not measure coverage. Reported, deliberately left
  unpatched pending a decision.
- ~30 hand-made modules in the runtime directory are not in the kit, and 15 of
  them have no alignments at all, so their profiles cannot be rebuilt.
- Every citation in the Claude-made modules is model-proposed and lives in
  `PMID_claude_generated`. None has been read and confirmed.
- Tobamovirus: 13 of 88 unassigned features are canonical names ("30K
  protein", "29K protein", "17 kDa protein") the binning rules should catch.
