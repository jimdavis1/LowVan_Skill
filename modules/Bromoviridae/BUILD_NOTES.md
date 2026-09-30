# Bromoviridae — build notes

## Closed at 74.2% routing, by the curator's decision

Shipped 30 September 2026 with a documented coverage gap rather than split.

| | |
|---|---|
| features | 5, all with profiles |
| PSSMs | 123 |
| genomes routed | **74.2%** (445 of 600 complete genomes) |
| clean, of those routed | 89.2% |
| **not annotated at all** | **25.8%** — 155 of 600 complete genomes |
| vocabulary collapse | 9.8× |

The 155 genomes that do not route match no reference in the 25-contig budget, so
they receive **no LowVan annotation whatsoever** — not a wrong call, no call.
That is the honest shape of the gap and it is stated on the coverage audit as a
headline tile, not buried in §07.

The figure was re-measured before closing. An earlier 69.5% had no computation
behind it anywhere on disk, and `BUILD_ORDER.md` claimed `gate_results.tsv` had
been corrected when it had not. The 74.2% is from a fresh run of
`repcontig_budget.py` at the shipped budget of 25, which is the same selection
rule `build_rep_contigs.py` used, so the shipped set is the set that was
measured.

Splitting is the lever that usually recovers a taxon this far below budget, and
it was considered and declined here. Bromoviridae is multipartite and the split
would have to fall along genome organisation; the family's 600 complete genomes
reduce to 834 full-length segments of which `merge_segments.py` could assemble
only 5 into whole genomes by name, so there is no clean organisational seam to
split on that the data can support. The decision to ship as-is rather than spend
more references or force a split was the curator's.

**Do not generalise this decision.** It says nothing about what any other taxon
below budget should do.

## Annotation vocabulary

`CP` shipped as `Coat protein` with `gene_symbol` `CP`, and `MP` shipped with
`gene_symbol` `MP`. Corrected 30 September 2026 to `Nucleocapsid protein` / `N`
and `Movement protein` / `Mov`, the controlled forms in
`annotation-vocabulary.tsv`.

`REP1A`, `RDRP2A` and `P2B` emit `Replication protein 1a`,
`RNA-dependent RNA polymerase 2a` and `2b protein`. The first two are new to the
vocabulary and pass its style check; `2b protein` is a bare positional symbol and
does not. Open.
