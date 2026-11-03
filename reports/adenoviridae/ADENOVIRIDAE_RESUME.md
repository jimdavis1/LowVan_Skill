# Adenoviridae — resume here

**Written 3 Nov 2026.** If the machine crashed, start at the first unticked
box. `PROGRESS.md` in this directory is the append-only log; its last line is
where work stopped.

## Where everything is

| | |
|---|---|
| scratch | `~/lowvan_scratch/adeno` (this directory) |
| kit | `$LOWVAN_KIT` = `.../CEPI/Alsuviricetes/LowVan_Skill` |
| env | `source .../CEPI/Retroviridae/env.sh` |
| repo to push to | `jimdavis1/LowVan_Skill` — **and only that one** |
| cost meter | `scripts/tally_cost.sh`, baseline = transcript record 8607 |

**Do not write to `Alsuviricetes/Viral_Annotation`** (perl libs only, curator's
instruction 7 Oct) or to `~/Viral_Annotation`.

## Scope, settled by the curator 3 Nov

- Build **Mastadenovirus** (4,244 genomes) and **Aviadenovirus** (447).
- **Measure then decide** on the next: Siadenovirus (61) and Barthadenovirus
  (55) are both candidates, roughly equal in size.
- Ichtadenovirus (2) and Testadenovirus (0 near-complete) are out of scope.
- Feature depth: **every known protein and mature peptide, including signal
  peptides.**

## Checklist

- [x] 1. BV-BRC dump. `dump.sh`, prefix `A`, `--min-len 20000`.
      5,015 genomes / 155,187 feature rows / 21,947 distinct sequences.
      `check_dump.py` PASS.
- [x] 1b. SKILL §1 intake rule: 139 `sp.` genomes dropped to
      `intake/DROPPED_sp.tsv`; 67 genus-blank held out to `intake/NO_GENUS.tsv`.
- [x] 2. Triage rules — `build_cols.py`. **The design decision to preserve:
      masses and ORF numbers do not bind by name.** Both go to a pool
      resolved by blastp (>=80% id / >=85% qcov) and then mmseqs clustering.
      Writing a rule for "34 kDa protein" by eye is the trap the skill names.
- [x] 2b. Literature probe list — `probes.json`. Every PMID resolved against
      PubMed by `scripts/pm.py` before it was written down.
- [ ] 3. Run `check_literature_gaps.py --verify-pmids` against the pass-3
      collections and resolve every `UNBOUND`.
- [ ] 4. Name the surviving UNC clusters, or record why they stay UNC.
- [ ] 5. Length windows per feature, read per species out of
      `<FEAT>.outliers.fasta` — never a pooled median.
- [ ] 6. Mature peptides. See `MATURE_PEPTIDES.md` — they are NOT in the dump
      and must be built from published adenain cleavage sites.
- [ ] 7. Alignments, PSSMs, rep contigs, module JSON.
- [ ] 8. Install, annotate, coverage. **Route against the FAMILY's rep contigs
      (SKILL §11), not the whole repository and not one module.**
- [ ] 9. Four pages per module.

## Defects already found and fixed (do not re-introduce)

1. Rules enumerated `( protein)?` by hand and missed "E3 14.7K protein" — 79
   features became accessory cluster UNC7. Fixed with `normalise()`, which
   strips the suffix and the hedges off the string once. **"hypothetical" is
   deliberately not stripped** — it would leave the bare word "protein" and
   `POOL` would stop matching 1,115 Aviadenovirus features.
2. `consolidate()` emptied `E1A_243R` into `E1A`. 243R is 289R minus an
   internal 46-residue block, so homology can never separate them. Fixed with
   the `PARALOG` exemption, which also protects FIBER1/FIBER2.
3. The accessory floor counted features only, so 28 Mastadenovirus "clusters"
   were a single sequence recurring in five genomes. Fixed: `ACC_SEQ_FLOOR 3`
   **and** `ACC_FLOOR 5`.
4. Clustered sequences were never removed from `pool`, so "still unplaced"
   and `UNASSIGNED_TRACKING.tsv` both overstated the unfinished work — ORF20
   sat there at 324 features while UNC2 was built from 72 of its sequences.
5. A reporting artefact, not a defect: "3,356 genomes" is genomes
   contributing a bound feature, against 3,362 Mastadenovirus genomes that
   have any annotated protein at all. The other 882 carry **no CDS in
   BV-BRC**. Per-genome binding is a median of 100%, mean 97.7%.

## Biology the triage found that the first rule set had missed

- **Mastadenovirus has two fiber genes** in the enteric species (HAdV-F40/41).
  PMID 7975882, "Human adenovirus type 41 contains two fibers". Pass 1 left
  them pooled as UNC2 (median 555 aa) and UNC5 (median 368 aa).
- **E4 ORF5** — UNC1 in pass 1, 59 sequences at a median of 151 aa.
- **E3 6.7K** — found by the probe pass, PMID 15507617.
