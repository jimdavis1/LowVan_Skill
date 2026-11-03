# Adenoviridae — progress log

Started 3 Nov 2026 01:46 UTC. Append-only; newest at the bottom. Every entry
is timestamped so a crash can be resumed from the last line.

Scope settled by the curator 3 Nov: build **Mastadenovirus** and
**Aviadenovirus**; measure the next-largest genus and decide. Feature depth:
**every known protein and mature peptide, including signal peptides.**

| | |
|---|---|
| cost meter baseline | transcript record 8607, 2026-11-03T01:46Z |
| tally | `scripts/tally_cost.sh` |

---
- 01:48Z Measured the family holding. 43,636 records; 5,351 >=5 kb; 5,015 >=20 kb.
- 01:48Z Genus holdings >=20 kb: Mastadenovirus 4,362 | Aviadenovirus 450 | blank 80 | Siadenovirus 66 | Barthadenovirus 55 | Ichtadenovirus 2 | Testadenovirus 0.
- 01:48Z CORRECTION to the scoping report: Atadenovirus is NOT absent. BV-BRC files it as Barthadenovirus (duck/ovine/bovine atadenovirus A-D, snake, lizard, bearded dragon), 55 near-complete. Four candidate genera, not three.
- 01:48Z Splicing feasibility settled on NC_001405 (human adenovirus 2): 30 single-exon CDS, 8 two-exon, ZERO with 3+ exons. The tripartite leader is 5' UTR, not CDS. One intron is what the splice caller handles.
- 01:48Z min-len set to 20000, measured: 25000 would cut Siadenovirus 66 -> 54 and that band (23,811-27,209 nt) is the genus's own size range. 20000 costs 40 partial records family-wide.
- 01:48Z DUMP LAUNCHED: dump.sh -> dump.log, prefix A, ~5,015 genomes.
- 01:52Z DUMP COMPLETE in 4 min and check_dump.py PASS. 5,015 genomes, 155,187 feature rows (139,560 protein-coding, 15,627 blank-md5 non-coding), 21,947 distinct sequences. Index alignment clean.
- 01:52Z INTAKE RULE (SKILL §1) applied: dropped 139 'sp.' genomes -> intake/DROPPED_sp.tsv. NOTE the character differs from the parvo case: 118 are 'Human adenovirus sp.' HAdV-B3 clinical isolates, i.e. a submitter naming artefact on an abundantly-represented type, not an unclassified bin. 5 'Siadenovirus sp.' and 3 'Aviadenovirus sp.' ARE the DensovirinaeAmbi case. Dropping all of them costs 3% of Mastadenovirus and nothing in diversity.
- 01:52Z 67 genomes carry NO genus (bat/raccoon/porcupine/Adenovirus capreolus etc.) -> intake/NO_GENUS.tsv. Named viruses but unassignable to a module; held out of training data, kept as a routing-evaluation set.
- 01:52Z IN SCOPE: Mastadenovirus 4,244 | Aviadenovirus 447 | Siadenovirus 61 | Barthadenovirus 55 | Ichtadenovirus 2.
- 02:01Z TRIAGE PASS 1 (build_cols.py, cols1.log). Design: masses and ORF numbers do NOT bind by name -- both go to a pool resolved by homology (blastp >=80% id / >=85% qcov) then mmseqs clustering. Result: Mastadeno 3,356 genomes / 42 named features; Avia 430 / 24; Sia 59 / 22; Bartha 54 / 29. Homology rescued 27,237 Mastadenovirus features out of mass-named strings.
- 02:01Z PASS 1 DEFECTS to fix: (a) Aviadenovirus leaves 7,163 features / 1,810 seqs unplaced = 46% of the genus -- the ORF pool is not being resolved. (b) Mastadenovirus bound only 3,356 of 4,244 in-scope genomes; 888 bound nothing. (c) E1A_243R emptied by consolidation -- 289R/243R are nested splice forms and consolidation must not arbitrate between them, same exemption parvo needed for VP1/VP2. (d) UNC floor too loose: 152 UNC clusters in Mastadeno, most are 1-2 sequences.
- 02:10Z Wrote RESUME_HERE.md (checklist + the 5 fixed defects) and MATURE_PEPTIDES.md.
- 02:10Z MATURE PEPTIDES: measured, not assumed -- ZERO mat_peptide features in all four reference genomes and ZERO mature-peptide deposits among the 21,947 distinct sequences (every 'protein VIII' deposit is 227 aa, the precursor length). They must be built from published adenain cleavage sites, not binned. The one signal peptide in the family is E3 gp19K (Mastadenovirus); SignalP 6.0 is not installed and the published site is the authority anyway.
- 02:10Z Wrote probes.json: 4 modules, every PMID resolved against PubMed first. Probe pass already paid for itself -- it found E3 6.7K (PMID 15507617), which no annotation string in the dump would have surfaced.
