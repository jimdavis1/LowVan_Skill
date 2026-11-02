## 2 Nov 2026 — state made durable; read `reports/PARVOVIRIDAE_RESUME.md` first

**The restart doc now lives in the repo, not just here.** This file is the
chronological record and keeps growing; the *current state of play* — status
table, the four open items, the traps, and the restart checklist — is now
`LowVan_Skill/reports/PARVOVIRIDAE_RESUME.md`, tracked and pushed to GitHub.
Point a fresh session at that, then come back here for the history.

This file and `PROGRESS.md` are also copied into the repo as
`reports/parvoviridae-RESUME_HERE.md` and `reports/parvoviridae-PROGRESS.md`,
so a lost disk costs nothing. Re-copy them whenever either changes:

```bash
K=~/Library/CloudStorage/Box-Box/1_Projects/CEPI/Alsuviricetes/LowVan_Skill
cp ~/lowvan_scratch/parvo/RESUME_HERE.md $K/reports/parvoviridae-RESUME_HERE.md
cp ~/lowvan_scratch/parvo/PROGRESS.md    $K/reports/parvoviridae-PROGRESS.md
```

Three copies of the resume doc now exist and they are allowed to drift:
`~/lowvan_scratch/parvo/` is authoritative, `CEPI/Parvoviridae/` is the Box
mirror, `reports/` is the pushed snapshot. Per `resume-docs-live-in-scratch`,
**check scratch first.**

**Noticed while gathering the status table, not investigated:** two
Erythroparvovirus genomes — `10798.2563` and `1511900.137` — each make exactly
one call whose feature name is the empty string. All six declared features are
named, so this is in the caller or the GTO, not the JSON. Logged as open item 4.

Repo state at the time of writing: HEAD `0e111472`, equal to `origin/main`,
vouchers **416**, nothing running.

---

## 1 Nov 2026 — the Amdoparvovirus "third splice product" is NS3, and NS3 is called short

**The name question is answered: there is no third protein.** The 87/90/92 aa
product filed in `UNASSIGNED_TRACKING` is AMDV's own NS3. The literature names
only NS1, NS2 and NS3 for AMDV and gives exact sizes:

> "NS3 is a small protein with only **87 amino acids**, in which 60 amino acids
> from the N-terminus are identical to those of the NS1 and NS2. Hence, NS3 only
> has a unique sequence of 27 aa."
> — Huang, Luo, Cheng, Best, Bloom & Qiu, *Virology* 2014;452-453:23-31,
>   PMID **24606679** (already cited on both NS2 and NS3)

The same paper fixes NS2 at 114 aa: its GST constructs are "the unique
C-terminus of the AMDV NS2 (**aa 61-114**) and NS3 (**aa 61-87**)".

**So the module has NS3 wrong for AMDV.** Measured from the collections:

| length | n | taxa |
|---|---|---|
| 114 aa (NS2, correct) | 63 | 28314 |
| 87/90/92 aa (**NS3, correct**) | 33 | 28314 only |
| 68/71 aa | 2 (AMDV) | plus 50 from the other amdoparvoviruses |

- The 87/90/92 and the 68/71 forms are **the same reading frame** — both read
  `...TTSSDLR | LCDCNKQH...` after the shared 60-residue exon. The short ones
  stop early; the 87 aa ones continue `QHNQSNRWMCGNK...` to the real terminus.
- **No AMDV genome deposits both** (long ∧ short = 0 of 69). 33 deposit the long
  form, 2 the short. The short AMDV pair are truncated annotations.
- 28 of the 32 long deposits are **labelled NS2**, which is why they landed in
  the NS2 collection and were never available to build NS3. Textbook
  `deposited-name-is-a-claim`.

**Consequences in the shipped module (commit `3030543f`/`47f8800d`):**

1. `NS3` declares `min_len 60, max_len 80`. AMDV's real NS3 is 87-92 aa, so it
   is **outside the window and cannot be called** — and min/max *is* enforced on
   a splice feature, unlike a PSSM call.
2. The 11 AMDV references in `Splice-Variants/Amdoparvovirus/NS3.fasta` translate
   to **66 and 68 aa**, i.e. they encode the truncated acceptor. Verified by
   translating every reference:
   `28314 {66: 1, 68: 10}` vs `2169703 {68: 1, 71: 14}`, `1961072 {71: 5}`.
3. `NS2`'s 227 references are 195× taxon 28314 and all translate to 114 aa, so
   NS2 itself is fine — the mislabelled deposits did not leak into its refs.

**The fix, in order:**

1. Re-bin: move the 87/90/92 aa deposits out of `held_pssm_features/NS2.collection.fasta`
   into the NS3 collection. They are `28314.*|NS2_116..144` and `28314.*|NS3_1..4`.
2. Widen `NS3` to `min_len 59, max_len 95` in `features_all.json` (spans the
   66 aa of *Carnivore amdoparvovirus 2/3* through AMDV's 92).
3. Rebuild `Splice-Variants/Amdoparvovirus/NS3.fasta` with
   `build_refs_from_deposits.py` against the re-binned deposits — that script
   places the junction so the spliced product *is* the deposit, which is exactly
   what is needed when two acceptors compete in one frame. Do **not** use
   `build_splice_refs.py` here; a canonical donor/acceptor search is what picked
   the short acceptor in the first place.
4. Re-run coverage against `VA_parvo_route`, re-run pass 5, republish the four
   Amdoparvovirus pages, re-commit.

Expected effect: AMDV NS3 goes from uncallable to callable in ~33 of the panel's
AMDV genomes; the audit page's NS3 rate and the §07 limits paragraph both change.

Vouchers were **416** at the time of writing, i.e. baseline — the 235,000 peak
reclaimed on its own without a reboot. There is headroom for this re-run.

**Also done this session:** Dependoparvovirus AAP now carries three verified,
categorised D-lits for the non-canonical CTG initiator — commit `0e111472`.
See the commit message for which paper supports coordinates and which function.

---

# ============================================================
# 1 Nov 2026, end of session  --  NINE OF NINE SHIPPED
# ============================================================

## ALL NINE MODULES SHIPPED AND PUSHED
  Dependoparvovirus 899e2548   Erythroparvovirus   7753b358
  Bocaparvovirus    1ca7b62b   TetraCopiparvovirus 4b040030
  Hamaparvovirinae  dc502b41   DensovirinaeMono    a9422bf3
  Amdoparvovirus    3030543f + 47f8800d (splice refs)
  Protoparvovirus   4fdc5cdc
  DensovirinaeAmbi  64839ed0
  All four pages exist for all nine. 36 pages, all lint-clean, all indexed
  in ARTIFACTS.md, all private until shared from each page's Share menu.

## WHAT THE THREE DIRECTIVES RESOLVED TO
  1. coverage_cutoff: the 26 Oct caller measures coverage against the
     PROFILE, the 23 Sep one against the HSP's own query span (always ~1.0,
     never fired). Synced into VA_spina. gen_module_json.py also hardcoded
     0.65 and ignored a per-feature value; fixed (40a1f8fd). Used on
     DensovirinaeAmbi VP at 0.80.
  2. NS2/NS3 to the splice module: done for Amdoparvovirus and
     Protoparvovirus, literature-backed. DensovirinaeAmbi's NS features are
     separate collinear ORFs and stay PSSM-called -- measured.
  3. 2,000 genomes: still does not fit in one boot. make_module_repo.sh
     (ce43dc80) cuts routing to the references a panel actually wins on and
     took the leak from ~515 to ~120 vouchers per genome. Three modules
     (~770 genomes plus pass-5 re-annotation) took ipc.vouchers from 344
     after a reboot to 233,000, against a 250,000 stop.

## RESULTS
  Amdoparvovirus   same-coordinate calls 274/300 -> 0; clean 7 -> 262
  Protoparvovirus  NS1+NS3 210/299 -> 0; 285 of 297 clean (96.0%)
  DensovirinaeAmbi NS1+NS2 11 -> 0; VP1+VP2 46 -> 40; VP1_SPLICED called in
                   34 genomes; 148 of 171 clean (86.5%)

## THE ONE PIECE OF WORK LEFT, AND IT NEEDS A REBOOT
  DensovirinaeAmbi VP1+VP2 at identical coordinates in 40 genomes.
  Raising the coverage cutoff closed the partial-match route; VP.12 then
  took over at full coverage. VP.12 is built from TWO 578 aa deposits of
  taxon 2050976, whose VP deposits run 536-942 aa in clusters at ~577, ~645,
  ~759, ~823 -- a nested set filed under one name.
  RECIPE: bin taxon 2050976's VP deposits by size per SPECIES (P.species.tsv
  has the species for 9,409 of 9,432 genomes), the way vp_split.json did for
  the rest of the family; rebuild the VP profiles; re-run the 171-genome
  panel against VA_only_DensovirinaeAmbi (35 references, ~120 vouchers per
  genome, so ~21,000). REBOOT FIRST.

## ALSO OPEN
  Dependoparvovirus AAP. 34 of 34 deposits start with M; 1 of 169 calls
  does. AAP initiates at a non-canonical CTG in a third frame inside the
  capsid gene and upstream_ext is not finding a start codon. The "too long"
  flag on 77 of 176 genomes is a symptom; widening max_len would hide it.
  Amdoparvovirus's third splice product, 87/92 aa, deposited 28 times as NS2
  and 4 times as NS3, all from AMDV. No agreed name, left undeclared.
  Bocaparvovirus VP1+VP2 14/299 and Protoparvovirus 7/299 -- capsid tail.

# ============================================================
# 1 Nov 2026  --  CURATOR STEER RECEIVED; SPLICE REBUILD UNDER WAY
# ============================================================

CURATOR, 1 Nov (verbatim):
  1. The NS1 problem, where the pssm hits a partial orf, can be fixed by
     dialing up the coverage cutoff of the pssm match. We discussed this and
     fixed it last time.
  2. If NS2 and NS3 are spliced proteins, and this is backed by the
     literature, then they should go into the splice variant module, and
     should be called that way.
  3. This computer should be able to handle annotating 2,000 genomes. In the
     previous session we were working on doing the coverage using only
     closely related references.

## WHAT THAT RESOLVED
(1) The "fix from last time" is in $LOWVAN_KIT/annotate/annotate_by_viral_pssm.pl
    dated 26 Oct: cov = (q_to - q_from + 1) / qlen, i.e. against the PROFILE,
    where the 23 Sep copy divided by the HSP's own query span and so was
    always ~1.0 and never fired. SYNCED into VA_spina (= VA_parvo_route via
    symlink) on 1 Nov; backups at *.bak_23sep. get_splice_variant_features.pl
    synced too (21,919 -> 26,972 b).
(3) Already done: VA_parvo_route is VA_spina with Viral-Rep-Contigs cut to the
    176 Parvoviridae references. Validated 30 Oct against all 1,271. Machine
    was rebooted; ipc.vouchers was 344 at the start of this session.

## MEASURED, with the synced caller
  Amdoparvovirus   NS2 and NS3 are ONE-INTRON spliced CDSs, and the genus
                   uses THREE acceptors off ONE donor, not two:
                     -> 114 aa   = NS2 (110 of 145 deposits, 9 taxa)
                     ->  66-71 aa = NS3 (50 of 57 deposits); 18 of the 37
                        references reproduce their own genome's deposit
                     ->  87/92 aa, deposited as NS2 in 28 and as NS3 in 4,
                        ALL from AMDV taxon 28314. No agreed name, so not
                        declared. -> UNASSIGNED_TRACKING.
                   CAUTION: intron lengths quoted anywhere for these are
                   relative to each build's OWN slice start, so a 1352 in one
                   build and a 1352 in another are not the same genomic
                   intron. Compare products, not intron lengths, across
                   builds.
  Protoparvovirus  NS2 splices at a GC--AG donor, intron 1472 nt. Brute force
                   over every donor/acceptor pair in the 24 genomes with a
                   deposit reproduced the deposited protein byte-for-byte in
                   22. build_splice_refs.py will never find it: it picks a
                   canonical GT that sits in the same window, which is why
                   every canonical build came out 10-11 residues short.
  DensovirinaeAmbi NS1/NS2/NS3 are SEPARATE COLLINEAR ORFs, called at full
                   profile coverage at distinct coordinates. NOT this defect.
                   Its vp1/vp2 collision is its own problem -- see below.

## TOOL WRITTEN
  parvo/build_second_acceptor.py -- builds the second product of a SHARED
  DONOR when exon B is too short for tblastn (Amdo NS3's exon B is ~10 aa, so
  the NS3 query gave 7 references from 300 genomes while NS2 gave 227).
  Scores exon B ONLY: whole-protein identity is ~0.77 for any acceptor
  because exon A is shared, and scoring on it chose a 428 nt intron in 135 of
  173 genomes. --deposits accepts outright any candidate that equals that
  genome's own deposited protein; --intron constrains to a MEASURED intron.

## ALSO LEARNED
  build_splice_refs.py --nudge 12 (the default) is too tight for this family.
  Amdo NS2: 14 references at 12, 227 at 24. Proto NS2: 1 at 12, 162 at 24.
  Raise it before concluding a feature is not spliced.

## RESULT, Protoparvovirus (coverage re-run complete, 297 tables)
  NS1+NS3 same-coordinate calls:  210 of 299  ->  0
  Calls across 300 annotated genomes: NS1 294, VP1 293, VP2 268, SAT 214,
  NS2 210, NS3 27, NP1 3.
  Only residual collision: VP1+VP2 in 7 genomes (was 8) -- the known capsid
  tail, not the NS defect.
  NS2 references 187 (GC--AG donor, intron 1472, 22 deposit-exact);
  NS3 references 8 (intron 1839, every one deposit-exact).

## RESULT, Amdoparvovirus (coverage re-run complete, 300 of 300)
  Same-coordinate calls, every feature pair:   274 of 300  ->  0
  Genomes clean with copy_num applied:           7 of 300  ->  262 of 300
  Call rate on 289 routed: VP 99.3, VP2 99.3, NS1 94.1, NS2 84.4, NS3 15.9
  copy_num on VP and NS1; VP2 blocked (0% standalone, all 287 calls inside
  VP); NS2/NS3 never, they are spliced.
  qc_duplicate_calls: EXACT none.
  11 genomes have no contig at all -- protein-only deposits.

## STATE
  Amdoparvovirus   NS2 (227 refs) + NS3 (43 refs) installed as special:splice.
                   Five-genome test: NS2 115-124 aa, NS3 67-72 aa, two exons
                   each, no NS1 collision. 300-genome coverage re-run RUNNING.
  Protoparvovirus  NS2 references built (176, 22 deposit-exact, median 0.99
                   identity to the query). NS3 junction discovery running.
  DensovirinaeAmbi NOT STARTED. vp1/vp2 collide in 46 of 171; 36 of those are
                   taxon 185639, where VP1 IS SPLICED -- ORF-A + ORF-B joined
                   by "splice II", 816 codons, PMID 21775445. VP.1 (qlen 816)
                   hits q1-266 at 5285..4488 and q241-816 at 4384..2663, cov
                   0.706 on the 3' piece, and is extended onto VP2's ORF.
                   Raising coverage_cutoff alone does NOT fix it: VP.12 (qlen
                   578, a VP2-length deposit sitting in the VP collection)
                   then wins at cov 1.000 on the same span. Needs BOTH a
                   higher coverage_cutoff AND retirement of the VP2-length
                   profiles from VP.

## VERIFIED PMIDs (eutils, 1 Nov)
  16378968  AMDV transcription profile, single promoter, J Virol 2006
  17428872  AMDV R2 mRNA is tricistronic (NS2, VP1, VP2), J Virol 2007
  24606679  AMDV small nonstructural proteins during infection, Virology 2014
  2942705   MVM splice junction map, J Virol 1986
  8212598   PPV NADL-2 transcription and translation map, Virology 1993
  21775445  Acheta domesticus densovirus expression strategy, J Virol 2011

# ============================================================
# 31 Oct 2026, 23:10  --  6 OF 9 SHIPPED; 3 HELD ON A CALLER DEFECT
# ============================================================

## SHIPPED + PUSHED (kit modules/, Box work/, runtime, pages live)
  Dependoparvovirus    899e2548    Erythroparvovirus   7753b358
  Bocaparvovirus       1ca7b62b    TetraCopiparvovirus 4b040030
  Hamaparvovirinae     dc502b41    DensovirinaeMono    a9422bf3
  Two pages each (registry, vocabulary saturation) -- 12 published, all
  lint-clean, all in ARTIFACTS.md.

## HELD, NOT SHIPPED: Amdoparvovirus, Protoparvovirus, DensovirinaeAmbi

MEASURED, same-coordinate calls per module (pass 4 tables):
  Amdoparvovirus     NS1+NS2+NS3  274/300 ; NS2+NS3 23 ; NS2 called 2x 288/300
  Protoparvovirus    NS1+NS3      210/299
  DensovirinaeAmbi   VP1+VP2       46/171 ; NS1+NS2 11 ; VP1 called 2x 47
With copy_num applied, Amdoparvovirus drops to 7 of 300 genomes clean (2.3%).

THE COLLECTIONS ARE CLEAN -- this is NOT the VP1/VP2 mislabelling again:
  Amdo  NS1 408-654 aa (n=167)  NS2  80-116 aa (n=145)  NS3 59-92 aa (n=57)
  Proto NS1 467-885 aa (n=802)  NS2 103-197 aa (n=247)  NS3 79-107 aa (n=26)

ROOT CAUSE, measured: **max_len is never enforced on a PSSM-called feature.**
In annotate_by_viral_pssm.pl, max_len is read only into $non_pssm_feat (line
~555) and used only at ~747, i.e. for special/leftover features. min_len is
used only as a read-through plausibility test (~424). So:
  Amdo  NS2 declares min 75 / max 122  and is called at 600 aa
  Proto NS3 declares 79-107 aa         and is called at 669 aa in 176/299
NS2 and NS3 are SPLICED: they share the N-terminal exon with NS1, so a short
profile necessarily hits inside the NS1 ORF, and upstream_ext/downstream_ext
then extend the call to NS1's own start and stop. The translation printed for
NS2 is byte-for-byte NS1's.

## SECOND FINDING: the whole coverage run used a STALE caller
  ~/lowvan_scratch/VA_spina/annotate_by_viral_pssm.pl   23 Sep  41,708 b
  $LOWVAN_KIT/annotate/annotate_by_viral_pssm.pl        26 Oct  45,365 b
95 differing lines. The kit version (a) computes coverage_cutoff against the
PROFILE length rather than the HSP's own query span, and (b) skips an
unreadable reference contig instead of dying on it. install_module.py does
not sync the perl, and run_gto_eval calls the repo copy by path, so the kit
fix never reached any parvo run.

TESTED on Amdoparvovirus 1093101.3, both callers, same repo and PSSMs:
  23 Sep: ns1 110-1909, ns2 110-1909, ns3 110-1909, ns2 1924-2151, vp1, vp2
  26 Oct: ns1 110-1909, ns2 110-1909, ns3 110-1909,                vp1, vp2
The newer caller removes the DUPLICATE NS2 -- but the one it removes is the
CORRECT 76 aa call at 1924-2151, and it keeps the 600 aa one. Neither caller
fixes the collision. Do not assume syncing the perl is the fix.

## CURATOR DECISIONS NEEDED
  1. Enforce max_len on PSSM calls (reject, or crop, or do not extend past
     it)? This changes behaviour for all 97 modules, so it is not mine.
  2. Should Amdo/Proto/DensoAmbi NS2 and NS3 be `special: splice_variant`
     with nucleotide references, like Astroviridae ORF1AB? They are not
     collinear, which is why no PSSM can place them.
  3. Sync the kit caller into the test repo and re-run coverage? That is
     ~2,200 genomes and will cross the voucher ceiling, so it needs a reboot
     and the watchdog.

# ============================================================
# 31 Oct 2026, 22:05  --  CONTINUOUS-SHIP MODE (curator instruction)
# ============================================================
Curator, 31 Oct: log progress consistently, do NOT auto-restart the session on
reboot, and push every artifact / GitHub update AS SOON AS IT COMPLETES so a
crash loses nothing.

## How state is kept now
  PROGRESS.md (here AND Box CEPI/Parvoviridae/PROGRESS.md) -- one timestamped
  line per milestone, written by progress.sh. READ ITS LAST LINES FIRST.
  state/            per-stage markers; resume_driver.sh skips what is done
  logs/resume_driver.log, logs/vouchers.tsv (the voucher curve)
  NO LaunchAgent is installed. After a reboot, re-run by hand:
      bash ~/lowvan_scratch/parvo/resume_driver.sh

## Per-module scripts (run in this order, per module)
  resume_driver.sh      -> pass5_pages.sh <M>   (copy_num, install, rescore,
                           dup-call QC, registry, rarefaction, coverage, pages)
  ship_parvo.sh <M>     -> kit modules/, Box work/, runtime; re-renders registry
  commit_parvo.sh <M> <msgfile> ["Title<TAB>artifactid" ...]
                        -> published.tsv, reports/, ARTIFACTS.md, commit, push

## Pushed so far
  c28462af  capsid_length_audit.py, SKILL.md step 2, dotted-key parser fix
  899e2548  Dependoparvovirus shipped; registry + saturation pages published;
            published.tsv backfilled with 17 rows (Astro + Parvo) that existed
            only in ARTIFACTS.md and would have been lost on regeneration
  b4314489  gen_facts.py --panel-note/--panel-label

## Defects found and fixed tonight
  pass5_pages.sh gave apply_copy_num.py features_all.json; it is keyed by
  MODULE so it needs <M>_Viral_PSSM.json. It also regenerated the module JSON
  AFTER copy_num, which strips it (gen_module_json builds from
  features_all.json). Both fixed; a quality RE-SCORE and qc_duplicate_calls.py
  were added, because the .qual.gto files were scored before copy_num existed.

## NOT DONE
  1. The nine coverage-audit pages. artifacts/generators/<slug>_prose.json does
     not exist for any module and pass5 SKIPS the page without it. The numbers
     are ready (make_g_files.py writes the per-module G.* dump files gen_facts
     needs; P.species.tsv was fetched from the BV-BRC genome API with curl).
     A first attempt at the Dependoparvovirus prose was stopped by a safety
     classifier mid-write and the partial file was deleted. NOT a data problem.
  2. Dependoparvovirus AAP: called 221-263 aa in AAV7/AAV8/avian/goose genomes
     against a 125-216 window learned from AAV2 (~209 aa). 75 of 176 genomes
     carry "Feature is too long". NOT investigated -- the window may be too
     narrow, or the long calls may be wrong. This is why only 87 of 176 are
     genuinely clean while 158 carry no genome/contig flag.
  3. Two Dependoparvovirus exemplars route to Hamaparvovirinae and
     Protoparvovirus (divergent, no close rep contig).
  4. Residual VP/VP2 same-coordinate calls: DensovirinaeAmbi 46/171,
     Bocaparvovirus 14/299, Protoparvovirus 8/299. Not investigated.

# ============================================================
# 31 Oct 2026, 21:35  --  RESUMED AFTER THIRD KERNEL PANIC (12:15)
# ============================================================

## What happened
Panic 12:15:46 today (report written 21:19 on restart), identical to 28 + 30
Oct: "Cannot grow ipc space beyond IVAC_ENTRIES_MAX. Some process is leaking
vouchers". Panicked task: bash. Uptime at panic 20.5 h (booted 30 Oct 15:48
from the previous panic -- so no clean reboot preceded this run).
The 30 Oct mitigation (VA_parvo_route: routing cut to the 176 Parvoviridae
references) reduced launches but did not stop the leak.

## The leak, MEASURED 31 Oct 21:2x (zprint is readable without root)
  100 launches each, delta in ipc.vouchers in-use:
    blastn +65, tblastn +50, makeblastdb +79, BV-BRC perl +91, conda python +80
    /bin/sh +0, /usr/bin/true x2000 +0
  5 real Hama genomes (run_gto_eval --jobs 1, 149 s): ipc.vouchers, bank_task,
  bank_account, Bank.ledger each +~2,580 => ~515 per genome, still held after.
  Non-Apple binaries leak; Apple platform binaries do not. CrowdStrike Falcon
  is in all three panic process lists (suspect, not proven).
  Gauge: zprint | awk '/^ipc.vouchers /{print $7}'. Curve: logs/vouchers.tsv.

## State at the panic (pass 4 = pass4.sh steps 0-2 + pass4_cov.sh)
  Coverage pass 4 COMPLETE for 8 modules (eval9 logs, 09:23-12:00).
  Hamaparvovirinae 162/300 -- all 162 outputs validated as parseable JSON.
  pass5_pages.sh had NEVER been run (no copynum/json5 logs).
  Generations: cov_precapsidfix/ (orig), cov_pass1/, cov_pass3/, cov/ = pass 4.

## MEASURED 21:45: VP and VP2 called at identical coordinates (genomes)
  module               precapsidfix  pass1    pass3    pass4(cov/)
  Dependoparvovirus    1/176         1/176    1/176    1/176
  Protoparvovirus      8/299         9/299    8/299    8/299
  Bocaparvovirus       16/300        14/299   14/299   14/299
  Erythroparvovirus    281/300       1/300    269/300  1/300
  Amdoparvovirus       0/300         0/300    0/300    0/300
  TetraCopiparvovirus  117/300       155/300  0/300    0/300
  DensovirinaeAmbi     49/171        52/171   46/171   46/171   <- largest residual
  DensovirinaeMono     -             16/40    0/40     0/40
  Hamaparvovirinae     -             0/297    0/297    0 so far
  Residuals (Boca 14, DensoAmbi 46, Proto 8) NOT investigated yet.

## RUNNING NOW: resume_driver.sh (nohup caffeinate, started 21:31)
  A  Hama coverage resume  ->  B  pass5_pages.sh per module (markers in state/)
  Watchdog: no new stage above 250k vouchers, kill + NEEDS_REBOOT above 350k.
  Log logs/resume_driver.log; state/STATUS; re-run is idempotent.
  LaunchAgent for auto-resume after reboot is WRITTEN, NOT INSTALLED (needs
  the curator): com.jdavis.lowvan.parvo.plist in this dir.

## NOT DONE (not started, not blocked)
  1. VP/VP2-at-identical-coordinates re-measure, all generations
  2. coverage-audit prose: artifacts/generators/<m>_prose.json do not exist,
     so pass5 SKIPS the 9 coverage-audit pages; write them after the data
  3. regenerate registry pages after copy_num; 27 pages to publish
     (only the 9 string-collapse pages are live + in ARTIFACTS.md)
  4. modules/<Mod>/ into the kit, ARTIFACTS.md, commit; kit has UNCOMMITTED
     SKILL.md, apply_copy_num.py, evaluate_coverage.py, capsid_length_audit.py

# ============================================================
# 31 Oct 2026, 03:10  --  PASS 3 LAUNCHED (bv82x8mwl)
# ============================================================

One more defect, found because TetraCopi's 394-485 aa VP2 cluster SURVIVED
the rebuild that was meant to remove it.

The N-anchored test compared each member against its species' LONGEST parent
deposit only. Porcine parvovirus 4's longest VP is a 914 aa outlier while the
species sits at 728, and its three fragments start at the 728 aa form, so they
tested clean. Against EVERY parent member of the species it finds nine more:

  Erythroparvovirus  B19 "VP2" 617, 608, 605 aa -> VP1 residue 1 (781 aa)
  TetraCopi          PPV4 485, 427, 394 aa; human parvovirus 4 596 aa
  Bocaparvovirus     porcine bocavirus 3 516 aa; carnivore bocaparvovirus 2 511

The three B19 ones are the 602-617 aa VP2 cluster that was looked at twice
and accepted, because 602 is near enough B19's real VP2 at 554 that no length
rule flags it. Truncated VP1.

Final moves: Dependo 2+1 | Proto 4+3+1 | Boca 40+26+5 | Erythro 91+6+3 |
Amdo 4+2 | TetraCopi 106+96+13 | DensoAmbi 12+1 | DensoMono 2+2+1 | Hama 0.

capsid_length_audit.py: all nine clean. TetraCopi VP2 is now three clusters
at 551-555 aa (was six spanning 394-713).

Both build_cols.py and capsid_length_audit.py carry the all-parents test with
the PPV4 case in the comment.

## RUNNING: pass3.sh, task bv82x8mwl
finalize -> install -> move cov/{anno,tbl,report.tsv} to cov_pass1/ -> coverage
on VA_parvo_route, 9 modules incl Hamaparvovirinae. ~4-5 h.

Coverage generations: cov_precapsidfix/ = original, cov_pass1/ = 30 Oct
length-class pass, cov/ = this one.

WHEN IT FINISHES: re-measure VP-and-VP2-at-identical-coordinates per module
(the awk is in the 03:5x block below) and compare all three generations.
Then: copy_num on fresh data (dotted-key parser fix), vocabulary rows, the 27
artifact pages, commit the kit.

UNCOMMITTED in the kit: SKILL.md (step 2, now synced from the installed copy
which was the only one carrying the 30 Oct work), capsid_length_audit.py
(new), apply_copy_num.py, evaluate_coverage.py.

# ============================================================
# 31 Oct 2026, 02:45  --  ADDENDUM: the second test (N-anchored)
# ============================================================

The anchor rule alone still left TetraCopiparvovirus with two bad VP2
clusters, 606-713 aa (porcine parvovirus 2) and 394-485 aa (porcine
parvovirus 4). Neither species has a genome depositing both forms, so no
anchor reaches them.

MEASURED: all six of those "VP2" deposits start at VP1 RESIDUE 1 of their own
species' VP1. They are truncated VP1, not an alternative start. A nested form
shares the parent's C-terminus; these share its N-terminus. Length cannot tell
the two apart -- both are just "shorter" -- which is why two passes missed it.
Those six profiles called VP2 at VP1's own coordinates in 132 genomes.

build_cols.py now sends them back to the parent form, where fasta-cluster-pssm's
existing substring-drop discards them. Moves: Proto 1, Boca 3, TetraCopi 4,
DensoAmbi 10, DensoMono 1.

capsid_length_audit.py: ALL NINE MODULES CLEAN (it gained the same test, and
now compares a form only against the forms it actually nests with -- connected
components, because Boca nests VP/VP2 and also NS1/NP1, and a flat set made
NS1 and VP1, both ~685 aa and unrelated, report each other).

SKILL.md: the step-2 section was only ever written to the INSTALLED copy
(~/.claude/skills/lowvan-module/SKILL.md). The kit copy was still at 28 Oct.
Both now carry the anchor rule, the N-anchored rule and the rejected 0.85
ratio, and are byte-identical (80,542 bytes).

RUN NEXT, one launch:   bash pass3.sh
  finalize -> install -> park cov/ into cov_pass1/ -> coverage on VA_parvo_route
Then re-measure VP-vs-VP2-at-identical-coordinates per module and compare
against the table in the 03:5x block below.

# ============================================================
# 31 Oct 2026, 03:5x  --  ANCHOR RE-BIN: second capsid pass
# ============================================================

## What the 30 Oct coverage re-run actually showed (cov/ vs cov_precapsidfix/)

Measured: genomes where VP and VP2 are called at the SAME start/stop.

| module              | before | after (30 Oct fix) |
|---------------------|--------|--------------------|
| Erythroparvovirus   | 281/300| **1/300**          |
| Bocaparvovirus      |  16/300|   14/299           |
| Dependoparvovirus   |   1/176|    1/176           |
| Amdoparvovirus      |   0/300|    0/300           |
| Protoparvovirus     |   8/299|    9/299  REGRESSED |
| TetraCopiparvovirus | 117/300|  155/300  REGRESSED |
| DensovirinaeAmbi    |  49/171|   52/171  REGRESSED |

So the 30 Oct fix was right where it was aimed (Erythro) and wrong elsewhere.

## Why it regressed, measured not guessed

TetraCopi's 155 collisions came from exactly three VP2 clusters:
  VP2.5 (n=3, 606-713 aa, porcine parvovirus 2)   -> 92 collisions at VP1's 1032 aa
  VP2.7 (n=3, 394-485 aa, porcine parvovirus 4)   -> 40 collisions at VP1's  728 aa
  VP2.4 (n=5, 804-880 aa, porcine parvovirus 3)   -> 23 collisions at VP1's 931 aa

PPV3 is the real mislabelling: genome 1126384.191 deposits VP1=883 and
VP2=513 and the 513 is an exact C-terminal substring of the 883. Six further
PPV3 genomes deposit "VP2" at 688-883 aa, i.e. VP1 length.

The 30 Oct rule MISSED it: groups() split PPV3's lengths into 4 classes
(925/883/880, 806/804, 688, 513) and the guard `2 <= len(gs) <= len(CAPSID)`
then skipped the whole species.

Proto and TetraCopi also LOST VP2 calls: 20 porcine parvovirus 5 and 22
porcine parvovirus 6 genomes. Those two species deposit NO VP2 anywhere in
the collections (every genome deposits one capsid, at 991 and 1189 aa), so
the VP2 they used to get came from a contaminated cross-species profile.
Losing it is correct.

## The rule that replaced it: MEASURED ANCHORS, not a ratio

Tried and REJECTED: "a VP2 member at >= 0.85 x its species' VP1 median is a
VP1." Simulated first -- it moves 143 genuine Amdoparvovirus VP2 deposits,
because the VP1->VP2 offset is family-specific: 227 aa in Erythroparvovirus,
~370 in Tetraparvovirus, ~42 in Amdoparvovirus (689 -> 647). No single
fraction separates the forms. DO NOT re-propose this.

What is in build_cols.py now: where one genome deposits two forms, the
shorter an exact C-terminal substring of the longer, AND the existing labels
agree with that order, that pair MEASURES both forms for that species. Every
member of that species then goes to the nearest anchor in RELATIVE distance.
Anchors outrank the length-class ranking, which stays as the fallback for
species with a nesting pair but no label-consistent one.

Moves: Dependo VP2->VP3 2, VP2->VP 1 | Proto 4+3 | Boca 40+26 | Erythro 91+6
| Amdo 4+2 | TetraCopi 106+96 | DensoAmbi 1 | DensoMono 2+2 | Hama 0.

Collections after (VP vs VP2 ranges): TetraCopi 602-1190 vs 250-713 (the
804-883 contaminants are gone); Amdo 642-702 vs 582-647 (42 aa offset kept);
Boca 417-724 vs 339-675; Erythro 670-947 vs 375-617.

capsid_length_audit.py: **all nine modules clean.**

## capsid_length_audit.py also had two bugs, both fixed

1. Without --forms, test 2 compared EVERY feature, so Bocaparvovirus NS1 and
   VP1 -- both ~685 aa in porcine bocavirus 3 and unrelated -- reported each
   other. It now audits only forms that actually nest.
2. A flat "forms that nest" set was still wrong: Boca nests VP/VP2 and also
   NS1/NP1 (2 deposits), which put NS1 back in the same set as VP. Test 2 now
   runs per CONNECTED COMPONENT of the nesting graph.

## WHERE THIS STOPPED / what is next

Collections are REWRITTEN (build_cols.py full run, logs/cols_anchor.log).
The installed PSSMs are therefore STALE for all 8 changed modules.

  1. bash rebuild_capsids.sh <8 modules incl DensovirinaeMono>   (sequential)
  2. bash finalize_capsids.sh <same 8>
  3. install, then a THIRD coverage run, then re-measure the table at the top
  4. then: copy_num on fresh data (dotted-key parser fix), vocabulary rows,
     the 27 artifact pages, commit the kit

build_cols.py needs mmseqs and blast on PATH:
  export PATH=/Users/jdavis/anaconda3/envs/lowvan/bin:/Applications/BV-BRC.app/runtime/bin:$PATH

cov/ now holds the 30 Oct (second) run; cov_precapsidfix/ the original.
Move cov/ aside again before the third run -- run_gto_eval resumes any genome
whose .ann.gto and table are both non-empty and would skip everything.

# Parvoviridae — state as of 2026-10-30, 22:05 EDT

## RUNNING NOW: cov_run7.sh — coverage RE-RUN on the 7 capsid-fixed modules

Log `logs/cov_run7.log`, per-module `logs/<Mod>.eval7.log`. 1,860 genomes at
~7.3/min with the routing wrapper, so expect ~02:20. Kill order if needed:
cov_run7.sh, then run_gto_eval.py, then `pkill -f annotate_by_viral_pssm.pl`.

Hamaparvovirinae (297/300) and DensovirinaeMono (40/40) moved NO capsid
sequence, so their cov_run6.sh results still stand and they are not re-run.

**Pre-fix coverage is held in `cov_precapsidfix/<Mod>/`** (anno, tbl,
report.tsv). It had to be moved aside: run_gto_eval RESUMES any genome whose
.ann.gto and table are both non-empty, so leaving it in place would have
skipped every genome. The GTOs and contigs were NOT touched, so nothing is
re-staged from BV-BRC.

### The measurement to make when it finishes
Per module, count genomes where VP and VP2 are called at the SAME
start/stop. That is the defect, and the before-figures are:

    Erythroparvovirus   281 of 295      TetraCopiparvovirus 117 of 269
    Bocaparvovirus       16 of 296      Protoparvovirus       7 of 293
    Dependoparvovirus     0 of 173      Amdoparvovirus        0 of 287

## DONE SINCE 21:30 — collections fixed, PSSMs rebuilt, 7 modules installed

1. `build_cols.py` re-bins nested capsid forms on the literature's size, per
   SPECIES, learned from same-genome nesting. Collections written
   (`logs/cols18_write.log`); `--no-write` reproduces them exactly
   (`logs/cols20_verify.log`).
2. `rebuild_capsids.sh` rebuilt VP/VP2/VP3/VP4 alignments + PSSMs + leftovers
   for the 7 modules that moved. Erythroparvovirus VP2 is now 554/558 aa
   profiles instead of one 781 aa VP1 profile.
3. `finalize_capsids.sh` re-derived cutoffs and rewrote each module JSON.
4. All 7 installed into VA_spina after `--check` returned 0 problems each.
   **Exactly 10 features changed in the runtime, all capsid, no feature set
   altered, 97 modules / 1132 features intact.** Backup:
   `VA_spina/Viral_PSSM.json.pre_capsidfix_20261030`.

        Bocaparvovirus      VP2  bit 110->309, max_len 734->715
        Protoparvovirus     VP2  bit  32->88,  max_len 779->628
        TetraCopiparvovirus VP2  bit 127->117, max_len 988->935
        Erythroparvovirus   VP2  max_len 827->601 | VP min_len 520->565
        Dependoparvovirus   VP2  max_len 768->634
        TetraCopiparvovirus VP   min_len 521->560
        Amdoparvovirus      VP 1033->1027, VP2 1061->1067

   The VP2 cutoffs ROSE because the VP1-length contaminants had been dragging
   each weakest-genuine-member down; the max_len drops are the same thing.

## KNOWN AND DELIBERATE, NOT FIXED

- **7 residual singletons** the audit still reports (Amdo 3, TetraCopi 1,
  DensovirinaeAmbi 2, DensoMono 1). Each lacks same-genome nesting evidence,
  and the generator does not guess without it. The two DensovirinaeAmbi ones
  are the Diciambidensovirus 588/573 aa proteins, which are NOT truncations
  of its real 202-205 aa capsid -- expected, leave them.
- **5 Bocaparvovirus VP2 members at 669-675 aa** (Aveparvovirus). For
  *Aveparvovirus galliform1* the pairs say VP1=675-677 and VP2=535-536, so
  these look like VP1 -- but each is `<Genus> sp. <isolate>` or a one-off
  isolate name, so the 3-token species key makes it a species of one with no
  pair. A depth-2 species fallback was written, measured, and REVERTED: it
  moved 1-2 members either way, did not resolve these, and merging by default
  would have put Porcine parvovirus 1 (VP1 729) with Porcine parvovirus 3
  (VP1 925). The re-run's VP2-at-VP1-coords count will say whether these 5
  actually cost anything.
- **Protoparvovirus VP bit_cutoff is 37** and was 38 before this work. Very
  low, pre-existing, NOT caused by the capsid fix. Worth a look.
- macOS bash 3.2 has no `wait -n`. The throttle in `rebuild_all.sh` and the
  first `rebuild_capsids.sh` spun instead of blocking -- a core burned and
  700 kB of errors in under a minute, and every module launched at once
  anyway. `rebuild_capsids.sh` is now sequential; **rebuild_all.sh still has
  the bug.**

## AFTER THE RE-RUN, IN ORDER
1. Measure VP2-at-VP1-coords per module (above). 2. copy_num on the fresh
data, with the dotted-key parser fix. 3. vocabulary rows. 4. 27 pages.
5. Commit the kit: capsid_length_audit.py is new and SKILL.md step 2,
apply_copy_num.py and evaluate_coverage.py are all modified, UNCOMMITTED.

---

# Parvoviridae — state as of 2026-10-30, 21:30 EDT

## THE CAPSID COLLECTIONS WERE CONTAMINATED. COLLECTIONS ARE FIXED;
## PSSMs ARE NOT YET REBUILT.

The deposited VP1/VP2/VP3 names are submitter claims and are frequently wrong.
Because the forms are ~100% identical where they overlap, a deposit under the
wrong name cannot be separated again by homology, self-recall or call rate.

Literature, verified this session: B19 VP1 781 aa, VP2 554 aa -- VP2 is VP1
without the 227-residue VP1u (PMID 11799199, Dorsch 2002 J Virol, resolved via
PMC135890 and checked against the PubMed title; an earlier GUESS of 11907233
was an African swine fever paper). RefSeq NC_000883 agrees (VP1 2624..4969,
VP2 3305..4969) and is corroboration, not the authority. Tetraparvovirus VP1
is 900-950 aa; PARV4 ORF2 = 914 aa.

What it had done:
  * Erythroparvovirus: 31 unique deposits named VP2 were 781 aa = VP1.
    fasta-cluster-pssm drops a cluster whose consensus is a SUBSTRING of
    another, and the real 554 aa VP2 is a substring of VP1 -- so both genuine
    VP2 clusters were dropped ("1 within 2, 3 within 2") and VP2 SHIPPED AS A
    VP1 PROFILE. VP2 was called at VP1's exact coordinates in 281 of 295
    genomes. Self-recall, call rate and the anno string all looked right.
  * TetraCopiparvovirus: same defect both ways. 22 genomes deposit
    "VP1"=555 aa and "VP2"=925 aa, inverted. 117 of 269 genomes affected.
  * Smaller counts in Bocaparvovirus (16 genomes) and Protoparvovirus (7).

The test is STRUCTURAL: VP2 is an exact C-terminal substring of VP1 from the
SAME GENOME. 34 of 34 TetraCopiparvovirus pairs nested exactly.

### Fixed in build_cols.py -- `nested capsid relabelled on length class`
Per SPECIES (never per genus: Tetraparvovirus alone holds VP1 at 914, 925,
931, 933 and 1032 aa). Lengths are grouped with a 5% gap -- 727/714/713/690
is ONE class, and ranking raw lengths moved 317 real Protoparvovirus VP1s.
Ranking applies ONLY to species that show the nesting in one of their own
genomes; without that guard it demoted the 202-205 aa Diciambidensovirus
capsid, which is real and complete. Species with no pair fall back to the
majority holder of their length class. Moves, all logged:

    TetraCopiparvovirus VP2->VP 105, VP->VP2 67 | Erythroparvovirus VP2->VP 91
    Bocaparvovirus VP2->VP 40 | Dependoparvovirus VP2->VP3 7 | Proto 5
    Amdoparvovirus 2 | DensovirinaeAmbi 1 | Hamaparvovirinae 0 (no VP2)

Every VP2 upper bound now sits below its module's VP median. `cols18_write.log`.

### New kit script: `capsid_length_audit.py` (and SKILL.md step 2 rewritten)
Independent check, three tests: same-genome nesting, minority label on a
length class, nested form at the parent's length. Clean on Dependo, Proto,
Boca, Erythro and Hama. **7 residual singletons** in Amdo (3), TetraCopi (1),
DensovirinaeAmbi (2, both the Diciambidensovirus 588/573 aa proteins that are
NOT truncations of its short capsid -- expected, leave them) and DensoMono
(1). The generator deliberately does not guess where there is no nesting
evidence; the audit reports them for the curator.

## WHAT MUST HAPPEN NEXT, IN ORDER

1. **Wait for cov_run6.sh.** DensovirinaeAmbi + DensovirinaeMono still to go.
   DO NOT install anything while it runs -- VA_parvo_route symlinks
   VA_spina's Viral-PSSMs and Viral_PSSM.json, so an install changes the repo
   under the running evaluation. That voided an earlier run on 29 Oct.
2. **Rebuild VP/VP2/VP3/VP4 alignments and PSSMs** for the 8 modules that
   moved. A contaminated collection does not become clean by renaming rows:
   the profiles were built from the mixed set. Then leftovers, cutoffs,
   bounds, module JSON, install (SERIALLY).
3. **Re-run coverage for every module whose capsid profiles changed.** The
   VP2 call rates measured before this fix are void. ~7.3 genomes/min with
   the routing wrapper, so all nine is about 5 h, not 20.
4. copy_num, vocabulary rows, 27 pages, commit.

## THE ROUTING WRAPPER (the reason step 3 is affordable)

`~/lowvan_scratch/VA_parvo_route` = every VA_spina entry symlinked EXCEPT
Viral-Rep-Contigs, which holds only the 176 Parvoviridae references. No perl
and no Viral_Annotation file was changed -- curator, 30 Oct.
Routing ran one blastn per reference file over 1,271 files: ~2,600 process
launches per genome, ~95% of runtime, and the trigger for both kernel panics.
Validated before use: a one-call screen of all 2,200 sampled genomes against
all 1,271 references reproduced the recorded routing on 1,982/1,982 finished
genomes; of the 524 left, 520 win on a Parvoviridae reference and 4 score
<150 everywhere. One finished genome re-run gave a byte-identical table.
Screen: `logs/route_screen.tsv`. Rate went 2 -> 7.3 genomes/min.

## KIT FIX ALREADY MADE: dotted feature keys

`apply_copy_num.py` and `evaluate_coverage.py` split family_assignments on "."
and took field 2, reading "Erythroparvovirus.7.5K.1" as feature "7". 7.5K and
9.5K reported 0.0% on 298 genomes that call them 283 and 284 times. Now strips
only a trailing numeric profile id (handles leftovers like "NP1.3_1"). Across
all 37,683 assignments in every scratch evaluation the old and new parsers
differ ONLY on those two features. Erythroparvovirus now: 9.5K 95.3%,
7.5K 95.0%, both earning copy_num.

## copy_num, MEASURED BUT NOT YET APPLIED -- and VP2 figures are now VOID

Re-measure after step 3. Pre-fix numbers, for reference only:
Dependo VP 100/AAP 97.7/NS1 94.3, Rep forms 50-59 (nested, expected low);
Proto NS1 99.3/VP 99.3; Boca NS1 100/NP1 99.7/VP 99.0; Erythro NS1 99.0/
VP 99.0/9.5K 95.3/7.5K 95.0/11K 63.4; Amdo NS1 100/NS2 100/NS3 99.0/VP 95.7;
TetraCopi NS1 99.3/VP 98.9/ORF3 14.4.
Denominator decision stands: genomes carrying >=1 call from the module.

---

# Parvoviridae — state as of 2026-10-30, 20:30 EDT (written after a kernel panic)

## THE MACHINE PANICKED AT 15:47:59 — nothing is running

Panic string, both times (28 Oct 05:28 AND 30 Oct 15:47):
`Cannot grow ipc space beyond IVAC_ENTRIES_MAX. Some process is leaking vouchers`
Victims were taskgated (28 Oct) and launchd (30 Oct) -- both on the exec path.
Both panics hit mid-coverage-run with perl/blastn workers live, after 5.0 and
2.5 days of uptime. Reports: /Library/Logs/DiagnosticReports/panic-full-*.panic

CAUSE (strong circumstantial, not proven): process-launch volume, not load.
Routing in annotate_by_viral_pssm.pl runs ONE blastn PER REFERENCE CONTIG
through a shell pipe -- 1,271 files in VA_spina/Viral-Rep-Contigs, so ~2,600
process launches per genome before a single PSSM is searched. At ~2
genomes/min that is ~90 execs/s sustained for days. trustd/taskgated/Falcon
CPU in the panic process table fit an exec storm. Cutting --jobs from 6 to 4
after the first panic did not help: the trigger is cumulative launches since
boot, not concurrency. The leak itself is an OS/security-agent bug (a user
process should not be able to panic the kernel); the coverage run triggers it.
Also: Activity Monitor was at 9.1 GB RSS / 53,000 CPU-s at panic time. Close it.

ROUTING IS ALSO THE RUNTIME. Measured on Hamaparvovirinae 1820046.201:
per-file loop 141 s single-core; ONE blastn with all 1,271 references
concatenated as queries: 1 s, same winner (Hamaparvovirinae.1.dna, bit 4985).
get_blastn_bit = max single-HSP bitscore per file; ties -> last file in sort
order (`unless ($match_bit < $best)`). A one-call version can reproduce that
exactly. NOT applied -- changing the runtime caller is the curator's call.

## COVERAGE STATE (cov_run5.sh died with the machine)

    Dependoparvovirus   176/176  done (124 rAAV held out)
    Protoparvovirus     300/300  done (299 tbl; 10788.3854 is a genuine no-call)
    Bocaparvovirus      300/300  done
    Erythroparvovirus   300/300  done (32 lost GTOs restored by prestage_gto.sh)
    Amdoparvovirus      300/300  done
    TetraCopiparvovirus 300/300  done
    Hamaparvovirinae      6/300  <- interrupted; run_gto_eval RESUMES
    DensovirinaeAmbi      0/184  GTOs staged
    DensovirinaeMono      0/40   GTOs staged

All GTOs are pre-staged. Restart = edit cov_run5.sh's module list down to the
last three and rerun (it resumes per genome). ~518 genomes, ~4.3 h at --jobs 4
with the per-file routing.

## SKILL / KIT SYNC (checked 30 Oct)

- Kit (Box LowVan_Skill) == origin/main == GitHub, 15e5d5ec, clean tree.
- ~/.claude/skills/lowvan-module is a 13 Oct copy -- 51 commits behind on
  skill/; SKILL.md 67 KB vs 75 KB. Sessions have been loading the OLD skill.
- VA_spina/annotate_by_viral_pssm.pl (the copy run_gto_eval actually runs,
  repo first on PATH) lacks the kit's empty-reference-file guard. Same
  routing logic otherwise.
- Pages: 9 string-collapse PUBLISHED (ARTIFACTS.md, 29 Oct). 27 not built.

## REMAINING, IN ORDER (unchanged from 29 Oct below)
1. finish coverage (3 modules)  2. copy_num per module, check each for the
rAAV-style contamination shape first  3. features_all -> module JSON ->
reinstall nine, SERIALLY  4. vocabulary rows  5. 9 registries + 9 saturations
+ 9 coverage audits, lint, publish (stop on first refusal)  6. commit, push.

---

# Parvoviridae — state as of 2026-10-29, 21:30 EDT

## WHAT IS RUNNING RIGHT NOW

`cov_run.sh` (pid 94639, launched 21:13), log `logs/cov_run4.log`, per-module
logs `logs/<Module>.eval.log`. Nine modules SERIALLY, `--jobs 4` throughout
(the curator set that after the machine kernel-panicked under a heavier load
-- do not raise it).

    Dependoparvovirus  176 genomes  <- in flight, ~125 done
    Protoparvovirus    300          queued
    Bocaparvovirus     300          queued
    Erythroparvovirus  300          queued
    Amdoparvovirus     300          queued
    TetraCopiparvovirus 300         queued
    Hamaparvovirinae   300          queued
    DensovirinaeAmbi   184          queued
    DensovirinaeMono    40          queued

Measured rate ~1.3-1.9 genomes/min. Expect completion Friday 30 Oct, late
morning to afternoon. An EARLIER estimate of 11 h was wrong -- it was never
measured. Measure before quoting.

`run_gto_eval.py` RESUMES: a genome whose `.ann.gto` and table are both
non-empty is reused. Killing and restarting the run is cheap. Kill in this
order or the children are orphaned and the load CLIMBS: cov_run.sh parent,
then run_gto_eval.py, then `pkill -f annotate_by_viral_pssm.pl`.

## DONE THIS SESSION (29 Oct)

- **124 rAAV vector genomes dropped** from the Dependoparvovirus coverage set
  on the curator's instruction. Held, not deleted, in
  `cov/Dependoparvovirus/held_rAAV_vectors/` with `DROPPED.ids` and a README
  giving the method. Identified by tblastn of this module's own NS1/VP/AAP
  collections against all 300 contigs: bimodal with an empty gap, every
  genome that calls anything is >=46.6% parvoviral protein, the other 124 have
  NO hit at all. NOT by name (all 300 are deposited as "Adeno-associated
  virus") and NOT by "called nothing" (circular).
- **TetraCopiparvovirus gained a feature: ORF3.** 16 unique sequences, every
  one exactly 204 aa, all Porcine parvovirus 4, 240/240 pairwise >=96%
  identity, zero tblastn hits to any other collection in the module. It had
  been swallowed by the global `^orf ?[0-9x]+( protein)?$` NOT_MODELLED rule.
  1 PSSM, 0 leftovers, bit_cutoff 378 (raw min = genuine min = 421, so the
  profile matches nothing but its own kind). Installed into VA_spina.
  `build_cols.py` carries a module-scoped exception so a rebuild reproduces
  it; the global discard is unchanged. Collections verified unchanged for
  NS1/VP/VP2 (515/559/104) before and after.
  Named "Uncharacterized lineage specific protein ORF3" -- no demonstrated
  function in the literature (Cheung 2010 PMID 20830210; RNA profile PMID
  21785854). It sits in the bocavirus NP1 genomic position but shares minimal
  identity with NP1, so it is NOT named NP1-like. PMIDs are in
  `PMID_claude_generated`, never `PMID`.
- **The other two dropped features stay dropped, now verified on clean data:**
  Bocaparvovirus NS2 is 4 unique sequences (all 577 aa) and TetraCopi SAT is
  1 unique (84 aa) -- both below the member floor of 5. DensovirinaeMono NS3
  is 2 unique in its own module. NOTE: synonyms.tsv counts DEPOSITS; the
  floor applies to UNIQUE sequences, and for these the two differ by 4-14x.
- DensovirinaeAmbi "ORF3" (8 unique, 533 aa) is NS1 under another name --
  96.8-99.4% identity over the full length to NS1_76. Already covered, left
  alone deliberately rather than rebinned.
- **Kit fix, committed `1e93ad08`:** `install_module.py` wrote the ~97-module
  master `Viral_PSSM.json` with `open(vp,"w")`, truncating it in place while
  annotators re-read it per genome. Now writes a sibling, fsyncs, renames.
- **Nine string-collapse pages published** and ARTIFACTS.md restructured: the
  eight Astroviridae rows had been appended to "Not tied to one taxon" and now
  sit in their own `### Avastrovirus` / `### Mamastrovirus` sections, and the
  nine Parvoviridae modules have sections of their own.

## STILL TO DO, IN ORDER

1. Wait for `cov_run.sh`. Watch `logs/cov_run4.log` for "COVERAGE RUNS DONE".
2. **copy_num.** Per-feature call rates from `cov/<mod>/tbl/*.feature.tbl`
   **column 17** (the PSSM name -- `report.tsv` carries flags only). The
   denominator is DECIDED: genomes carrying >=1 call from the module, not all
   routed genomes. Default rule `called/routed >= 0.85`; `--require` overrides
   on biology.
3. Regenerate `features_all.json` -> module JSON -> reinstall all nine.
4. Regenerate the vocabulary rows. The "Used for Genome Quality" column
   derives from copy_num and is still 0 on all 50 committed rows. ORF3 adds a
   51st; `gen_vocab.py` reads the INSTALLED JSON so it picks ORF3 up for free.
5. **27 remaining pages**: 9 PSSM registries (self-recall BLAST -- do not
   start while the coverage run has the box), 9 vocabulary saturations and 9
   coverage audits (both need the coverage data).
6. Commit and push.

### Predictions recorded BEFORE the measurement, to be checked against it

- Protoparvovirus NS1 should FALL from 99.7% -- part of that rate was calling
  genomes that are not protoparvoviruses.
- Dependoparvovirus should RISE, its Rep forms especially, now that the
  waterfowl clade reaches a module that models Rep78/68/52/40 and AAP.
- If Protoparvovirus does NOT move, the mechanism above is misread somewhere.

## OPEN, RECORDED, NOT BEING CHASED

- 565 unrouted unclassified genomes.
- DensovirinaeAmbi routes at 77.7% (curator decision to ship).
- Astroviridae ORFX: 3 deposits, below the floor of 5.
- Astroviridae strings such as `Nonstructural polyprotein nsP1a` predate the
  28 Oct functional-naming rule. Flagged for the curator, NOT changed.

---

# ---- PREVIOUS STATE (2026-10-28) ----

# Parvoviridae — state as of 2026-10-28, late session

FIRST DNA VIRUS FAMILY IN THE PROJECT. 95 existing modules are all RNA
viruses. The pipeline is PSSM-over-ORFs with blastn rep-contig routing, which
in principle does not care how a genome replicates — that is the assumption
this family tests. Nothing has broken it so far.

## DONE AND MEASURED

- Dump: 9,432 genomes / 26,892 features / 12,895 sequences / 469 strings.
  `check_dump.py` PASS. `--min-len 3500` is the measured knee, not a round
  number (see the git history of dump.sh for the per-genus cost of 4000).
- Routing gate: nine modules, seven routing >=94%, Hamaparvovirinae 90.7%,
  DensovirinaeAmbi 77.7% shipping with the tail open ON THE CURATOR'S
  DECISION: *"if they are truly diverse, we won't catch them with
  representative genomes."* Rep contigs already exist at
  `gate/<Module>/Rep-Contigs` and `gate2/Densovirinae{Ambi,Mono}/Rep-Contigs`,
  built with the same `build_rep_contigs.py` call the ship step uses.
- Collinearity measured across four genera: everything collinear, single HSP,
  full query, same frame. **Parvoviridae needs NO splice machinery.**
- Collections FINAL: 51 features over 9 modules. `UNASSIGNED_TRACKING.tsv`
  has no row at or above the member floor of 5 in ANY module — the skill's
  own finished test passes.
- `features_all.json` written for all nine.
- `annotation-vocabulary.tsv` 652 -> 703 rows, 51 Parvoviridae rows,
  generated FROM features_all.json by `gen_vocab.py` so the two cannot drift.

      Dependoparvovirus    8  VP VP2 VP3 NS1 REP68 REP52 REP40 AAP
      Protoparvovirus      7  VP VP2 NS1 NS2 NS3 SAT NP1
      Bocaparvovirus       5  VP VP2 NS1 NP1 ORFX
      Erythroparvovirus    6  VP VP2 NS1 7.5K 9.5K 11K
      Amdoparvovirus       5  VP VP2 NS1 NS2 NS3
      TetraCopiparvovirus  3  VP VP2 NS1
      Hamaparvovirinae     5  VP NS1 NP1 NS2 NS3
      DensovirinaeAmbi     7  VP VP2 VP3 VP4 NS1 NS2 NS3
      DensovirinaeMono     5  VP VP2 NS1 NS2 NS3

## THE RULE THAT CHANGED THIS SESSION

Curator, 28 Oct: *"the NP1 NP2 and nonstructural etc, are all the reason why
we need literature based functions and names that are unambiguous"*.
Saved as memory `annotation-names-must-be-functional`.

The feature KEY stays short and positional — it is only a filename. The
`anno` string states a FUNCTION from the literature, or states explicitly
that the function is unknown, in the curator's own phrasing
("Uncharacterized lineage specific 7.5 kDa protein"). All 51 strings were
rewritten on this rule; the earlier draft said things like "Replication
initiator protein" and "Nonstructural protein NS2", which is the thing being
objected to. PMIDs are in `PMID_claude_generated` in features_all.json and
NOT in the vocabulary table, which has no provenance column.

## WHY THAT RULE IS NOT COSMETIC — the measurement that proved it

Clustering every module's non-nested features at 30% identity / 70% coverage
found name/homology disagreements everywhere:

    Hamaparvovirinae  ONE cluster of 62 proteins = 43 "NS2" + 19 "NP1"
    Bocaparvovirus    44 "NS2" deposits sit inside the NS1 clusters
    Bocaparvovirus    17 "NS1" + 11 "NS2" deposits are actually NP1
    Amdoparvovirus    69 "NS3" deposits sit in the 114 aa NS2 group
    Protoparvovirus   4 "NS2" deposits belong with the NP1-like group

`build_cols.py` now runs `consolidate()`: one homology cluster, one feature,
majority label wins — **unless OVERRIDE says literature does**. The
Hamaparvovirinae case is exactly why the override exists: the majority said
NS2 and the majority was WRONG. Literature settled it (cachavirus NP = 210
aa; feline chaphamaparvovirus NP1 gene 1871..2434 = 188 aa; both inside the
cluster's 137–295 aa range). The nested capsid/Rep forms are excluded from
consolidation — they are nested by design and homology cannot separate them.

## THINGS DROPPED, WITH THE MEASUREMENT

- **Bocaparvovirus NS2** — all 4 remaining unique sequences are 100%
  identical to NS1 residues 284–860. A strict truncation of NS1, which the
  curator's checklist forbids as a profile and which no bit_cutoff can
  separate. HBoV NS2 is a SPLICED product that "shares the C terminus with
  NS1", so these deposits are the shared part annotated alone, not an
  alternative start like VP2. **Open item: revisit with splice machinery if
  wanted.**
- **TetraCopiparvovirus SAT** — 1 unique sequence (14 identical deposits).
  Cannot cluster. Open item.
- Two UNC discovery groups fell under the floor after consolidation absorbed
  them into NS1/NS3. Nothing above the floor is left anywhere.

THIN, built at `-m 2` (recorded in BUILD_PARAMS as a real departure):
DensovirinaeAmbi VP3 (4) and VP4 (4), DensovirinaeMono NS3 (2).

## PIPELINE FIXES MADE IN build_cols.py THIS SESSION

1. Homology rescue BEFORE clustering, at >=80% id over >=85% of the query —
   the same two numbers `rescue_unassigned.py` and `filter_unchar.py` use,
   which MUST agree. Caught "major caspid protein" (a typo no name rule can
   ever match) and ~60 others. Rescued capsids are re-routed BY LENGTH,
   because VP and VP2 are ~100% identical and best-bitscore cannot choose.
2. Accessory discovery on the leftover pool -> UNC1..UNCn at the floor of 5.
3. `consolidate()` — above.
4. The floor is re-applied AFTER consolidation, and it is a DISCOVERY floor:
   it binds on UNC groups only. A named, literature-supported feature with 4
   members is built at -m 2, not deleted.
5. **Stale collection files are cleared before writing.** A feature that
   stopped being declared left its fasta behind and would have been built
   into a PSSM anyway. This bug silently contaminated one QC run already.
6. `non ?structural` matched a space or nothing but NOT a hyphen. 84 misses
   in Erythroparvovirus alone.

## NEXT

1. PSSM builds — `build_all.sh`, 4 modules at a time, logs in `logs/`.
   **RUNNING / CHECK `logs/<Module>.build.log` FIRST.**
2. leftover profiles (`build_leftover_pssms.py`) — the divergent members, and
   a direct readout of what each module will miss in the wild.
3. `qc_cross_feature.py` on every module (the consolidation should have left
   it with nothing to say — that is the point of running it).
4. derive cutoffs, module JSON, install (SERIALISE: install rewrites the
   runtime Viral_PSSM.json), coverage runs, copy_num.
5. place the 872 unclassified genomes by homology (`Contigs.unassigned/`).
6. 36 artifact pages (4 x 9). Lint before publishing.
   **Publishing caution: failed publishes land in ANL security logs as
   suspected malware callbacks. Stop on the first refusal, do not retry.**

## STILL TRUE FROM EARLIER

- `~/Viral_Annotation` is READ-ONLY (perl libs only). Runtime target is
  `$LOWVAN_DATA_DIR`.
- copy_num default rule: called/routed >= 0.85, curator may override with
  `--require`. Not yet applied to any Parvoviridae module.
- `viral_genome_quality.pl` line 112 loads min_len/max_len only for features
  carrying copy_num, so non-essential features are never length-checked.
  Reported, not changed.

## FRAGMENT FILTER — added after the first PSSM run failed

The first build FAILED to produce anything for Protoparvovirus NP1: its 28
sequences made 26 clusters and none reached the member floor. That was not
biology. The sequences span 101–290 aa, and mmseqs needs 80% COVERAGE as
well as 80% identity, so the length spread — not divergence — broke it.

Measuring every collection showed the cause was family-wide: deposits that
are truncated without SAYING so. Protoparvovirus VP held a 26 aa member
against a median of 727; NS1 a 36 aa one against 665; TetraCopiparvovirus
NS1 a 29 aa one against 652. They pass the `, partial` discard, and they
pass the homology rescue because a fragment aligns over nearly all of
ITSELF, so query coverage is high by construction.

`build_cols.py` now applies a per-genus length window (0.6x–1.6x of that
genus's own median) and writes what it removes to `<mod>/outliers/<FEAT>.fasta`.

**The window is PER GENUS, and genera with under 5 members are exempt
entirely.** Both rules were learned the hard way in this session:

  * a single module-wide window wiped Diciambidensovirus out of
    DensovirinaeAmbi VP, Ichtchaphamaparvovirus out of Hamaparvovirinae VP,
    and Chaphamaparvovirus out of Hamaparvovirinae NS2.
  * a LOOSE absolute bound for small genera still wiped Diciambidensovirus.
    Its capsid is 202–205 aa against a module median of 759 and it is REAL:
    two independent species (Bactericera trigonica and Diaphorina citri
    densovirus) agree on the size, they are 58% identical to each other over
    113 residues, and they align to a 242 aa VP2 member. Complete short
    proteins, not truncations. They stay in, fail -mi 0.8 as they should,
    and the LEFTOVER pass is their designed home.

Final check: **0 genus wipe-outs across all 51 features.**

## THE 872 UNCLASSIFIED GENOMES — measured, and mostly a gap

Routed each against the nine modules' rep contigs (megablast, e<=1e-10,
best bitscore wins):

    routed    307 / 872  (35.2%)
    unrouted  565

    DensovirinaeAmbi 115 | Protoparvovirus 79 | TetraCopiparvovirus 42
    Hamaparvovirinae  40 | Bocaparvovirus 14  | DensovirinaeMono 8
    Dependoparvovirus  5 | Erythroparvovirus 4

565 genomes — 6% of the family — reach no module. Routing is nucleotide
BLASTn, so this is the divergent-tail signature and MORE REFERENCES IS THE
WRONG INSTRUMENT; see memory `divergent-tail-is-not-a-reference-count`.
Recorded as an open gap, not chased. Placement table: `unassigned_place/placed.tsv`.

## 2026-10-29 — COVERAGE RUNS, and the rAAV vector finding

Coverage runs are task `b9hq4le4v`, ~4 h per 300-genome module at --jobs 4.
Done: Dependoparvovirus (01:22), Protoparvovirus (05:27). Bocaparvovirus
running. Six to go. Logs `logs/<Mod>.{gto,eval}.log`, output under `cov/<Mod>/`.
Per-feature call rates are read from `cov/<Mod>/tbl/*.feature.tbl`, column 17
(`<Module>.<FEAT>.<n>.pssm`) -- NOT from report.tsv, which only carries flags.

### Dependoparvovirus: 124 of 300 sampled genomes are RECOMBINANT AAV VECTORS

They route to the module correctly and then call nothing:

    "Routed to Dependoparvovirus; no profile above threshold"

That is the right answer, not a defect. Measured:

  * `10804.1490` opens `ttggccactccctctctgcgcgctcgctcgctcactgagg...` -- the AAV2 ITR.
  * no-call genomes: longest ORF median **1482 aa** (max 1594) against 735 aa
    for Rep78 in the called genomes. A payload, not a viral gene.
  * 20 ORFs >=200 aa from 10 no-call genomes, blastp vs ALL Dependoparvovirus
    collections, e<=1e-5: **zero hits**.
  * names do not separate them -- they are deposited as "Adeno-associated
    virus", "Adeno-associated virus - 2", "adeno-associated virus 2", the same
    strings the real genomes carry.

**Consequence for copy_num:** the denominator must be genomes that are actually
dependoparvoviruses, not genomes that routed. On the routed denominator every
Dependoparvovirus feature falls under 0.85 and NOTHING gets copy_num, which
means `viral_genome_quality.pl` can never flag a bad dependoparvovirus genome.
Rates on the 176 genomes carrying any viral call:

    VP 100.0  VP2 99.3  VP3 97.8  AAP 97.8  NS1 92.8
    REP40 73.9  REP52 63.0  REP68 68.1      (nested forms, expected lower)

Protoparvovirus, for contrast, is clean: 294/299 routed, NS1 99.7, VP 99.7,
VP2 98.6, NS3 73.1, SAT 72.8, NS2 2.4, NP1 0.3.

**Check every other module for the same shape before applying copy_num** -- any
taxon used as a vector or a cloning backbone will carry the same contamination.

## 2026-10-29 09:10 — ARTIFACT PAGES, 9 of 36

String collapse: all nine rendered into `artifacts/<module>-string-collapse.html`,
`check_artifact_html.py` clean, NOT yet published.

Three defects fixed to get there; each would have put a wrong number on a
published page.

1. **synonyms.tsv was in a private schema.** build_cols.py wrote
   module/feature/string/count. The shared reader wants
   count/annotation/genus/bins_to with the UNBOUND rows present. No genus
   column means a string that means different things in two genera collapses
   to one row; no unbound rows means the collapse ratio has no denominator.
   Now emitted per module through `synonyms_table.write`.
2. **`collect_synmap.py` read from `workdir` and wrote to the CWD.** Nine
   modules in a row left ONE set of agg/unbinned/typos files and eight silent
   overwrites. Patched IN THE KIT to write beside its inputs
   (`skill/scripts/collect_synmap.py`) -- uncommitted.
3. **The vocabulary could disagree with the collections.** `consolidate()`
   moved sequences between features and left their strings behind, so five
   phantom collapse targets appeared: UNC1 in Hamaparvovirinae and
   DensovirinaeAmbi, plus the three features dropped on measurement
   (Bocaparvovirus NS2, TetraCopiparvovirus SAT, DensovirinaeMono NS3).
   `hits` entries are now `(genome, seq, annotation, genus)` and synonyms.tsv
   is DERIVED from hits at the end instead of accumulated beside it. Dropped
   features bind to NOT_MODELLED. All nine agg.json target sets now equal
   their features_all.json key sets exactly.

`build_cols.py --no-write` recomputes every count and leaves the collections
alone -- the PSSMs are built from those files, so regenerating a report must
not touch them. Verified by md5 before and after.

Collapse measured: Protoparvovirus 118 strings -> 7 features (16.9x),
Bocaparvovirus 77 -> 5 (15.4x), TetraCopiparvovirus 41 -> 3 (13.7x),
Hamaparvovirinae 66 -> 5 (13.2x), Erythroparvovirus 64 -> 6 (10.7x),
Dependoparvovirus 71 -> 8 (8.9x), DensovirinaeAmbi 59 -> 7 (8.4x),
Amdoparvovirus 40 -> 5 (8.0x), DensovirinaeMono 30 -> 4 (7.5x).

### NOT DONE
- registry pages (9) -- `collect_registry.py` does self-recall BLAST per
  module. HELD so it does not compete with the coverage runs for CPU.
- vocabulary saturation (9) and coverage audit (9) -- both need the coverage
  runs.
- copy_num. Denominator decision 29 Oct: genomes carrying >=1 call from the
  module, NOT all routed genomes (see the rAAV section above).
- PUBLISHING. Nothing published yet. Stop on the first refusal.

## 2026-10-29 12:40 — THE MISFILING FIX, and what it invalidated

### What was wrong

`build_cols.py` assigned a protein to a module by the **BV-BRC genus of its
genome**. A genome filed under the wrong genus damages the module twice:

  1. its proteins join the wrong module's COLLECTIONS (wrong training data);
  2. it becomes eligible as a REP CONTIG, and a rep contig is a BLASTn
     subject, so it wins routing for the clade it really belongs to.

Both are invisible downstream because the profiles and the references end up
wrong in the SAME direction and corroborate each other.

**20 genomes, 17 of them filed as Protoparvovirus.** Worst case
`Protoparvovirus sp. HNXY21` (taxid 1987106, BV-BRC's catch-all
"Protoparvovirus sp."): a waterfowl dependoparvovirus, 98.5% identical over
5,073 nt to goose parvovirus, chosen as a Protoparvovirus rep contig. It took
**37 of 300 sampled dependoparvoviruses** with it. Each was annotated with
three generic Protoparvovirus features instead of the eight Rep/VP/AAP
features Dependoparvovirus models. No error anywhere; the calls simply came
from the wrong module.

Not caught by the routing gate because each module was measured against its
OWN rep contigs in isolation (Protoparvovirus scored >=94% there). The error
only exists when all 97 modules compete in one repo.

### What was NOT rebuilt, and why — CURATOR'S CALL, 29 Oct

**No PSSM was rebuilt.** Every profile inside a feature carries that feature's
single annotation string, so a Protoparvovirus NS1 profile built from goose
parvovirus can only ever emit Protoparvovirus's own NS1 annotation -- it
cannot leak a dependoparvovirus name. Once routing is correct those profiles
are inert. Verified with the real caller: all 38 reassigned sequences are
already called by the EXISTING destination profiles, every one far above its
cutoff (Dependoparvovirus VP scored 1583 against a cutoff of 420).

Protoparvovirus therefore still holds ~14 profiles trained wholly or partly on
foreign proteins (NS1 17/20/22/23.fa, VP 17/21/23/25.fa, VP2 15/17/19/20/21.fa,
NP1 3.fa). Harmless for annotation. They will look odd in the registry page's
self-recall. Deleting just those is cheap if ever wanted. **Open, deliberate.**

### Tools: new and changed

- `find_misfiled.py` NEW. All-vs-all blastp on the NS1 anchor; a genome is
  misfiled when its >=90%-identity relatives overwhelmingly live in another
  module. **Best-single-hit with a margin DOES NOT WORK** -- a cluster of
  misfiled genomes vouches for itself; that test found 8 of 18 and missed
  HNXY21, the worst one. `--from <snapshot>` re-runs against pre-fix
  collections. **It MERGES into misfiled.json, never replaces**: once applied,
  the evidence is gone, and one overwrite dropped 18 fixes down to 2.
- `converge.sh` NEW. Rebuild -> detect -> repeat to a fixed point. Needed:
  pass 1 found 16 in Protoparvovirus, and removing them exposed 2 more that
  those 16 had been propping up. Converged at 20.
- `gate_fix.sh` NEW. Rebuilds rep contigs with assignment going through
  misfiled.json. Output in `gate3/`. Supersedes gate.sh / gate2.sh.
- `fix_bounds.py` NEW. **Reproduces gen_module_json.py's 6% PAD.** A first
  draft set raw collection min/max and proposed 100 changes, ALL narrowings,
  including in modules that only GAINED sequences -- that asymmetry is the
  tell. Shipped, it would have flagged correctly-modelled proteins as too
  short (the Velarivirus failure).
- `build_cols.py` — `modof()` resolves genome->module through misfiled.json,
  and WARNS LOUDLY if the file is absent. Also `--no-write` (recompute counts
  without touching collections) and synonyms.tsv now DERIVED from hits.
- KIT, uncommitted: `collect_synmap.py` wrote to the CWD while reading from
  `workdir` (8 of 9 modules silently clobbered); `gen_module_json.py` dropped
  `replicon_geometry` on the `--segments` path -- 187 segments in the runtime
  carry "linear" and the only 9 without it were these nine modules.

### State

- collections rebuilt, converged, 20 differ from the pre-fix snapshot
  (Protoparvovirus -17 NS1 / -12 VP / -10 VP2 / -2 NP1; the rest gain).
- rep contigs rebuilt, all nine. **0 misfiled genomes serve as a wrong-module
  reference in the runtime; HNXY21 is in no .dna file.**
- 35 of the 38 hijacked genomes now route to Dependoparvovirus. The 3 that do
  not are the divergent tail: two hit NOTHING in any module, one
  (Yunnan rodent dependoparvovirus 2) leans weakly to Protoparvovirus on real
  but poor similarity. Pre-existing gap, not a regression.
- bounds: 3 changed, each traceable to a moved genome. Protoparvovirus
  cutoffs re-derived (NS1 136 -> 148).
- all nine install with **0 problems**; 97 modules in the runtime.
- **coverage runs RESTARTED** on the corrected repo, same deterministic
  genome sample (contigs copied from the old run, not re-drawn, so the two
  are comparable). Old results preserved in `cov_old/`. ~20 h at --jobs 4.

### PREDICTIONS, recorded BEFORE the measurement

- Protoparvovirus NS1 should FALL from 99.7% -- part of that rate was calling
  genomes that are not protoparvoviruses.
- Dependoparvovirus should RISE, its Rep forms especially, now that the
  waterfowl clade reaches a module that models Rep78/68/52/40 and AAP.
- If Protoparvovirus does NOT move, the mechanism above is misread somewhere.

### STILL OPEN
- copy_num. Denominator: genomes carrying >=1 call from the module, NOT all
  routed genomes (rAAV vectors, see the 29 Oct section above).
- 27 of 36 artifact pages. 9 string-collapse pages are rendered and lint
  clean in `artifacts/`, NOT published. Nothing published yet.
- re-check the three dropped features (Bocaparvovirus NS2,
  TetraCopiparvovirus SAT, DensovirinaeMono NS3) against the CORRECTED
  collections -- all three were judged on contaminated data and
  TetraCopiparvovirus gained 15 sequences.
- vocabulary commit held for curator review of the 50 naming strings.

## 2026-10-29 13:05 — FOREIGN PROFILES DELETED (curator-approved, with a limit)

Curator, 29 Oct: *"You can delete any profile that is in the wrong module and
was built from the wrong organisms ... and doesn't cause a loss in the
matching. My comment was mainly on the fact that if a profile was mixed, there
was no sense rebuilding it."*

So: delete only profiles built ENTIRELY from misfiled genomes, and only where
deleting costs no match. Both halves bind.

  * 14 Protoparvovirus profiles contain a misfiled genome; **13 are wholly
    foreign**, 1 is mixed (VP2/19.fa, 1 of 2) and stays.
  * Deleting all 13 WOULD have cost 4 matches: VP2 1511906.642, 10788.3822,
    1511906.641 and NP1 10788.3858 drop below bit_cutoff. Those profiles are
    the sole coverage for genuinely divergent Protoparvovirus members that
    happen to resemble them.
  * **10 deleted, 0 sequences lost.** Kept as sole coverage: NP1.3, VP2.20,
    VP2.21.

      deleted   NS1 17,20,22,23 | VP 17,21,23,25 | VP2 15,17

Alignments and PSSMs archived to `deleted_profiles/`, not destroyed. Repo
Protoparvovirus 91 -> 81 profiles; `install_module.py --check` still 0
problems. The test that decides this is in /tmp/protoscores.json (every
Protoparvovirus collection sequence scored against every profile); regenerate
it with psiblast against collections/Protoparvovirus if needed.

**Coverage restarted 13:05 on the FINAL repo** (logs/cov_run3.log). The two
earlier attempts are void: the first measured the pre-fix repo, the second was
killed 15 min in because the profile deletion would have changed the repo
underneath it.

### NOTED, NOT YET CHASED
**11 of 787 Protoparvovirus NS1 collection sequences are called by NO profile**,
before or after the deletion. Unrelated to the misfiling. This is the
leftover-coverage signature from the skill's Alpharhabdovirinae P case, where
the uncovered population and the uncalled genes were the same set. Check
whether build_leftover_pssms.py reached them.
