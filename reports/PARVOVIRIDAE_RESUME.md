# Parvoviridae — nine modules, shipped; four things open

State of play as of **2 November 2026**. Point a fresh session at this file.

This is the current state, not a history. The chronological record — every
dated decision, every failed build and why — is `RESUME_HERE.md` in the build
directory, newest block first, mirrored to
`CEPI/Parvoviridae/RESUME_HERE.md`. The append-only run log beside it is
`PROGRESS.md`.

**This is the project's first DNA family.** Everything before it was RNA
viruses, so the splice machinery, the ambisense handling and the nested
alternative-start forms were all exercised here for the first time.

## Where things are

| thing | path |
|---|---|
| the kit (the repo that matters) | `LowVan_Skill/`, remote `jimdavis1/LowVan_Skill` |
| per-module JSON, PSSMs, rep contigs | `LowVan_Skill/modules/<Module>/` |
| splice references | `LowVan_Skill/Splice-Variants/<Module>/<FEAT>.fasta` |
| published pages, index | `LowVan_Skill/reports/`, `LowVan_Skill/ARTIFACTS.md` |
| build scratch, collections, measurements | `~/lowvan_scratch/parvo/<Module>/` |
| the test repo coverage runs against | `~/lowvan_scratch/VA_parvo_route/` (176 refs) |
| cut-down repos for cheap re-runs | `~/lowvan_scratch/VA_only_<Module>/` |

**Push only to `jimdavis1/LowVan_Skill`.** `~/Viral_Annotation` and its remote
are off-limits. `Alsuviricetes/Viral_Annotation` is referenced for its perl
libs only — per the curator, 7 October: **do not write to it**.

## Status

**All nine modules are built, evaluated, installed, shipped and pushed**, each
as its own commit, with four pages each — 36 of 36, none missing, all indexed
in `ARTIFACTS.md`. Last commit `0e111472`.

| module | feats | specials | refs | exemplars | routed here | unrouted |
|---|---|---|---|---|---|---|
| Amdoparvovirus | 5 | NS2, NS3 | 2 | 42 | 41 | 1 |
| Bocaparvovirus | 5 | — | 26 | 74 | 73 | 1 |
| DensovirinaeAmbi | 8 | VP1_SPLICED | 26 | 82 | 53 | 11 |
| DensovirinaeMono | 4 | — | 12 | 24 | 24 | 0 |
| Dependoparvovirus | 8 | — | 26 | 54 | 52 | 0 |
| Erythroparvovirus | 6 | — | 15 | 30 | 30 | 0 |
| Hamaparvovirinae | 5 | — | 26 | 108 | 95 | 3 |
| Protoparvovirus | 7 | NS2, NS3 | 26 | 42 | 38 | 4 |
| TetraCopiparvovirus | 4 | — | 26 | 41 | 41 | 0 |

Five splice-reference files ship and are tracked:
`Amdoparvovirus/{NS2,NS3}.fasta`, `Protoparvovirus/{NS2,NS3}.fasta`,
`DensovirinaeAmbi/VP1_SPLICED.fasta`.

### What the three curator directives bought

1. **Coverage cutoff made live.** The 23 September caller computed
   `cov = length(hseq)/(q_to-q_from+1)` — the HSP's own span, always ~1.0, so
   the threshold was inert. The 26 October kit caller uses the **profile**
   length as the denominator. The test repo was running the old copy;
   syncing it is what made the gate work at all. Separately,
   `gen_module_json.py` had `coverage_cutoff` hardcoded at 0.65 and silently
   discarded any per-feature value — fixed in `40a1f8fd`. Used on
   DensovirinaeAmbi VP at 0.80.
2. **NS2/NS3 moved to the splice module** where the literature backs it.
   Amdoparvovirus same-coordinate calls 274/300 → **0**, clean genomes
   7 → **262**. Protoparvovirus NS1+NS3 210/299 → **0**, 285/297 clean.
   DensovirinaeAmbi NS1+NS2 11 → **0**. DensovirinaeAmbi NS1/NS2/NS3 were
   measured to be separate collinear ORFs and correctly stay PSSM-called.
3. **Routing cut to close relatives.** ~515 → ~120 vouchers per genome.
   Three modules (~770 genomes plus pass-5 re-annotation) in one boot with no
   panic. 2,000 genomes in a single boot still does not fit.

### The machine

The panic is `Cannot grow ipc space beyond IVAC_ENTRIES_MAX. Some process is
leaking vouchers @ipc_voucher.c:573`. Watch it with:

```bash
zprint | awk '/^ipc.vouchers /{print $7}'
```

Non-Apple binaries (blastn, the perl caller) leak vouchers that are **not**
reclaimed while the run is in flight. Three panics before this was understood.

**New, 1 November: they do reclaim once the runs end.** The counter peaked at
235,000 during the three-module day, settled to 198,925, and was back to
**416** — baseline — the next morning **without a reboot**. So a re-run no
longer requires a restart first; it requires the previous run to have finished.
Budget roughly 120 vouchers per genome against a 245,000 stop.

## Open, in priority order

### 1. Amdoparvovirus NS3 is called short in AMDV — **new, 1 November**

**This is a correctness bug in a shipped module, and the evidence is settled.**

Huang, Luo, Cheng, Best, Bloom & Qiu, *Virology* 2014;452-453:23-31, **PMID
24606679** — already cited on both features — gives the sizes outright:

> "NS3 is a small protein with only **87 amino acids**, in which 60 amino acids
> from the N-terminus are identical to those of the NS1 and NS2. Hence, NS3
> only has a unique sequence of 27 aa."

and builds its antibodies against "the unique C-terminus of the AMDV NS2
(**aa 61-114**) and NS3 (**aa 61-87**)". Qiu 2006 (PMID 16378968) names only
NS1, NS2 and NS3 for AMDV. **There is no fourth small nonstructural protein.**

The 87/90/92 aa product previously logged as an unnamed "third splice product"
is therefore AMDV's own NS3, at exactly the published length. Measured from the
collections:

| length | n | taxa |
|---|---|---|
| 114 aa — NS2, correct | 63 | 28314 |
| 87/90/92 aa — **NS3, correct** | 33 | 28314 only |
| 68/71 aa | 52 | 2 AMDV (truncated), 50 other amdoparvoviruses |

- The long and short forms are **the same reading frame** — both read
  `...TTSSDLR | LCDCNKQH...` after the shared 60-residue exon. The short ones
  stop early; the 87 aa ones continue `QHNQSNRWMCGNK...` to the real terminus.
- **No AMDV genome deposits both** — 0 of 69. The two short AMDV deposits are
  truncated annotations.
- 28 of the 32 long deposits are **labelled NS2**, so they binned into the NS2
  collection and were never available to build NS3.

Consequences in the shipped module (`3030543f`, `47f8800d`):

1. `NS3` declares `min_len 60, max_len 80`, so the real protein is **outside
   the window and uncallable** — and unlike a PSSM call, min/max *is* enforced
   on a splice feature. Measured on the shipped 300-genome run: NS3 is called
   46 times and **every single call is 59–72 aa** — the longest anywhere in the
   module is 72, and not one is 87–92. In AMDV itself it fires in **10 of 242
   genomes (4%)**, against 19 of 25 for taxon 2169703 and 11 of 12 for 1513315.
   NS2 by contrast is healthy: 225 of its calls are 114 aa, the published
   length.
   *(An earlier revision of this file said "10 of 41 exemplars". That came from
   `coverage/Amdoparvovirus.json`, which is built from the 42-exemplar
   **pre-splice** run of 31 Oct and does not describe the shipped module. The
   figures above come from `cov/Amdoparvovirus/tbl`, the post-splice run the
   module shipped on. Check which run a number came from before quoting it.)*
2. All 11 AMDV references in `Splice-Variants/Amdoparvovirus/NS3.fasta`
   translate to **66 or 68 aa** — they encode the truncated acceptor. Verified
   by translating every reference: `28314 {66: 1, 68: 10}` against
   `2169703 {68: 1, 71: 14}` and `1961072 {71: 5}`.
3. NS2 itself is **fine**. Its 227 references are 195× taxon 28314 and all
   translate to 114 aa; the mislabelled deposits did not leak in.

**The fix, in order:**

1. Re-bin: move the 87/90/92 aa deposits out of
   `held_pssm_features/NS2.collection.fasta` into the NS3 collection. They are
   `28314.*|NS2_116..144` and `28314.*|NS3_1..4`.
2. Widen `NS3` to `min_len 59, max_len 95` in `features_all.json` — spans the
   66 aa of *Carnivore amdoparvovirus 2/3* through AMDV's 92.
3. Rebuild `Splice-Variants/Amdoparvovirus/NS3.fasta` with
   **`build_refs_from_deposits.py`**, which places the junction so the spliced
   translation *is* the deposit. Do **not** use `build_splice_refs.py`: a
   canonical donor/acceptor search is what picked the short acceptor in the
   first place.
4. Re-run coverage against `VA_parvo_route`, re-run pass 5, republish the four
   Amdoparvovirus pages, re-commit.

**The page's explanation matters more than its figure.** The published audit
currently reads the low NS3 rate as biology: "NS3's 15.9% is the honest figure.
It is expressed at much lower levels than NS2, it is deposited by only a
minority of submitters… Where NS3 is called, the junction was found; where it
is not, nothing was invented." The expression claim is true of the protein and
false as an explanation of *this number*. The rate is low because AMDV's NS3 is
87 aa and the feature declares `max_len 80`. **A measurement artefact is
currently published as a finding about the virus**, and that paragraph has to be
rewritten along with the bar.

**Related defect, same run:** the splice caller writes gene symbol
`nonstructural` for both NS2 and NS3 rather than `ns2`/`ns3`, so the two are
indistinguishable by symbol in the feature table and can be told apart only by
the annotation string in column 18. `evaluate_coverage.py` keys on the symbol,
which is why `cov/Amdoparvovirus` shows 290 calls named `nonstructural` and no
`ns2` or `ns3` at all.

### 2. DensovirinaeAmbi VP1 and VP2 at identical coordinates, 40 of 171

Down from 46, not solved. Raising VP1's coverage cutoff to 0.80 did what it was
meant to — it rejected the partial match that was being extended over the cap2
frame — and a different profile then won at full coverage. That profile,
`VP.12`, is built from two 578 aa deposits of taxon **2050976**, whose VP
deposits run 536–942 aa in clusters near 577, 645, 759 and 823: a nested set
filed under one name.

**The "bin by size per species" recipe cannot be done — measured 2 November.**
The only two taxa with any capsid spread are the two whose species label is
*incertum*, meaning uncertain, so there is no species structure inside them to
bin on:

| taxon | n | range | spread | species label |
|---|---|---|---|---|
| 2050976 | 30 | 536–942 | **406** | *Protoambidensovirus* **incertum3** |
| 2030127 | 21 | 508–1168 | **660** | *Protoambidensovirus* **incertum5** |
| every other taxon | 1–14 | — | **0–110** | properly classified |

All 50 genomes of taxon 2050976 carry the single label *Protoambidensovirus
incertum3*. The two together are **51 of 86 VP deposits** and **45 of the 71
exemplars that carry features**, so excluding them shrinks the module rather
than fixing it.

Three real options, and the choice is a taxonomy call for the curator:

1. **Bin by size within the taxon** and name by size rank — clusters sit near
   577, 645, 759 and 823. This assigns capsid names on an assumption, which is
   the same move that [[deposited-name-is-a-claim]] warns about, in the other
   direction.
2. **Collapse them into one uncharacterized capsid bin**, claiming no VP1/VP2
   distinction where the data does not support one.
3. **Exclude the two incertum taxa** from the VP profiles and accept a much
   smaller, better-grounded module.

The superseded recipe was: re-bin by size per species using `P.species.tsv`,
rebuild the VP profiles, re-run the 171-genome panel against
`~/lowvan_scratch/VA_only_DensovirinaeAmbi` (35 refs, ~21,000 vouchers). **The
ambisense genera is a curation call and wants the curator's eye.**

### 3. Dependoparvovirus AAP N-termini — cause recorded, caller unchanged

34 of 34 deposits start with M; **1 of 169 calls does**. The cause is that AAP
does not start at an ATG. Three verified citations are now in the module JSON
under `PMID_claude_generated`, categorised (commit `0e111472`):

- **coordinates** — **20479244** (Sonntag, Schmidt & Kleinschmidt, *PNAS*
  2010): AAP is a nested alternative ORF of *cap* and "translation is
  initiated at a nonconventional translation start site". The same paper
  warrants declaring AAP genus-wide: "sequence analysis suggests that also all
  other species of the genus *Dependovirus* encode a homologous protein in
  their cap gene".
- **coordinates** — **39205189** (Kuz, McFarlin & Qiu, *Viruses* 2024): names
  the codon — "both AAP and MAAP translations are initiated at non-canonical
  start codons, **CUG** for AAV2 and AAV5 AAP".
- **function** — **21917944** (Sonntag et al., *J Virol* 2011): AAP promotes
  capsid assembly across AAV1, 2, 5, 8 and 9.

**Decided and applied, 2 November, on the curator's instruction: turn the
scan-to-Met-start off.** AAP is now `upstream_ext: 0`. That flag gates exactly
one thing and nothing else — the caller's own comment says so:

```perl
# When the pssm match doesn't start with an M, it will search to the left until it finds the first one.
# This must be turned off if there is a non-AUG start.
if (($hseq !~ /^M/i) && ($upstream_ext))
```

339 features across the repo already set it to 0, so this is well-trodden. The
"too long" flag on 77 of 176 was `upstream_ext` walking on looking for a start
codon that is not there; widening the window would have hidden that rather than
fixed it.

**Not** also set: `start_to_met`, which rewrites the first residue to M and is
arguably right for a CUG initiator, since a CUG used as a start is decoded as
Met and all 34 deposits begin with M. No feature anywhere in the repo uses it.
Measure the re-run first, then decide.

### 4. Smaller, undiagnosed

- ~~Two Erythroparvovirus genomes make a call with an empty feature name~~ —
  **diagnosed and fixed in the tool, 2 November, commit `e519511d`.** Not a
  caller bug: the annotator correctly wrote a `no_features_called` sentinel row
  for `10798.2563` and `1511900.137` ("Routed to Erythroparvovirus; no profile
  above threshold") and `evaluate_coverage.py` counted it as a call. Chasing it
  turned up two more defects in the same few lines — `routed_no_calls` was not
  a category, and `refs_used` was `len(flags)` reading the wrong column.
  **Five published pages state a wrong routing number as a result** and have
  not yet been republished:

  | module | page says ships/used | measured |
  |---|---|---|
  | Amdoparvovirus | 1 / 0 | 1 / **1** |
  | DensovirinaeAmbi | 25 / 1 | 25 / **24** |
  | Dependoparvovirus | 25 / 15 | 25 / **13** |
  | Hamaparvovirinae | 25 / 28 | 25 / **24** |
  | Protoparvovirus | 25 / 0 | 25 / **20** |

  Hamaparvovirinae published more references used than it ships, which cannot
  be true of anything. The nine coverage JSONs are already backfilled from the
  surviving feature tables, so regenerating these pages needs **no
  re-annotation** — just `pass5_pages.sh` and a republish.
- **Residual capsid tails**: Bocaparvovirus VP1+VP2 at identical coordinates in
  14 of 299, Protoparvovirus 7 of 299.
- **Two Dependoparvovirus exemplars misroute**, one to Hamaparvovirinae and one
  to Protoparvovirus.
- **124 of 300 "Dependoparvovirus" genomes are rAAV vectors.** They route and
  call nothing, and they inflate every denominator on that module's pages. Not
  a defect — a property of what BV-BRC holds — but do not read the rates
  without it in mind.

## Traps this family taught, the hard way

Each of these cost a wrong result before it was understood.

- **Contig FASTA is soft-masked lower case.** An upper-case codon table turns
  every translation into `XXXX…`, silently. A builder reported "0 references,
  8 no junction reproduces the deposit" and it read as a real negative. One
  `.upper()` in `read_fasta` turned it into 8 references.
- **`build_splice_refs.py --nudge 12` (the default) is far too tight.**
  Amdoparvovirus NS2 gave 14 references at 12 and 227 at 24; Protoparvovirus
  NS2 gave 1, then 162.
- **A canonical GT can beat the real donor.** Protoparvovirus NS2 built 10–11
  residues short through three canonical builds because the true donor is
  **GC–AG** with a canonical GT inside the same window. `--noncanonical` does
  not help — it only fires when no GT..AG exists at all. Comparing the spliced
  product against that genome's own deposit is what found it, 22 of 24 exact.
- **Score the exon you are choosing, not the whole protein.** With a shared
  donor, whole-protein identity is ~0.77 for any acceptor, so it cannot
  discriminate. Amdoparvovirus NS3 picked a 428 nt intron in 135 of 173 genomes
  against the 1,333 nt the two-HSP build measured. Scoring **exon B only**
  (`k = sd // 3`) fixed it.
- **Ambisense genomes put the capsid on the minus strand.** Flip the contig
  when the minus-strand bit score wins, or the builder reports "fewer than two
  plus-strand HSPs" and looks like missing data.
- **Overlapping exons need slack.** Scindoambidensovirus VP1's two exons
  overlap by 26 query residues, putting the junction 78 nt outside the default
  `--slack 18`. At `--slack 90`: 13 references, all 816 aa, intron 175 nt.
- **`commit_parvo.sh` must add `Splice-Variants/<M>/` explicitly** — it sits
  outside `modules/<M>/`, so Amdoparvovirus first shipped without the
  references its special features are called from. Fixed in `47f8800d` and in
  the script.
- **Killing `run_gto_eval` orphans its children.** After a `TaskStop`, also
  `pkill -f annotate_by_viral_pssm.pl` and `pkill -f New-annotate-viral-taxon`,
  or load climbs instead of falling.
- **Git on Box: never run a bare `git status`.** The full-worktree refresh
  hangs, and so does a path-scoped one over `modules/`. Add by explicit path
  and commit through plumbing — `write-tree`, `commit-tree`, `update-ref`.
  `git rev-parse HEAD` and `origin/main` are instant and are enough to confirm
  a push landed.

## Restarting after a crash

```bash
# 1. where did it stop?
tail -20 ~/lowvan_scratch/parvo/PROGRESS.md
head -60 ~/lowvan_scratch/parvo/RESUME_HERE.md     # newest block first

# 2. is the repo whole? (never a bare `git status` here)
cd ~/Library/CloudStorage/Box-Box/1_Projects/CEPI/Alsuviricetes/LowVan_Skill
git rev-parse --short HEAD origin/main             # equal == everything pushed

# 3. is anything still running?
pgrep -fl 'run_gto_eval|annotate_by_viral_pssm' | grep -v pgrep

# 4. how much voucher headroom?
zprint | awk '/^ipc.vouchers /{print $7}'          # vs a 245,000 stop
```

Work is logged as it happens with `~/lowvan_scratch/parvo/progress.sh "<line>"`
and each module is committed and pushed the moment its pages publish, so a
crash costs at most the module in flight. Nothing restarts automatically on
reboot — that is deliberate, at the curator's instruction.
