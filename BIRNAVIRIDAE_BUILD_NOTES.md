
## Length-window audit against the literature, before copy_num

Run 20 Nov 2026 on the 47-exemplar coverage set (`coverage_b4`), after the POLY
single-sequence profiles. Every call's translated length was compared against
the literature window in `length_refs.json` — **not** against the generator's
measured window in the JSON, which is wider on every feature and would have
hidden all of this.

```
feat    n    min   max    literature   JSON (generated)
PEP46  42     34    46      36- 52        38- 48
POLY   43    750  1143     940-1160      883-1207
PVP2   47    395   521     480- 560      456- 577
VP1    46    833  1045     780-1060      737-1060
VP2    47    336   461     395- 475      391- 495
VP3    36    175   264     215- 310      186- 301
VP4    38    150   243     165- 295      199- 257
VP5    25    113   148      92- 160       98- 163
X       2     71   106      58- 120       66- 112
```

VP1, VP5 and X are clean against the literature. The other six have calls below
their literature floor, and **the adjacency test separates them into four
different causes.** Length alone would have mixed all four.

### 1. Partial genomes — benign, no action

| genome | segA length | POLY call | verdict |
|---|---|---|---|
| `10995.6189` | 2261 nt | 1..2259 | runs the whole contig |
| `10995.6191` | 2251 nt | 1..2250 | runs the whole contig |
| `1230254.7`  | 2778 nt | 3..2777 | runs the whole contig |
| `2920476.3`  | 2354 nt | — | pVP2 starts at nt 3 |

A full birnavirus segment A is ~3,100-3,500 nt. These four are incomplete
deposits, and the short POLY/pVP2/VP2 calls follow from the contig, not from
the profile. The genome is what is truncated. When `copy_num` goes on these
features these genomes will be flagged Poor, which is the correct reading.

### 2. A genuine lineage length — widen pVP2

`2583366.4` and `2926591.3` both call pVP2 at **474 aa**, six residues under
the literature floor of 480. Adjacency confirms it is genuine: VP2 ends at
1418/1424, pep46 begins at 1419/1425 — **flush** — and pVP2 ends 66 nt beyond
pep46, which is exactly the span the remaining C-terminal peptides occupy. The
480 floor comes from IBDV (512) and IPNV (508) only; it was never measured on
this lineage. **Window 480 -> 470.**

### 3. Clipped calls — the collection, not the window

Each of these has a nucleotide gap to its neighbour, so widening the window
would bury the defect:

| genome | feature | called | gap to neighbour | implied true length |
|---|---|---|---|---|
| `321302.3`  | pep46 | 34 aa | 21 nt to VP2 end | ~41 aa |
| `595652.4`  | VP4   | 150 aa | 174 nt before, 201 nt after | ~275 aa |
| `2556101.4` | pVP2  | 475 aa | 195 nt to VP4 | ~540 aa |
| `2556101.4` | VP3   | 175 aa | 201 nt to VP4 end | ~242 aa |
| `2920479.4` | pVP2  | 504 aa | 252 nt to VP4 | ~588 aa |
| `2920479.4` | VP3   | 212 aa | 393 nt to VP4 end | ~343 aa |

Every implied length lands inside the literature window, which is what makes
these clipping rather than lineage variation. `2556101.4` (Tetranychus) and
`2920479.4` (Entomobirnavirus-like) are the two most divergent genomes in the
set and are the same two flagged elsewhere in these notes — the leftover
profiles for PVP2, VP4 and VP3 have not been built, and that is the fix.

### 4. A pVP2/VP2 coordinate collision — a real defect

`2920474.4` and `2920476.3` each call **pVP2 and VP2 at identical
coordinates** (201..1385 and 3..1187, both 395 aa). pVP2 must be longer than
VP2 by construction — VP2 is pVP2 minus its C-terminal peptides — so an
identical pair means one profile matched the other's extent. At 395 aa the
call is in VP2's window and 85 residues under pVP2's, so **the pVP2 call is the
wrong one**. This is the nested alternative-END trap the skill documents for
capsids. In `2920474.4` pVP2 also starts 165 nt after the POLY call, so it is
clipped at both ends.

**Not yet fixed.** The window must not be widened to admit it.

### Leftover profiles built, 20 Nov 2026

`build_features.py` prints `re-run build_leftover_pssms.py for this key
afterwards` and it had never been run for this module. The leftover pool is not
a remainder — it is the set of sequences that **failed** the feature's own
clustering identity floor, which is to say the divergent lineages, which is
exactly the population the clipped calls in the window audit above belong to.

```
PVP2    3 profiles      VP4     3
VP2     3               PEP46   3
VP3     3               VP5     7
                        X       2
```

24 profiles over 25 sequences, so almost all are single-sequence and their
cutoffs must be measured rather than inherited (the builder says so on every
one). Installed: 69 PSSMs. Re-measuring on the 47-exemplar set as `coverage_b5`.

### What the leftover profiles did (coverage_b4 -> coverage_b5)

Same 47 exemplars, same contigs, diffed on (contig, type, feature, start, end,
strand) -- never on the genome id in column 1, which is a per-run random value
and reports every genome as changed.

```
feature   genomes_b4  genomes_b5      calls_b4  calls_b5
PEP46            42          43             42        43
POLY             43          43             43        43
PVP2             47          47             47        47
VP1              46          46             46        46
VP2              47          47             47        47
VP3              36          37             36        37
VP4              38          43             38        43
VP5              25          27             25        27
X                 2           2              2         2
```

**+18 calls, 0 lost, 15 genomes changed at coordinate level.** VP4 is the big
one: 38 -> 43 of 47.

Three of the six clipped calls from the window audit were repaired, and the
repaired lengths are the ones the adjacency test predicted:

```
321302.3   pep46   34 aa -> 41 aa     (predicted ~41)
2583366.4  pVP2   474 aa -> 486 aa, VP4 549 -> 612 nt, VP3 start -42 nt
2920479.4  VP3    212 aa -> 235 aa     (start moved 69 nt upstream)
```

#### A structural check that came out clean

```
VP3 end  -> polyprotein stop     35 / 35   (100.0%)
polyprotein start -> pVP2 start  38 / 43   ( 88.4%)
```

VP3 terminating exactly at the polyprotein stop in every genome that calls both
is the literature layout (pVP2-VP4-VP3, VP3 C-terminal) confirmed independently
of any length window.

#### And one that did not

```
                 b4                      b5
pVP2 -> VP4   23/38 flush (60.5%)     23/43 flush (53.5%)
VP4  -> VP3   21/35 flush (60.0%)     23/36 flush (63.9%)
```

**The flush COUNT did not fall — it is 23 in both runs.** The percentage fell
because five genomes gained a VP4 call they did not have before and all five
are gapped. That is the same trade as the POLY singletons: absent became
present-but-clipped, which is better but is not clean, and saying "VP4 38 -> 43"
without this line would be misleading.

VP4 -> VP3 improved on both measures: two more flush and 486 nt less total gap.

**20 genomes still have a pVP2 -> VP4 gap**, 318 nt in the worst case. These
are cleavage junctions, so `downstream_ext` is correctly 0 and the call is the
profile's match extent -- more profiles alone will not close them. The fix is
the one used on the Arteriviridae ORF1b cassette: project the junction from
reference polyproteins and place it on the protease's own chemistry (VP4 is a
Lon-like serine protease with documented (S/A)-X-A cut sites). Not yet done.

#### The pVP2/VP2 collision: one of two resolved

`2920476.3` is fixed -- a leftover profile, `VP2.3`, now calls VP2 at 351 aa
against pVP2's 395, so the two features are distinct again.

`2920474.4` is NOT fixed. `PVP2.1` and `VP2.1`, both Avibirnavirus-trained,
still call 201..1385 identically. Both are clipped to the same conserved core,
and pVP2 also starts 165 nt after the polyprotein start. Open.

#### New short calls created by the new profiles

```
PEP46   2777470.3 = 29 aa   (literature floor 36)
VP4     663281.3 = 130, 2725373.3 = 152, 1230254.4 = 155, 1230254.7 = 155
PVP2    663281.3 = 479
VP3     2777470.3 = 302 aa  -- inside the literature window (<=310) but
                               over the generated JSON max of 301
```

All of these are in genomes that previously had no call for that feature. They
are gains, and they are clipped. Both halves go in the record.

### VP4 junction chemistry, measured rather than cited

The 20 remaining pVP2 -> VP4 gaps cannot be closed with more profiles. The
Arteriviridae ORF1b cassette had the same problem and was fixed by projecting
the junction from reference polyproteins and then placing it on the protease's
own chemistry. Before reusing that method here it has to be shown that birnavirus
chemistry is sharp enough to place a cut, and **it is only sharp enough at one
of the two junctions.**

Measured on the 46 junctions that are already flush in `coverage_b5` --
Avibirnavirus, Aquabirnavirus, Blosnavirus and Entomobirnavirus:

```
pVP2 | VP4   P2: W=13, L=10        P1: A=22, T=1   P1': S=13, A=10
VP4  | VP3   P2: T9 M5 L3 N3 I1 G1 Q1   P1: A=23   P1': S=14, A=9
```

P1 is Ala in 45 of 46 and P1' is Ser or Ala in 46 of 46. P2 is conserved at the
pVP2|VP4 junction ([WL] in 23 of 23) and carries no signal at VP4|VP3.

**The uniqueness test -- how many motif matches fall inside the window the
alignment localises to:**

```
                       +/-4     +/-6     +/-8    +/-12
A|[SA]  both          71.7%    71.7%       -     47.8%
[WL]A|[SA]  pVP2|VP4  95.7%    95.7%    95.7%    95.7%
A|[SA]      VP4|VP3   47.8%    47.8%    47.8%    43.5%
```

For comparison, the Arteriviridae nsp4 motif E/Q|[GSA] was unique in **98.1%**
of +/-12 windows. **Alanine followed by Ser/Ala is far too common to place a cut
on its own** -- the bare motif is worse than a coin toss at +/-12, and quoting
the published (S/A)-X-A|(A/S/G) consensus without this measurement would have
produced a confidently wrong cassette.

With P2 included the pVP2|VP4 junction is as clean as the arteriviral one:
unique in 22 of 23 even at +/-12. The single exception is `1661392.4`, the one
P1=Thr junction (`DSWT|SGAQ`).

VP4|VP3 stays ambiguous, but **the ambiguity is a near-miss, not a real one.**
The offset of the competing site from the true cut:

```
   +1  : 8   Avibirnavirus   ...YDLAMA | ASEFKE...   true A|A, rival A|S one over
   -3  : 4   Aquabirnavirus  ...LMASNA | SGMDEE...
  -10  : 1   311176.4        ...WFYCGA | ADEETI...
```

Twelve of thirteen competitors are within 3 residues. So the design is:

- **pVP2 | VP4** -- localise by alignment, place by chemistry on [WL]A|[SA].
- **VP4 | VP3** -- localise by alignment, place by chemistry on A|[SA], and
  break ties by distance from the aligned position. Worst-case error ~3
  residues, against gaps of 35-108 codons.

Neither junction is placed on chemistry alone and neither is placed on
alignment alone.

## Projected junctions merged, and three windows corrected (20 Nov 2026)

### What went into the collections

The 33 junctions the projector placed produced 37 collection members. Four
were already present by genome id (`311176.4` in PVP2 and VP3, `2926591.3` in
VP3, `321302.3` in VP4) and were skipped, so **33 members were appended**:

| feature | before | appended | after |
|---|---|---|---|
| `PVP2` | 294 | 13 | 307 |
| `VP3`  | 216 | 11 | 227 |
| `VP4`  | 179 |  9 | 188 |

They carry a header in the same shape as the four members already constructed
on literature coordinates, so provenance is readable without a side file:

```
>fig|2583366.4.projected.9002|Largemouth_bass_birnavirus_LBBV-SS-2015|514-827 projected from the VP4 junction chemistry
```

The lineages they add are the ones with no flush junction at all before:
Entomobirnavirus, porcine birnavirus, Largemouth bass birnavirus, BSNV, TV-1,
Rocky Mountain birnavirus and Rotifer birnavirus.

### Three length windows were wrong, and the collections proved it

**`VP4` max 295 -> 320.** The 295 was drawn from Avibirnavirus and
Aquabirnavirus, which are the short end of the family. Four genomes in two
other genera place longer on junctions that BOTH validated: BSNV `311176.4`
305 aa, TV-1 `321302.3` 302 aa, Largemouth bass `2583366.4` and `2926591.3`
314 aa each. Note that these exceed the ORIGINAL reports' coordinates for BSNV
(234) and TV-1 (212) -- the published numbers are narrower than what the
junctions in those genomes support, so this ceiling is measured, not cited.

**`PVP2` min 480 -> 470.** `2583366.4` and `2926591.3` carry a genuinely short
pVP2 at 474 aa. **Adjacency settles it, not length**: the call abuts VP4 with
no gap, so it is a short lineage rather than a clipped profile.

**`VP1` max 935 -> 1050, and this one resolved twelve standing findings at
once.** The 935 came from BSNV (PMID 12477876) and the two vertebrate genera.
The ICTV Birnaviridae profile gives VP1 as **844-1045 aa across the family**
(PMID 30484762), and Drosophila X virus VP1 is 977 aa, the largest birnavirus
VP1 reported at the time (PMID 12069523). The twelve deposits above the old
ceiling were Entomobirnavirus (Culex Y 998, Port Bolivar 1022, Mosquito X
995/998), Ronavirus (Rotifer birnavirus 1045 -- the ICTV maximum exactly) and
six unclassified birna-like viruses at 974-1010. **Every one falls inside the
ICTV range.** The 1050 is 1045 with a small margin; it is not a per-genus
ceiling, because no VP1 length is published for the birna-like lineages.

Both PMIDs were verified through NCBI eutils before use. A direct fetch of
`pubmed.ncbi.nlm.nih.gov/<pmid>/` returns HTTP 203 and a 249-byte "cookies must
be enabled" page, which confirms nothing; `esummary.fcgi` returns the record.

### The three findings that were signed off rather than fixed

**`VP4` SPLIT, 14x 212-243 aa + 3x 305-314 aa.** Deliberate and expected: the
short mode is Avibirnavirus/Aquabirnavirus/TV-1 and the long mode is BSNV plus
Largemouth bass. Every member of both modes begins at a placed pVP2|VP4
junction and ends at a placed VP4|VP3 junction, so this is genus-level length
variation in one protein, not two proteins sharing a bin.

**`PEP11`, `PEP7A`, `PEP7B` SPLIT.** The pVP2 C-terminal maturation peptides
(IBDV pVP2 442-512 releases pep46, pep7a, pep7b and pep11; PMID 11836417,
15269363). The SPLIT finding is an artefact of the length: at 7-13 aa a single
substitution moves a member into its own cluster, so 11 sequences scatter into
5 groups at any identity floor a longer protein would use. **None is modelled**
-- no profile is built and all three are absent from the module JSON, because a
7- to 13-residue PSSM cannot be made specific. The bins are retained as
literature provenance for the pVP2 cleavage cassette.

After these, the length audit and the bin audit both read **unresolved: 0**.

## copy_num applied to four features, withheld from five (21 Nov 2026)

`copy_num` is what makes a feature essential: `viral_genome_quality.pl`
evaluates only features that carry it, and a genome missing one is Poor. It
also switches on the `min_len`/`max_len` check, so the windows had to be right
first -- they were corrected on 20 Nov and re-audited clean before this step.

The rule applied: **copy_num only where every uncalled genome is demonstrably
an incomplete genome.** One complete genome that fails to call a feature is
enough to withhold it, because the flag would then report a profile gap as a
genome defect.

| feature | called | copy_num | why |
|---|---|---|---|
| pVP2 | 47/47 | **1** | no miss |
| VP2 | 47/47 | **1** | no miss |
| VP1 | 46/47 | **1** | sole miss `2920476.3`, a 5139 nt partial with no segment B |
| POLY | 45/47 | **1** | `2920476.3` partial; in `2777470.3` pVP2 *and* VP1 both start at coordinate 1, so the contigs are 5'-truncated and the polyprotein has no start codon |
| VP3 | 43/47 | omitted | three misses are partials, but `2920474.4` is complete at 6740 nt -- the documented pVP2/VP2 collision |
| VP4 | 43/47 | omitted | `3081665.5` is complete and its VP4 slot measures 265 aa inside the 165-320 window, yet the polyprotein carries **no A\|[AS] site anywhere in residues 630-810** |
| pep46 | 43/47 | omitted | `663281.3` is complete and an A\|FG site sits exactly where pep46 should begin, 46 aa long, but its C-terminal chemistry does not resolve |
| VP5 | 27/47 | omitted | genus-restricted; copy_num would flag every genus lacking it |
| X | 2/47 | omitted | Blosnavirus/Telnavirus only |

### Verified end to end, not assumed

A `copy_num` that is wrong is silent: annotation is unchanged and only quality
scoring moves. So it was measured through the real GTO pipeline
(`New-annotate-viral-taxon.pl` -> `viral_genome_quality.pl`) on five genomes
chosen to span the decision:

| genome | state | verdict |
|---|---|---|
| `10995.1006` | complete IBDV | Good |
| `3081665.5` | complete, VP4 uncalled | **Good** -- confirms withholding copy_num from VP4 |
| `663281.3` | complete, pep46 uncalled | **Good** -- confirms the same for pep46 |
| `2777470.3` | 5'-truncated | **Poor**, `Missing VP2-VP4-VP3 polyprotein` |
| `2920476.3` | 5139 nt partial | **Poor**, `Feature is too short` x2 (pVP2 412, VP2 350) |

Had VP4 and pep46 carried copy_num, the two complete genomes would have been
reported Poor for a profile gap rather than a genome defect.

### Two coverage limits, deliberately not repaired

Both were found by the bracketed-slot test: a mature peptide that nothing
called, in a genome where the neighbouring products bound a gap of the right
size. Both were rejected as projected members because the junction chemistry
does not place them, and a member projected onto an unvalidated junction is the
exact defect this build removed from NSP2 and from the pVP2 deposits.

- **`3081665.5` VP4.** POLY 162-3356, pVP2 ends 1664, VP3 starts 2397, so the
  slot is 265 aa, inside 165-320. A scan of polyprotein residues 630-810 finds
  **no A|[AS] site at all**, so neither end of the product can be placed.
- **`663281.3` pep46.** VP2 ends 1446 and pVP2 ends 1704; an `A|FG` site sits
  two residues past the VP2 call and yields a 46 aa product, but the downstream
  cut is `A|A`, not the `A|[FW]G` maturation chemistry, so the C-terminus is
  unplaced.

A third genome, `2556101.4`, has a 72 aa pep46 slot -- outside the 36-52 window
-- and `10995.6189`/`10995.6191` have 103/100 aa VP3 slots, both far below the
215 aa floor, in genomes whose polyprotein runs to the contig edge. None is
recoverable.

Noted separately: in `2777470.3` pep46 is called at **29 aa**, below the 36-52
window. It is unevaluated because pep46 carries no copy_num, and the genome is
Poor on other grounds, but it is a real short call.

## The VP4→VP3 junction, measured genome by genome (21 Nov 2026)

`check_cleaved_ends.py` over the 47 coverage exemplars: **30 of 42 genomes
carrying both features put VP3's start exactly flush with VP4's end**, and 12
do not, with offsets from −39 to +273 nt. Length alone cannot say which side is
wrong, so each junction window (VP4's last 8 residues + the gap + VP3's first
8) was translated and scanned for VP4's own cut motif, `A|[AS]`. That motif
places the boundary independently of either profile, and it splits the 12 into
three groups that have nothing to do with each other:

**(a) One side is short and the chemistry says which — 4 genomes.**

| genome | gap | junction window | reading |
|---|---|---|---|
| `2583366.4` | +12 | `EEQVDTYA·SFPK·SKDVTT` | VP4 ends on `A|S` — correct. **VP3 starts 4 aa late.** |
| `2926591.3` | +12 | `EEPVDTYA·SFPK·SKDVTT` | same lineage, same 4 aa |
| `48981.6` | +18 | `TPKQTLPHWTAGSA·SMNPFM` | VP3 starts on `A|S` — correct. **VP4 ends 6 aa short.** |
| `2920479.4` | +273 | …`VAELALSA·ATPGKA` | VP3 starts on `A|A` — correct. **VP4 ends 91 aa short.** |

Not repaired. Each is one or two genomes against a profile that is flush
everywhere else (VP3.5 calls 2, VP4.1 calls 12), so re-deriving the collection
on this adjacency — the NSP7B remedy — would rebuild a profile on two members.
The boundaries are recorded here so a later pass can act on them without
re-measuring.

**(b) Both sides are wrong and the motif is ambiguous — 2 genomes.**
`1230254.4` and `1230254.7` carry `A|[AS]` at two positions 3 residues apart
(`…HWTAS·SASVPFSN…`), and neither is adjacent to either call. `2725373.3` is
the same lineage and its downstream site lands one residue off VP3's start.
Nothing here settles the junction; it needs a published BSNV-lineage cleavage
map, which does not exist for these isolates.

**(c) There is no `A|[AS]` site anywhere in the window — 3 genomes.**
`2556101.4`, `3043168.4` and `595652.4` have a 141–234 nt gap containing no
instance of the motif at all. This is the more interesting result: it is not a
profile defect but evidence that **VP4's specificity is not `A|[AS]` outside
the genera where it was established**. The module should not be forced to
produce a flush junction in lineages whose protease has not been characterised.

**(d) The one overlap.** `3065069.4`, −39 nt: VP4 (`…ELEDMF`) runs 39 nt past
VP3's start (`ADQELD`), with no `A|[AS]` at either end. VP4.5 over-reaches into
VP3 in this genome alone. Unresolved, and now explained rather than merely
listed.

Neither VP4 nor VP3 carries `copy_num`, so none of this reaches genome quality
scoring — which is the reason the copy_num decision was made the way it was.

## qc_duplicate_calls over the 308-genome panel: clean

Every nesting the tool found is declared genome architecture — pVP2, VP2, VP4,
VP3, VP5 and pep46 inside POLY at 100% of calls — and there are **no EXACT
duplicates and no undeclared NESTED calls**. No feature is firing inside
another by accident, which is the failure mode that would have earned a feature
copy_num it should not have.

## qc_cross_feature found a VP4 member that was 71 aa of the X peptide (21 Nov 2026)

`qc_cross_feature.py` reported one MISLABEL, twice over (clusters 1 and 3):

```
X   fig|311176.4.mat_peptide   71  collides VP4  100.0%  fig|311176.4.projected.900
```

The projected VP4 member for Blotched snakehead virus was 305 aa, and its
**first 71 residues were byte-identical to the deposited X peptide**,
`fig|311176.4.mat_peptide.6`. In BSNV the polyprotein is pVP2-**X**-VP4-VP3
(PMID 12477876); the projection had been placed on the pVP2|VP4 junction as if
X were not there, so it started one peptide too far upstream.

**Why this one mattered more than its member count suggests.** It sat in a
single-sequence cluster and therefore had its own profile, `VP4.8`, which is
the only profile covering the BSNV lineage. Uncorrected, VP4 would be called
starting at X's start codon in every BSNV-lineage genome — the feature would
have been called on another feature's locus, which is the failure mode 2a
exists to prevent. Nothing else caught it: the length audit had the member
inside the widened 165-320 window, the profile recalled its own training
sequence, and VP4's call rate was 99.4% on the panel.

**Two independent confirmations of the corrected boundary.** Trimming the 71 aa
prefix leaves **234 aa**, which is the published BSNV VP4 size; and X's own
C-terminus reads `...LFKLPCLRPQA|ADLPISLL...`, an `A|A` — VP4's own cleavage
motif, placing the junction from the chemistry without reference to either
profile.

**Verified on the genome, not just in the collection.** After the rebuild and
install, `311176.4` annotates:

| feature | nt | length | junction |
|---|---|---|---|
| pVP2 | 166-1623 | | |
| X | 1624-1836 | 71 aa | flush with pVP2 |
| VP4 | 1837-2538 | 234 aa | flush with X |
| VP3 | 2539-3375 | | flush with VP4 |

Every junction in the pVP2-X-VP4-VP3 cassette is now exactly adjacent, which is
what the published organisation requires and what the module produced for no
BSNV genome before this fix.

VP4 rebuilt: 188 members -> 3 alignments, 5 PSSMs, 6 leftover -> 5 more
profiles (10 total, down from 11; `VP4.10` was retired at install).

## The `copy_num` set re-tested against the panel, and POLY's one miss explained (21 Nov 2026)

The Arteriviridae build produced a sharper test for `copy_num` -- *is there a
complete genome, routed to this module, in which an essential feature is never
called?* -- and Birnaviridae was put through it rather than left on the earlier
reasoning.

Over the **308-genome held-out panel in 13 taxa**, three of the four essential
features (`PVP2`, `VP2`, `VP1`) are called in every genome.  `POLY` has exactly
one miss: **`2777470.6`, Aedes birnavirus**, which calls pVP2, VP2, pep46, VP3
and VP1 but not the polyprotein they are cut from.

That looks like the defect this test exists to find, and it is not.  The
annotator's own trace shows **all fourteen POLY profiles scoring `0` with
`stop-free 0`** -- not a near-miss under the 284-bit floor, but no stop-free
ORF for any of them to match.  The reason is in the deposit: segment A is
3,717 nt and **pVP2 is called from position 1**, so the 5' UTR and the
polyprotein's own start codon are not there.  VP3 runs to 3717, the contig's
last base.

The genome is 5'-truncated on segment A, and `POLY`'s flag is therefore
**correct**: it says this genome does not carry a complete polyprotein, which
is true.  Lowering the floor would not call it and would weaken a floor that
is already below some of its own training set.  **The essential set stands
unchanged: POLY, pVP2, VP2, VP1.**


## Swept for single-sequence contaminants, clean (21 Nov 2026)

The Arteriviridae build turned up a defect class this module had never been
tested for: **one mislabelled deposit becoming its own profile.**
`audit_bins.py`'s foreign-subset test reports a group only at five members or
more, but `build_leftover_pssms.py --min-seqs 1` builds a profile from a single
sequence -- deliberately, and Birnaviridae `X` is exactly why -- so a single
mislabelled deposit gets the same treatment, with no conservation information
to temper it.

`qc/singleton_contaminants.py` asks, for every member of every collection, for
its best hit in every *other* collection, with no count floor.  Over all 12
Birnaviridae collections and 2,170 sequences: **no member matches a foreign
feature** at 90% identity over 85% of both sequences.

**One pair had to be declared rather than detected.**  The sweep separates a
contaminant from a product by subject coverage -- a contaminant is the same
protein as something else, a product is part of something else -- and that
fails when a product is nearly as long as its precursor.  **pVP2 is VP2 plus
the C-terminal structural peptides, 512 aa against 441 in IBDV, so VP2 covers
86% of pVP2** and clears an 85% floor.  All 398 VP2 members were reported on
the first run.  The right answer is that VP2 *is* the matured pVP2, which is
what the module models, so the pair is declared (`VP2:PVP2`) rather than the
threshold moved.
