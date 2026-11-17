# Picobirnaviridae: routing, the reference budget, and a lesson about
# whole-genome references (17 November 2026)

Extracted from the module's BUILD_NOTES.md (full notes in the Box build
archive). Kept here because the finding is not Picobirnaviridae-specific:
**"the 25 references should be whole genomes" is the right default and is not
always achievable, and the way to tell is to measure co-deposition in the dump
before selecting anything.**

In this family only **101 of 5,860 genome names carry both segments**, so a
whole-genome-only reference set is drawn from 1.7% of the dump and routes
*less* of segment 1 than a mixed set does. Reoviridae and Orthomyxoviridae do
not hit this, because their segments are routinely deposited together. The
merge rule still does the real work -- see "The merge rule is what rescues it"
below -- it just cannot be applied to genomes the dump never co-deposited.

## Routing: why segment 1 does not route alone, and why that is not a defect

The step-11 coverage audit ran 3,100 exemplars (95% identity clustering of the
5,920 contigs, >=1,100 nt) against the 67-reference family routing repo. It
reported 1,621 routed / 1,479 unrouted -- 52.3%. Split by segment, which is the
only way that number means anything:

    seg2 (RdRp)    1,376 / 1,643  = 83.7%
    seg1 (capsid)    232 / 1,408  = 16.5%
    unseg             12 /    41  = 29.3%

Feature calls follow exactly: RdRp 1,380, CP 211, ORF1 199.

### What was ruled out

**It is not reference selection, and whole-genome references make it worse.**
Only **101 of 5,860 genome names in this dump carry both segments.** All 101
were built as merged references and 25 greedy-selected against the exemplars:

    reference set                     seg1     seg2    overall
    current mixed 25                 16.4%    83.7%     52.0%
    best 25 whole genomes            10.9%    81.4%     48.5%
    all 101 whole genomes (ceiling)  14.1%    82.5%     50.5%

Every whole genome available routes *less* segment 1 than the mixed set, because
restricting to co-deposited genomes discards the 2,165-contig seg1 pool the
seg1-only references were drawn from. Reoviridae and Orthomyxoviridae do not hit
this: their segments are routinely deposited together.

**It is not a coverable gap.** All-vs-all among the 1,408 seg1 exemplars at the
annotator's own BLASTn settings: **488 hit nothing but themselves** at 150 bits.
The greedy curve does not saturate -- 10 refs 14.7%, 20 refs 20.7%, 40 refs
29.8%, 60 refs 36.3%, still climbing ~5 genomes per reference. Picobirnavirus
segment 1 has no nucleotide-level core to cluster on.

### The merge rule is what rescues it

Measured on the 83 both-segment genome names that have exemplars:

    routed as a MERGED whole genome   70 / 83  = 84.3%
    their seg1 contigs, ALONE         44 / 325 = 13.5%

A capsid segment presented on its own routes 13.5% of the time; the same capsid
delivered attached to its RdRp partner routes 84.3%, and once the genome routes
the capsid is called by PSSM -- protein-level search, far more sensitive than
nucleotide routing. The held-out whole-genome panel agrees: 77/89 = 87%.

**So the 47.7% unrouted figure is an artifact of auditing loose contigs.** It is
the right way to audit and the number is reported honestly, but the module's
real-world unit is the merged genome, and that is where it performs.

### The routing floor, measured and deliberately left alone

`-mcb` (annotate_by_viral_pssm.pl:109) defaults to 150 and is the binding
constraint. Swept against the 67-reference family repo:

    mcb    seg1     seg2    misrouted to another module (of 3,100)
    150   16.4%    83.7%     4
    100   21.0%    86.4%     5
     60   40.6%    89.7%     6
     40   72.3%    90.9%    10
     30   86.9%    92.5%    81   <- breaks

At mcb 40 segment 1 would go 16.4% -> 72.3% and segment 2 would also gain, at a
cost of 10 misroutes in 3,100 (0.3%). The cliff between 40 and 30 is sharp.

**`-mcb` is a global command-line flag with no per-module field** -- the JSON
block carries only `close_genomes` and `segments`, so this family cannot declare
the floor it needs, and lowering it globally would change every other module in
the same run. On 17 Nov 2026 the decision was taken (JD) to **leave the floor at
150** and ship the module as measured. The ~1,180 standalone capsid records that
would become annotatable at mcb 40 are a known, quantified gap, recorded here
and on the coverage-audit page rather than bought with a global change.

### Was the reference selection itself optimal? Yes, within 1.6 points.

An all-vs-all of the 3,100 exemplars (1,867,408 hits >= 150 bits) allows the 25
references to be chosen by greedy set cover instead of by hand. Enforcing the
budget rule -- one file = one genome, at most one contig per segment, merged
where the dump provides a partner:

    reference set            seg1     seg2   overall
    current 25 (hand)       16.4%    83.7%    52.0%
    greedy-optimal 25       17.0%    85.5%    53.6%

The hand-picked set was already within 1.6 points of the best any 25 files can
do, so the references are **not** the constraint and were left unchanged rather
than churning close_genomes.json, the alignments and the held-out validation for
under two points.

**A trap worth recording.** Greedy over genome *names* reports seg1 39.0% and
overall 64.8% -- a large, entirely false gain. "Marmot picobirnavirus" is one
name covering 298 separate contigs in this dump, and treating a name as a
reference file builds a 298-contig mega-reference that the 25-file budget does
not permit. Only one name in the whole dump has more than 3 exemplar contigs,
which is exactly why the error is easy to miss: 2,726 of 2,727 candidates behave
themselves. Candidates must be genomes, not names.

Independently, the optimal set lands on **7 merged multi-segment files out of
25 -- the same count the hand-picked set has.** That number is not a curator's
choice; it is what co-deposition in this dump allows.

### The unbound strings, checked against the member floor

`gen_collapse.py` reports 10 strings unbound, four of which reach the member
floor of 5: `putative ORF1` 63, `hypothetical protein` 37, `putative capsid` 8,
`polyprotein` 7. The standing rule is that an unassigned row reaching the floor
means the module is unfinished, so each was checked by blastp against the three
collected features rather than waved through.

**661 of 740 sequences carrying those strings hit an existing feature**, nearly
all adoptably (>=80% identity over >=60% of the query):

    putative capsid       -> CP     268 adoptable, 11 weak
    putative ORF1         -> ORF1   197 adoptable, 27 weak
    hypothetical protein  -> ORF1   114 adoptable, 23 weak
    hypothetical protein  -> CP       6 adoptable, 10 weak
    polyprotein           -> RDRP     1 adoptable,  3 weak

These are CP and ORF1 deposited under uninformative names. They are not a
missing feature; they are the reason the homology rescue exists.

**79 sequences hit nothing.** Clustered at 30% identity / 0.6 coverage, exactly
one cluster reaches the floor: **five `polyprotein` sequences of 2,486-2,857 aa,
all from genome 3004038, "Army ant associated picorna-like virus 1".**
Picobirnaviruses do not encode a polyprotein. These are picorna-like virus
records misfiled into the family bin -- the same genome whose contigs route to
Sapelovirus in the coverage audit. They are correctly left unbound, and the
defect is in the source taxonomy, not in this module's rules.

No other cluster reaches 3 members. **The floor gate is satisfied.**
