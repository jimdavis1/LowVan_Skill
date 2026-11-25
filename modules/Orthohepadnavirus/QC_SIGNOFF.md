# Cross-feature QC: why nothing was moved

Both `qc_cross_feature.py` and `audit_contamination.py` report findings on
these modules, and in a family of nested alternative starts **every one of them
is the architecture, not a defect**.  Neither tool takes a `--protect`
argument, so the structural check has to be applied by hand.  It was, and the
answer is the same each time: the collections are clean and nothing is moved.

## qc_cross_feature.py — Avihepadnavirus, 4 MISLABEL rows

    PRECORE cluster 2   304 aa   collides with C        100.0%
    C       cluster 1   262 aa   collides with PRECORE   99.6%
    PRECORE cluster 1   305 aa   collides with C         99.2%
    S       cluster 1   167 aa   collides with L         99.4%

preC is core plus a 29-residue N-terminal extension and L is S plus the preS
domain, so each pair is ~100% identical over the shorter member's whole length
**by construction**.  A tool comparing cluster masters cannot distinguish that
from a mislabel.

## audit_contamination.py — Avihepadnavirus, 13 proposed C -> PRECORE moves

Not applied.  The measurement that settles it:

| collection | n | lengths |
|---|---|---|
| `C` | 88 | **262 aa, every single one** |
| `PRECORE` | 73 | 305 (68), 304 (2), 295 (3) |

DHBV core is 262 aa and precore 305 aa (RefSeq NC_001344.1).  Both collections
hold exactly the right protein at exactly the right length, and **35 of the 88
C members are an exact C-terminal substring of a precore member** — which is
what "precore is core with an N-terminal extension" means.  A C sequence scores
near-perfectly against the precore canonical because it *is* the back of it.

Applying `--write` here would have moved 13 genuine 262 aa core proteins into
the precore collection, where they would be the short members of a bin whose
profile is built at 305 aa — and `fasta-cluster-pssm` discards a cluster whose
consensus is a substring of another's, so the likely outcome is losing them
altogether.

The 14 further `S -> L` margin hits were blocked by the tool's own length
guard, which is the same architecture seen from the surface ORF.

## The general point

For a taxon whose features are nested alternative starts, cross-feature
identity is **not** evidence of mislabelling and cannot be.  What separates a
real mislabel from the architecture is the length class and which terminus the
sequence shares, both of which are checked in `LENGTH_AUDIT.md` and by
`resolve_surface_forms.py`.  Those two came back clean, and they are the ones
with the power to answer the question.

---

# UNASSIGNED_TRACKING: the one row still above the member floor

After `resolve_surface_forms.py` re-attributes the generic surface and core
strings, one Orthohepadnavirus row still clears the floor and one
Avihepadnavirus row does: **the empty annotation string**, 376 and 7 unique
sequences respectively.  These are deposits BV-BRC carries with no product
string at all, so no name rule can reach them.

The gate exists to catch *a protein the taxon encodes and the module does not
call*.  These are not that, and the measurement says so.  Blasting all 376
against the seven feature canonicals:

| best canonical hit | n |
|---|---|
| L | 121 |
| C | 78 |
| X | 69 |
| S | 49 |
| P | 46 |
| M | 4 |
| preC | 2 |
| no hit | 7 |

Every one of them is a protein the module already models, and 269 of 376 clear
the homology-rescue threshold (>=80% identity over >=60% of the query).  The
features they would join hold 4,248 to 11,119 members each, so they add no
diversity to any profile.

**They are deliberately not rescued.**  174 of the 376 are surface forms, and
homology cannot place a surface form: L, M and S are 100% identical where they
overlap, so a rescue would assign them arbitrarily and contaminate exactly the
collections the structural test was written to keep clean.  Leaving 376 unique
sequences unbound costs nothing; rescuing them by homology would cost the three
surface profiles.

---

# audit_bins.py did not complete, and what was done instead

**It does not scale to bins of this size.** The script caps the within-bin
coherence comparison at `--max-members` (default 400) but does **not** cap the
cross-bin query: every member of every collection is blasted against a database
of all the others. For Orthohepadnavirus that is ~60,000 queries against a
~60,000-sequence database. It ran for 44 minutes and was still on the first of
seven collections, projecting to five hours or more on a contended machine.

It was stopped. What it exists to find is a bin that **should not exist** -- a
pile of partials or a second copy of a named protein that the accessory-
discovery pass declared as `UNC1..UNCn`. **These two modules have no such bin.**
All twelve features across both of them are named, literature-anchored proteins
with a verified length reference; `ACCESSORY_CLUSTERS.tsv` is empty in both.
There is no discovered bin for the audit to judge.

What *was* run in its place is stricter on the question that actually bites
here. Every one of the 59,913 Orthohepadnavirus collection members was blasted
against the modal-length canonical of all seven features **as the module then
stood**, and every member whose best hit was a different feature was examined:

| from | to | n | verdict |
|---|---|---|---|
| S | L | 10,855 | architecture -- S is an exact substring of L, ties break arbitrarily |
| M | L | 4,021 | architecture |
| C | preC | 2,255 | architecture |
| S | M | 21 | architecture |
| **S** | **P** | **64** | **real: polymerase fragments** |
| **L** | **P** | **11** | **real** |
| **M** | **P** | **9** | **real** |
| **X** | **P / preC** | **2** | **real** |

The 86 cross-ORF members were removed and the four affected features rebuilt;
they are logged in `Ortho/CROSS_ORF_REMOVED.tsv`.

**That catch is the argument for having done it.** The 64 in `S` were 225 aa --
inside the 210-250 window, so no length rule could see them -- and they are
polymerase: the cluster's consensus begins `GEHHIRIPRTPARVTGGVFLVDKNPHNTKESRL`
and carries `GLSPFLLTQFTSAICS`, the reverse-transcriptase motif. They had formed
their own 35-member cluster and **built one of the two S profiles**, which would
have called `Small surface antigen` on the polymerase locus. The only thing that
surfaced them was `check_met_start.py`: the cluster's modal first residue was G,
not M, in 100% of its members. After removal, S builds one profile instead of
two and every CDS profile begins at its protein's start.

**Recommended kit change, for review rather than applied here:** make
`--max-members` apply to the cross-bin query as well as the self-comparison. Its
help text already describes it as sampling the bin, `audit_contamination.py`
provides the exhaustive per-sequence pass, and the bin-level verdicts
(DUPLICATE / FRAGMENT / INCOHERENT / PARALOGUE) are properties of the bin that
400 representatives estimate well. That one line is the difference between a
five-hour gate and a five-minute one on an over-sequenced taxon.


## The 24 Nov 2026 re-run, over twelve collections

The pass above was made when the module held seven collections.  It was re-run
after `HBwl`, `p22`, `P17_NONA` and `P17_A` were built, over all twelve, and
reported `**Unresolved: 0`.

The re-run needed the nested pairs declared rather than discovered.  Three
findings came back -- S/M, M/L and PRECORE/C -- and none is a defect: they are
the in-frame isoform families this file has described from the start.  The
pairs created that day had been protected and the pre-existing ones had not,
which is why they surfaced only then.  All seventeen now live in
`collections/<Module>/protected_pairs.json`, which ships with the module:

    C/P17_A, C/P17_NONA, C/P22, C/PRECORE, L/M, L/S, L/WHOLEL, M/S,
    M/WHOLEL, P17_A/P17_NONA, P17_A/P22, P17_A/PRECORE, P17_NONA/P22,
    P17_NONA/PRECORE, P22/PRECORE, S/WHOLEL, WHOLEX/X

A protected pair **suppresses a finding**, so it is a curation decision and the
file is the record of it.  Without it the audit is not reproducible: re-run
anywhere else it reports defects that were deliberately waived.  The envelope
family is settled independently by length -- S 210-249, M 255-299, L 360-431,
disjoint -- so no sequence is misfiled by the waiver.

`BIN_AUDIT.md` and `LENGTH_AUDIT.md` ship beside this file for the same reason.
