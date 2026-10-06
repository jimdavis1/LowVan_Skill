#!/usr/bin/env python3
"""Record the 6 October leftovers in section 07 of each coverage audit.

Section 07 ("What it cannot do") is where a documented gap gets named, so it is
where these belong. Two kinds go in:

  * what the projection rebuild changed, including where it COST calls -- a
    page that reports only the gains is not an audit;
  * which genera with no module of their own route into this one, measured
    6 Oct on all 404 genomes >=6 kb of the 36 uncovered genera.

Edits the hand-written prose JSON, never the generated audit JSON, so the note
survives the next regeneration.
"""
import json, os, sys

R = os.path.join(os.environ["LOWVAN_KIT"], "reports", "generators")
B = "<b>%s</b>"

#  module -> (extra p07 paragraphs, replacement verdict07 or None)
NOTE = {}

NOTE["pico_ldr_vp4"] = ([
    B % "The projection rebuild this section called for was done on 6 October, "
        "and it closed the gap." +
    " Hunnivirus annotates every product <i>except</i> VP1 on its two reference "
    "genomes and tiles the whole polyprotein, leaving a 234-residue hole between "
    "the end of VP3 and the start of 2A. Deriving VP1 from that hole &mdash; the "
    "references' own coordinates, nothing predicted &mdash; produced three "
    "training sequences, and three were enough: "
    + B % "VP1 went 0/76 to 75/76, 2B 0/76 to 76/76, 3A 0 to 71/76 and the leader 0 to 54/76."
    " Tottorivirus gained 2B and VP1 at 6/6, Mischivirus the leader at 12&rarr;21/42."
    " The module gained 309 calls in total.",

    B % "Every transferred cut had to clear two tests it did not face before." +
    " The polyprotein must be at least 80% identical to the reference it is "
    "copied from, and the cleavage site must be conserved: P1 is measured across "
    "all the references that place a feature, and the feature is projected only "
    "if that residue is conserved and every projected cut lands on it. Nineteen "
    "feature/genus combinations in this module failed that second test and were "
    "left alone rather than filled in &mdash; Bopivirus VP4 at 88%, Erbovirus 3C "
    "at 86%, Anativirus VP1 at 50% on two references.",

    B % "The module also annotates a genus it does not declare." +
    " Of the 36 Picornaviridae genera no module claims, Boosepivirus routes here "
    "on 15 of 15 exemplars and receives "
    + B % "13.9 features per genome" +
    " &mdash; near-complete annotation for a genus with "
    + B % "zero mat_peptide records of its own" +
    ". Malagasivirus and Mosavirus route here too. That is the partitioning rule "
    "paying off: these genomes share the organisation, so profiles trained on "
    "other genera fit them.",
],
    "Both genus-conditional splits work and neither feature ever fires on the "
    "wrong genus. The projection rebuild of 6 October closed the Hunnivirus gap "
    "this page previously reported &mdash; VP1 0/76 to 75/76, and 309 calls "
    "gained across the module &mdash; using only cut coordinates carried from "
    "genomes at least 80% identical across a conserved cleavage site. What "
    "remains unprojected is refused on evidence: nineteen feature/genus "
    "combinations whose P1 residue is not conserved enough across their "
    "references to place the cut with confidence.")

NOTE["pico_ldr_vp0"] = ([
    B % "The 6 October projection rebuild helped Shanbavirus and cost Salivirus." +
    " Shanbavirus gained VPg (21&rarr;25 of 53) and the leader (0&rarr;2); the "
    "module is net +1 call. But "
    + B % "Salivirus 3A fell from 13/13 to 4/13" +
    ", and the cause is worth naming precisely because it is counter-intuitive: "
    "the 3A collection grew from 26 to 57 sequences, every new one at &ge;80% "
    "identity to a reference that annotates the cut, and with the larger set the "
    "main pass clusters differently. The profile that used to cover Salivirus no "
    "longer exists.",

    B % "That is a re-clustering effect, not bad data, and not the length window." +
    " The window was held at 120&ndash;202 by the rule that a window may widen but "
    "never narrow; widening it alone recovers only three of the nine lost calls. "
    "Recovering the rest means splitting 3A by genus or forcing the leftover pass "
    "to cover it. "
    + B % "It should be done before this module is called finished." ,

    B % "Two genera with no module of their own route here." +
    " Gallivirus lands on this module for 9 of 11 exemplars and Ailurivirus for "
    "5 of 5, measured 6 October across all 404 genomes of the 36 unclaimed "
    "genera. Gallivirus gets 7.9 features per genome; Ailurivirus gets 2.8, which "
    "is routing without much annotation &mdash; it reaches the module and "
    "collects little beyond the polyprotein.",
],
    "Four of the five declared genera are called completely or nearly so, and "
    "Shanbavirus remains limited by contributing one sequence per protein to the "
    "training data. The 6 October projection rebuild left the module net +1 but "
    "moved Salivirus 3A from 13/13 to 4/13 through re-clustering, which is an "
    "open regression, not a data limit. Gallivirus and Ailurivirus route here "
    "without being declared.")

NOTE["pico_noldr_vp4"] = ([
    B % "Tremovirus is in the wrong module, and VP4 is the symptom." +
    " VP4 is called on 54/54 Cosavirus and 16/19 Dicipivirus exemplars and on "
    + B % "0 of 10 Tremovirus" +
    ", unchanged by the 6 October projection rebuild. The cause is size, not a "
    "missing profile: Cosavirus VP4 has a median of 68 residues, Dicipivirus 44, "
    "and Tremovirus "
    + B % "20" +
    ". The feature declares 54&ndash;82 with a bit cutoff of 31, both computed "
    "from the pooled median.",

    B % "A second VP4 feature is not the fix." +
    " The two forms are the same protein, so it would need a second annotation "
    "string for <i>Capsid protein VP4 (1A)</i>, and one annotation string means "
    "one feature. The real statement is that a ~20-residue VP4 is a different "
    "genome organisation &mdash; it is the <i>hepatovirus</i> form (Hepatovirus "
    "VP4 median 23, declared 18&ndash;28 in its own module), not the cosavirus "
    "one. By the rule that a module is defined by genome organisation rather "
    "than phylogeny, "
    + B % "Tremovirus does not belong here" +
    ". Moving it is a partitioning change and is left for a curator; Tremovirus "
    "routes here meanwhile and gets its other twelve features.",

    B % "Rosavirus routes here on 48 of 48 exemplars and has no module of its own." +
    " It receives 7.1 features per genome despite contributing "
    + B % "zero mat_peptide records" +
    " to any collection in this family &mdash; it is annotated entirely by "
    "profiles trained on Cosavirus and Dicipivirus. Measured 6 October across "
    "all 404 genomes &ge;6 kb of the 36 unclaimed genera.",
], None)

NOTE["kobuvirus"] = ([
    B % "2A and VP3 were deliberately left out of the 6 October projection "
        "rebuild, because the source annotation contradicts itself." +
    " Counting the residues either side of every annotated start, the same "
    "genomic position carries two incompatible conventions:"
    "<br><code>VP3&nbsp;&nbsp;LRYVAP|QHWKT&nbsp;&nbsp;61"
    "<br>VP3&nbsp;&nbsp;RHVTKQ|HWKTR&nbsp;&nbsp;46"
    "<br>2A&nbsp;&nbsp;&nbsp;FIVKVQ|RPTYV&nbsp;&nbsp;&nbsp;3"
    "<br>2A&nbsp;&nbsp;&nbsp;AQRPTY|VHWAL&nbsp;&nbsp;55</code><br>"
    "VP3 differs by one residue, 2A by five.",

    B % "Pooled, that reads as an unconserved cleavage site &mdash; P1 'P 53%' "
        "for VP3 and 'Q 48%' for 2A &mdash; and the site test refused both." +
    " That is the correct outcome: a projected cut inherits whichever convention "
    "the nearest reference happened to use, so projecting would have mixed both "
    "forms into one collection. The minority form is the right one in each case, "
    "since 3C cleaves <i>after</i> P1, which puts VP3 at H and 2A at R. "
    + B % "Normalising the majority to the motif would make both projectable" +
    " but rewrites boundaries on 116 source-annotated genomes, so it is a "
    "curation decision and was not taken here. Both features keep their full "
    "collections (VP3 116, 2A 96) and are unaffected at runtime.",

    B % "Three genera with no module of their own are annotated by this one." +
    " Ludopivirus routes here on 7 of 7 exemplars at 10.1 features per genome, "
    "Sakobuvirus on 4 of 5 at 9.8, and Rabovirus on 9 of 12. None is declared by "
    "this module; all three were measured on 6 October across the 404 genomes "
    "&ge;6 kb of the 36 unclaimed genera.",
], None)

NOTE["hepatovirus"] = ([
    B % "VPg and 2B were refused by the cleavage-site test on 6 October and are "
        "no longer projected." +
    " For VPg, P1 is E in 79% of the references that place it &mdash; but the "
    "remaining 21% are four different residues (D, T, G and C), so the cut "
    "cannot be placed with confidence on a target that resembles them. For 2B, "
    "P1 is Q in 58%. Both keep their pre-projection collections (9 and 13 "
    "sequences) and both still build two profiles.",

    B % "That cost the module 13 calls of 1,316, and the trade was deliberate." +
    " 2B fell 102&rarr;97, 3A 108&rarr;102, 2A 106&rarr;105 and VPg 99&rarr;98. "
    "Against that, "
    + B % "480 of this module's 1,602 projected sequences &mdash; 30% &mdash; "
          "had been carried across below 80% identity" +
    " and are gone. A 1% loss in calls in exchange for removing a third of the "
    "training data that could not be trusted residue-by-residue is the right "
    "direction for an annotation module.",

    B % "Gruhelivirus and Crahelivirus route here without being declared." +
    " Three of three and one of one exemplars respectively, at 8.3 and 9.0 "
    "features per genome, measured 6 October across the 404 genomes of the 36 "
    "genera no module claims.",
], None)

NOTE["sapelovirus"] = ([
    B % "2A is not projected: its P1 residue is conserved at only 54%." +
    " The 6 October rebuild added 182 sequences to this module at &ge;80% "
    "identity, and held the leader gain at 171 of 209, but 2A was left alone "
    "because the cut cannot be placed. The module lost 4 calls of 2,603 (2B "
    "&minus;5, VPg &minus;2, against 2A +1 and 3A +1) while dropping every "
    "projected sequence that had come across below 80%.",

    B % "VP2's cut is conserved on a lysine, and that is not an error." +
    " The VP0&rarr;VP4+VP2 maturation cleavage is autocatalytic rather than a 3C "
    "product, and in this genus it sits on a K at 99%. An earlier version of the "
    "test required P1 to be Q or E &mdash; the canonical 3C residues &mdash; and "
    "wrongly refused it. "
    + B % "Conservation is the criterion, not the identity of the residue" +
    ".",

    B % "Parabovirus routes here on 12 of 12 exemplars at 9.8 features per "
        "genome, and has no module of its own." +
    " Three Rabovirus exemplars and one Felipivirus land here as well. Measured "
    "6 October across all 404 genomes &ge;6 kb of the 36 unclaimed genera.",
], None)

NOTE["avihepatovirus"] = ([
    B % "The projected boundaries reproduce the published motif independently." +
    " After the 6 October rebuild, 2A1 ends in <code>NPG</code> on "
    + B % "18 of 18" +
    " sequences (most common tail <code>GVEPNPG</code>) and 2A2 begins with "
    "<code>P</code> on 50 of 51 &mdash; the StopGo <code>NPG|P</code> site "
    "documented for DHAV 2A1. 2A3 begins <code>HLPR</code> on 34 of 34, the "
    "motif the split was made at. None of this was fitted; the cuts were carried "
    "from reference genomes and landed on the motif.",

    B % "Aalivirus and Grusopivirus route here without being declared." +
    " Five of five exemplars each, at 9.8 features per genome, measured "
    "6 October across the 404 genomes &ge;6 kb of the 36 genera no module claims.",
], None)

NOTE["cardiovirus"] = ([
    B % "Mupivirus routes here on 8 of 8 exemplars and has no module of its own." +
    " It receives 7.0 features per genome. Measured 6 October across all 404 "
    "genomes &ge;6 kb of the 36 Picornaviridae genera that no module claims; of "
    "those 36, eighteen route into an existing module and collect between 3 and "
    "14 features with no further work.",
], None)

n = 0
for slug, (paras, verdict) in NOTE.items():
    p = os.path.join(R, "%s_prose.json" % slug)
    if not os.path.exists(p):
        print("  MISSING %s" % p); continue
    d = json.load(open(p))
    if slug == "pico_ldr_vp4":
        #  this section described the gap that the rebuild closed; replace it
        d["lede07"] = ("The gap this page reported has been closed, and what is "
                       "still missing is refused on evidence.")
        d["p07"] = paras
    else:
        d["p07"] = list(d.get("p07", [])) + paras
    if verdict: d["verdict07"] = verdict
    d["date"] = "6 October 2026"
    json.dump(d, open(p, "w"), indent=1, ensure_ascii=False)
    print("  %-18s p07 now %d paragraphs" % (slug, len(d["p07"])))
    n += 1
print("  %d prose files updated" % n)
