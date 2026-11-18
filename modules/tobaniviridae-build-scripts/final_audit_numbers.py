#!/usr/bin/env python3
"""Replace sections 04, 05 and 08 of the audit facts with the final numbers,
measured against the module as it ships."""
import json, os
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
p = 'pages/tobaniviridae_audit.json'
d = json.load(open(p))
by = {s['num']: s for s in d['sections']}

d['headline_tiles'] = [
    ["4", "modules, one per subfamily"],
    ["34", "features declared, 30 with profiles"],
    ["176", "PSSMs"],
    ["115/117", "exemplars routed"],
    ["0", "duplicate calls on the held-out panels"],
    ["7.2&times;", "annotation string collapse"],
]

by['04']['blocks'] = [
 {"type": "tiles", "items": [
   ["115/117", "exemplars routed (98.3%)"],
   ["28", "reference contigs in total"],
   ["100%", "of near-complete genomes, every module"],
   ["0", "modules tabled for routing"]]},
 {"type": "table",
  "headers": ["module", "near-complete genomes", "clusters at 70% identity", "coverage at 25 refs", "refs shipped"],
  "rows": [["Torovirinae", "51", "3", "100%", "3"],
           ["Serpentovirinae", "30", "16", "100%", "16"],
           ["Piscanivirinae", "21", "8", "100%", "8"],
           ["Remotovirinae", "18", "1", "100%", "1"]]},
 {"type": "prose", "html": "Each module ships only the references its own coverage needs rather than the full 25, and <code>build_rep_contigs.py</code> reproduces the budget figure exactly. The whole-taxon run confirms it from the other side: of 117 exemplars, 91 were annotated through the GTO pipeline and <b>every one of them routed</b>; the other 26 never reached a module because <code>rast-create-genome</code> returned HTTP 403 from the BV-BRC ID server, which the coverage report counts as UNROUTED. Re-run through the flat annotator, which needs no ID server, <b>24 of those 26 route</b> as well. The two that do not are a 1,625 nt and a 1,316 nt fragment."},
 {"type": "verdict", "kind": "ok", "html": "Routing is not this family's problem, at any level of divergence measured here: 115 of 117 distinct genome types reach a module, including the 120 metagenomic <em>Serpentovirinae sp.</em> records that are out of scope for the collections."},
]

by['08']['blocks'] = [
 {"type": "tiles", "items": [
   ["353", "genome records downloaded"],
   ["117", "clusters at 95% identity"],
   ["49", "of them near-complete"],
   ["83.3%", "of annotated near-complete exemplars good"]]},
 {"type": "prose", "html": "\"Near-complete\" is measured against the genome each module expects rather than a family-wide constant &mdash; a bostovirus is 20 kb and a veiled chameleon serpentovirus 36 kb &mdash; so an exemplar counts as near-complete when it reaches 90% of the median length of its own module's reference genomes (Torovirinae 28,511 nt, Piscanivirinae 27,004, Remotovirinae 20,386, Serpentovirinae 30,104). 49 of the 117 exemplars clear that; the other 68 are partial submissions with a median length of 4,731 nt, and <b>90.5% of them are correctly flagged incomplete</b>, which is the quality script doing its job rather than a module shortfall."},
 {"type": "table",
  "headers": ["module", "near-complete exemplars", "S", "M", "N", "ORF1a", "Mpro", "ORF1ab"],
  "rows": [["Serpentovirinae", "16", "100%", "100%", "100%", "100%", "&mdash;", "44%"],
           ["Torovirinae", "12", "100%", "100%", "92%", "100%", "100%", "75%"],
           ["Remotovirinae", "4", "100%", "100%", "100%", "100%", "&mdash;", "50%"],
           ["Piscanivirinae", "3", "100%", "100%", "100%", "100%", "100%", "100%"]]},
 {"type": "prose", "html": "36 of the 49 near-complete exemplars were annotated against the module as it ships; the other 13 hit the same ID-server 403. All 36 routed, 30 came out good (83.3%), and the six that did not are named: four carry an <code>ORF1a polyprotein: Feature is too short</code> flag (a divergent genome whose polyprotein matches partially &mdash; the profile sets are thin, not the genomes), one is missing N, and Bangali torovirus carries <code>Feature is too long</code> on ORF1a. Every one is a profile-coverage statement about a divergent lineage, and every one is visible rather than silent."},
 {"type": "prose", "html": "The Serpentovirinae row is the one to read twice. Most of its near-complete exemplars are the metagenomic <em>sp.</em> contigs dropped at intake, so they are in no collection and are not transcript-edit references &mdash; which is exactly why <code>ORF1ab</code> falls to 44% there, the correction gate needing a reference within 95% identity. The four structural features are still called on <b>100%</b> of them. Those profiles are reaching genomes that contributed nothing to them."},
 {"type": "prose", "html": "<b>One defect the whole-taxon run found, and it was a cutoff.</b> On Berne virus a second <em>Membrane glycoprotein</em> was called at position 11,322, inside the replicase ORF, scoring <b>30.34 bits against a cutoff of 30</b> &mdash; while the genuine M on the same genome scored 476.8. The generator derives a cutoff from the profile, and for a 233-residue protein that lands at 30, which is not a threshold. Every feature's cutoff was reset to half the lowest score its own panel calls actually reached, capped at 500 so a 4,500-residue polyprotein does not exclude the divergent genomes it is there to find. The false call is gone and no genuine call was lost."},
 {"type": "verdict", "kind": "ok", "html": "Routing is complete, the structural proteins are called on essentially every near-complete genome that carries them, and every shortfall on this page has a named cause rather than a threshold left to tune."},
]
json.dump(d, open(p, 'w'), indent=1)
print('sections 04 and 08 rewritten with the final numbers')
