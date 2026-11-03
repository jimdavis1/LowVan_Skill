#!/usr/bin/env python3
"""Bin the Adenoviridae dump into per-module feature collections.

FOUR MODULES, partitioned on ACCESSORY SET. Every adenovirus genus carries the
same ~17-gene core -- IVa2, pol, pTP, 52K, pIIIa, penton, pVII, pX, pVI, hexon,
protease, DBP, 100K, 33K, 22K, pVIII, fiber -- and what separates them is what
sits outside it:

    Mastadenovirus    E1A, E1B 19K/55K, pIX, pV, and the whole E3 cassette
                      (12.5K, gp19K, CR1-alpha/beta/gamma/delta, RID-alpha,
                      RID-beta, 14.7K, ADP) and E4 (orf1,2,3,4,6,6/7)
    Aviadenovirus     NO E1, NO E3, NO E4, no pV. Instead: GAM-1, dUTPase, a
                      parvovirus Rep-like gene, a lipase, TWO fiber genes
                      (fiber-1 and fiber-2), and ~40 lineage ORFs
    Siadenovirus      smallest genome (24-27 kb), a sialidase, one E3, no pIX
    Barthadenovirus   p32K, LH1/LH2, RH2/RH4/RH5, E4.1/E4.2/E4.3, core
                      protein 1/2  (this is BV-BRC's label for ICTV
                      Atadenovirus: duck/ovine/bovine atadenovirus, snake
                      and lizard adenoviruses)

Four genuinely different peripheral sets on one shared core is the textbook
case for four modules under the partitioning rule, and the vocabulary says so
independently: 1,867 distinct annotation strings in Mastadenovirus against 452,
104 and 136 in the others, with almost no overlap outside the core.

THE ONE THING THIS FILE DOES DIFFERENTLY FROM EVERY EARLIER TAXON.

    A MASS IS NOT A NAME, AND NEITHER IS AN ORF NUMBER -- SO NEITHER ONE
    BINDS. BOTH GO TO HOMOLOGY.

The adenovirus vocabulary is dominated by masses. In Mastadenovirus the top of
the tail reads "6.8 kDa protein", "24.6 kDa protein", "34 kDa protein",
"12.1 kDa glycoprotein", "10.1 kDa protein", "54.9 kDa protein", "16.1 kDa
protein"; in Aviadenovirus it reads ORF0, ORF1, ORF1A/1B/1C, ORF2, ORF8 ...
ORF67. Together those two shapes are thousands of features.

Some of those masses ARE the conventional names of real adenovirus proteins --
E1B 19K and 55K, L1 52/55K, L3 23K protease, L4 22K/33K/100K, E3 gp19K, E3
14.7K, E4 34K. Writing a rule for the rest of the masses by eye is exactly the
trap `references/annotation-triage.md` names: "p22" was four proteins in
Crinivirus. Here "34 kDa protein" sits at a median of 292 aa (E4 34K, 294 aa)
but "33 kDa protein" sits at 228 aa and "23 kDa protein" at 208 aa, and nothing
in the string says which locus any of them came from.

The ORF numbers are worse, because they are positional by construction and the
position is not stable across the genus: FAdV-1 (CELO), FAdV-4 and FAdV-8
number their right-hand ORFs differently, so "ORF19" is not one protein. One of
them is 1,233 aa, which is DNA-polymerase sized.

So the rule set below binds ONLY strings that name a protein functionally or by
a universally-used adenovirus designation. Every mass-named and every
ORF-numbered deposit goes into the pool, where two things can happen to it, in
this order:

    1. HOMOLOGY RESCUE. blastp against the collections the named rules
       already built; adopted at >=80% identity over >=85% of the query.
       This is what recovers the real E3 and E4 proteins hiding behind
       "12.1 kDa glycoprotein", and the Aviadenovirus core proteins hiding
       behind an ORF number.
    2. ACCESSORY CLUSTERING. Whatever still has no home is clustered by
       mmseqs, and every cluster at or above the member floor is declared as
       an uncharacterized lineage-specific protein. Only a true singleton is
       dropped, which is the rule `filter_unchar.py` enforces and which
       Velarivirus shipped without.

The deposited name is used for ONE thing after that: naming the cluster, if a
majority of its members agree and the literature supports the name.

MATURE PEPTIDES ARE NOT BUILT HERE. The curator's scope is every known protein
AND mature peptide, including signal peptides. Adenovirus mature peptides are
adenain cleavage products of six precursors -- pIIIa, pVI, pVII, pVIII, pX(mu)
and pTP -- and NOT ONE of the 21,947 sequences in this dump is a mature
peptide: there are zero mat_peptide features in the reference genomes and zero
in BV-BRC. They cannot be binned from deposits because they were never
deposited. They are built from the precursor collections by applying published
cleavage sites, in build_mature.py, after this runs. Same for the gp19K signal
peptide. See MATURE_PEPTIDES.md.
"""
import collections, json, os, re, shutil, subprocess, sys, tempfile

W = os.path.dirname(os.path.abspath(__file__)) or "."
T = "A"

#  The member floor. Below it a cluster cannot train a profile; at or above
#  it, the protein is real and the module must call it or say why not.
ACC_FLOOR = 5
#  ...and at least this many DISTINCT sequences. One sequence cannot train a
#  profile however many genomes carry it.
ACC_SEQ_FLOOR = 3

GENERA = {
 "Mastadenovirus":  ["Mastadenovirus"],
 "Aviadenovirus":   ["Aviadenovirus"],
 "Siadenovirus":    ["Siadenovirus"],
 "Barthadenovirus": ["Barthadenovirus"],
}

# --------------------------------------------------------------------------
#  THE SHARED CORE
# --------------------------------------------------------------------------
#  Every genus carries these. The alternates below are the real BV-BRC
#  vocabulary, read out of vocab/<genus>.tsv -- not invented.
#
#  NOTE on the precursor names. pIIIa/pVI/pVII/pVIII/pX are PRECURSORS and
#  IIIa/VI/VII/VIII/mu are the adenain products. BV-BRC uses the two
#  interchangeably on what is always the full-length precursor translation
#  (every "protein VIII" deposit is 227 aa, the same as every "pVIII"), so
#  both forms bind to the precursor here. The mature forms are built
#  separately and are not in this dump at all.

HEXON  = r"^(l3 ?hexon|hexon|hexon protein|hexon capsid protein|major capsid protein|" \
         r"capsid protein ii|protein ii|ii|hexon protein ii|major core protein hexon|" \
         r"l3 hexon protein|hexon \(ii\))$"
PENTON = r"^(penton|penton base|penton base protein|penton protein|penton capsid protein|" \
         r"capsid protein iii|protein iii|iii|penton base protein iii|l2 penton base|" \
         r"penton-base protein)$"
FIBER  = r"^(fiber|fibre|fiber protein|fibre protein|capsid protein iv|protein iv|iv|" \
         r"l5 fiber|fiber capsid protein|fiber head protein|fibre head)$"
PROT   = r"^(protease|proteinase|l3 protease|adenain|endopeptidase|" \
         r"adenovirus protease|cysteine protease|23 ?kda protein[ae]se|" \
         r"l3 23k protease|protease \(23k\)|adenoviral protease|23k protease)$"
POL    = r"^(dna polymerase|pol|adpol|dna pol|e2b dna polymerase|dna pol protein|" \
         r"dna-directed dna polymerase|dna polymerase protein|polymerase|" \
         r"e2b polymerase|dna polymerase e2b)$"
PTP    = r"^(ptp|ptp protein|pre-?terminal protein|preterminal protein|" \
         r"terminal protein precursor|terminal protein precursor ptp|" \
         r"dna terminal protein|dna terminal protein precursor|terminal protein|" \
         r"e2b preterminal protein|e2b pre-?terminal protein|" \
         r"precursor terminal protein|ptp \(pre-?terminal protein\))$"
DBP    = r"^(dbp|dbp protein|dna binding protein|dna-binding protein|" \
         r"single-?stranded dna-?binding protein|ssdna-?binding protein|" \
         r"e2a dbp|e2a dna-?binding protein|early e2a dna-?binding protein|" \
         r"e1a dna binding protein dbp|72k dna binding protein|" \
         r"dna-?binding protein dbp|early protein e2a)$"
IVA2   = r"^(iva2|iva2 protein|encapsidation protein iva2|maturation protein iva2|" \
         r"maturation protein|protein iva2|iva-?2)$"
K52    = r"^(52k|52 ?kda protein|l1 52k|late l1 52 ?kda protein|p52k|52k protein|" \
         r"52/55 ?kda protein|52-55k|52/55 ?k( protein)?|encapsidation protein 52k|" \
         r"l1 52/55k|52k encapsidation protein|52/55 kda protein)$"
K100   = r"^(100k|100 ?kda protein|hexon assembly protein 100k|l4 100k|" \
         r"100 ?kda hexon assembly-?associated protein|hexon-?assembly protein|" \
         r"hexon assembly-?associated protein|hexon-?associated protein|" \
         r"100 ?kda hexon-?assembly associated protein|100k protein|" \
         r"l4 100 ?kda protein|100 kda protein hexon assembly)$"
K33    = r"^(33k|33 ?kda protein|33kda protein|protein 33k|l4 33k|l4 33 ?kda protein|" \
         r"33k protein|late 33 ?kda protein|late 33 ?kda protein fragment)$"
K22    = r"^(22k|22 ?kda protein|protein 22k|l4 22k|l4 22 ?kda encapsidation protein|" \
         r"22k protein|l4 encapsidation protein|encapsidation protein 22k)$"
PIIIA  = r"^(piiia|iiia|piiia protein|iiia protein|protein iiia|capsid protein piiia|" \
         r"capsid protein precursor piiia|capsid protein iiia|protein iiia precursor|" \
         r"iiia capsid protein|l1 piiia|l1 peripentonal hexon associated protein iiia|" \
         r"peripentonal hexon-?associated protein|piiia capsid protein|" \
         r"minor capsid protein piiia|hexon-?associated precursor piiia)$"
PVI    = r"^(pvi|vi|pvi protein|protein vi|capsid protein pvi|" \
         r"capsid protein precursor pvi|protein vi precursor|vi capsid protein|" \
         r"pre-?protein vi|minor capsid protein precursor|minor capsid protein pvi|" \
         r"l3 pvi|minor capsid protein vi|capsid protein vi)$"
PVII   = r"^(pvii|vii|pvii protein|protein vii|l2 ?pvii|core protein pvii|" \
         r"core protein precursor pvii|major core protein precursor pvii|" \
         r"major core protein|core protein vii|l2 pvii protein|" \
         r"major core protein vii|core protein precursor vii)$"
PVIII  = r"^(pviii|viii|pviii protein|protein viii|capsid protein pviii|" \
         r"capsid protein precursor pviii|l4 ?pviii|hexon-?associated protein viii|" \
         r"minor capsid protein pviii|capsid protein viii|protein viii precursor)$"
PX     = r"^(px|mu|px protein|protein x|l2 ?px|core protein px|" \
         r"core protein precursor px|core protein l2 px|l2 ?pmu|pmu|" \
         r"core protein mu|mu core protein|core protein precursor x|protein mu)$"

CORE = {"HEXON": HEXON, "PENTON": PENTON, "FIBER": FIBER, "PROTEASE": PROT,
        "POL": POL, "PTP": PTP, "DBP": DBP, "IVA2": IVA2, "K52": K52,
        "K100": K100, "K33": K33, "K22": K22, "PIIIA": PIIIA, "PVI": PVI,
        "PVII": PVII, "PVIII": PVIII, "PX": PX}

#  U exon protein. Described in human adenovirus as a late l-strand product
#  that overlaps E3 antisense (PMID_claude_generated 17881437), and deposited
#  across Mastadenovirus, Aviadenovirus, Siadenovirus and Barthadenovirus
#  alike, so it is core rather than genus-specific.
UXP = r"^(u exon|u exon protein|u-?exon|u-?exon protein|u protein|uxp|" \
      r"protein u|u exon protein uxp)$"

# --------------------------------------------------------------------------
#  PER-MODULE PERIPHERALS
# --------------------------------------------------------------------------
MAST = dict(CORE, **{
 "UXP": UXP,
 "PV":   r"^(pv|v|pv protein|protein v|core protein v|minor core protein v|"
         r"protein v precursor|v core protein|l2 ?pv|core protein pv)$",
 "PIX":  r"^(pix|ix|pix protein|ix protein|protein ix|capsid protein ix|"
         r"l1 ?pix|minor capsid protein ix|capsid protein pix)$",
 "E1A":  r"^(e1a|e1a protein|control protein e1a|e1a control protein|"
         r"early e1a|early protein e1a|e1a 13s|e1a 12s|e1a transforming protein)$",
 #  E1A 243R and 289R are the 12S and 13S SPLICE products, not separate
 #  genes, and they are N-terminally and C-terminally identical to each
 #  other -- the 12S form is the 13S form minus an internal 46-residue
 #  block. They are therefore declared as splice features, not as two
 #  PSSMs that would fire on one locus. Bound here so the deposits are
 #  kept, then handed to the splice build.
 "E1A_243R": r"^(control protein e1a 243r|e1a 243r|e1a 12s product|"
             r"e1a 243 ?r protein)$",
 "E1B19K": r"^(e1b 19k|e1b 19 ?kda protein|e1b 19k protein|control protein e1b 19k|"
           r"e1b small t antigen|e1b small t-?antigen|19 ?kda e1b protein|"
           r"e1b-?19k|e1b 19 ?kda)$",
 "E1B55K": r"^(e1b 55k|e1b 55 ?kda protein|e1b 55k protein|control protein e1b 55k|"
           r"e1b large t antigen|e1b large t-?antigen|e1b-?55k|"
           r"e1b protein, large t-?antigen|55 ?kda e1b protein|e1b 55 ?kda)$",
 "E3_12_5K": r"^(e3 12\.5k|control protein e3 12\.5k|early e3 12\.5 ?kda glycoprotein|"
             r"e3 12\.5 ?kda protein|e3 12\.5 ?kda immune modulating protein|"
             r"e3 12\.5 ?kda glycoprotein)$",
 #  gp19K. Named for its mass but universally so, and the name carries its
 #  function: it is the MHC class I antigen-binding ER-resident glycoprotein
 #  (PMID_claude_generated 8350398), which is why "19.3 kDa MHC class I
 #  antigen-binding" binds here and a bare "19 kDa protein" does not.
 "E3_GP19K": r"^(e3 gp19k|gp19k|membrane glycoprotein e3 gp19k|e3 gp19 ?kda|"
             r"19\.?[0-9]? ?kda mhc class i antigen-?binding( glycoprotein| protein)?|"
             r"e3 19 ?kda glycoprotein|mhc class i antigen-?binding glycoprotein)$",
 "E3_CR1A": r"^(e3 cr1-?alpha|membrane glycoprotein e3 cr1-?alpha|"
            r"e3 cr1-?alpha protein|cr1-?alpha)$",
 "E3_CR1B": r"^(e3 cr1-?beta|membrane glycoprotein e3 cr1-?beta|"
            r"e3 cr1-?beta protein|cr1-?beta)$",
 "E3_CR1G": r"^(e3 cr1-?gamma[0-9]?|membrane glycoprotein e3 cr1-?gamma[0-9]?|"
            r"e3 cr1-?gamma protein|cr1-?gamma[0-9]?)$",
 "E3_CR1D": r"^(e3 cr1-?delta|membrane glycoprotein e3 cr1-?delta|"
            r"e3 cr1-?delta protein|cr1-?delta)$",
 "E3_RIDA": r"^(e3 rid-?alpha|membrane protein e3 rid-?alpha|rid-?alpha|"
            r"e3 rid-?alpha protein|e3 10\.4k|10\.4 ?kda protein rid-?alpha)$",
 "E3_RIDB": r"^(e3 rid-?beta|membrane protein e3 rid-?beta|rid-?beta|"
            r"e3 rid-?beta protein|e3 14\.5k|14\.5 ?kda protein rid-?beta)$",
 "E3_14_7K": r"^(e3 14\.7k|control protein e3 14\.7k|e3 14\.7 ?kda protein|"
             r"14\.7 ?kda protein e3)$",
 #  ADP, the adenovirus death protein, is E3 11.6K and is deposited as
 #  "adenovirus death protein" and as "E3 11.6K".
 "E3_ADP": r"^(adp|adenovirus death protein|e3 11\.6k|e3 adp|"
           r"adenovirus death protein adp|11\.6 ?kda protein adp)$",
 "E4ORF1": r"^(e4 ?orf ?1|e4 orf1 protein|control protein e4orf1|"
           r"e4 control protein orf1|early 4 orf1 protein|e4 hypothetical protein 1)$",
 "E4ORF2": r"^(e4 ?orf ?2|e4 orf2 protein|control protein e4orf2|"
           r"e4 control protein orf2|early 4 orf2 protein|e4 hypothetical protein 2)$",
 "E4ORF3": r"^(e4 ?orf ?3|e4 orf3 protein|control protein e4orf3|"
           r"e4 control protein orf3|early 4 orf3 protein|e4 hypothetical protein 3)$",
 "E4ORF4": r"^(e4 ?orf ?4|e4 orf4 protein|control protein e4orf4|"
           r"e4 control protein orf4|early 4 orf4 protein|e4 hypothetical protein 4)$",
 "E4ORF6": r"^(e4 ?orf ?6|e4 orf6 protein|control protein e4orf6|e4 34k|"
           r"e4 control protein orf6|early 4 orf6 protein|control protein e4 34k|"
           r"e4 control protein 34 ?kda|e4 34k protein|34k protein e4)$",
 "E4ORF67": r"^(e4 ?orf ?6/7|e4 orf6/7 protein|control protein e4orf6/7|"
            r"e4 control protein orf6/7|early 4 orf6/7 control protein|"
            r"e4 34k/orf6-7)$",
 #  The L1 i-leader protein and the L1 13.6K are the same locus in some
 #  species and distinct in others; both names are deposited and both are
 #  real. Kept apart; consolidation decides if they are one cluster.
 "I_LEADER": r"^(i-?leader|i leader|l1 i-?leader|leader i)$",
 "L1_13_6K": r"^(protein 13\.6k|l1 13\.6 ?kda protein|13\.6 ?kda protein|l1 13\.6k)$",
 "AGNO": r"^(agnoprotein|15\.3 ?kda agnoprotein|agno protein)$",
})

AVI = dict(CORE, **{
 "UXP": UXP,
 #  TWO fiber genes. FAdV-1 (CELO) and FAdV-4 carry a long and a short
 #  fiber at the same locus, 710 aa and 410 aa in CELO, and they are not
 #  alternative starts of one another -- they are separate genes with
 #  separate promoters. FIBER above would pool them into one bimodal
 #  collection, which the engine forbids, so the explicit names bind first.
 "FIBER1": r"^(fiber-?1|fiber ?1|fiber-?1 protein|long fiber protein|"
           r"long fiber|fibre-?1|fiber 1 protein)$",
 "FIBER2": r"^(fiber-?2|fiber ?2|fiber-?2 protein|short fiber protein|"
           r"short fiber|fibre-?2|fiber 2 protein)$",
 #  GAM-1, the CELO anti-apoptotic protein (PMID_claude_generated 9060680).
 "GAM1": r"^(gam-?1|gam-?1 protein|gam1|anti-?apoptotic protein gam-?1)$",
 "DUTPASE": r"^(dutpase|deoxyuridine triphosphatase|dutp pyrophosphatase|"
            r"deoxyuridine 5'-?triphosphate nucleotidohydrolase|putative dutpase)$",
 "REPLIKE": r"^(rep-?like protein|rep protein|parvovirus rep-?like protein|"
            r"putative rep protein|rep-?like)$",
 "LIPASE": r"^(lipase|putative lipase|lipase-?like protein)$",
})

SIA = dict(CORE, **{
 "UXP": UXP,
 "SIALIDASE": r"^(sialidase|neuraminidase|sialidase protein|putative sialidase)$",
 "E3": r"^(e3|e3 protein|siadenovirus e3)$",
 "HYD": r"^(hyd|hyd protein|hydrophobic protein)$",
 #  SHUTOFF WAS A FEATURE HERE AND IT WAS WRONG. The 8 "shutoff protein"
 #  deposits are 703 aa and consolidation put every one of them into K100
 #  -- correctly. The late 100K protein IS the host shut-off protein; the
 #  string is a synonym for a feature already declared, not a gene of its
 #  own, and declaring it twice would have had two profiles firing on one
 #  locus. It is listed under K100 below instead.
})
assert SIA["K100"].endswith(")$"), SIA["K100"][-10:]
SIA["K100"] = (SIA["K100"][:-2] + r"|shutoff protein|host shutoff protein|"
               r"shut-?off protein)$")

BARTHA = dict(CORE, **{
 "UXP": UXP,
 "P32K": r"^(p32k|p32k protein|32 ?kda protein p32k)$",
 "LH1": r"^(lh1|lh1 protein|left-?hand orf 1)$",
 "LH2": r"^(lh2|lh2 protein|left-?hand orf 2)$",
 "LH3": r"^(lh3|lh3 protein|left-?hand orf 3)$",
 #  RH1 and RH6 were MISSING from the feature set entirely -- 7 RH1 deposits
 #  went to UNASSIGNED and the rest were swept into one accessory cluster
 #  that mixed RH6, RH1 and RH0, three different genes in one would-be
 #  profile.
 #
 #  RH0 AND RH3 ARE NOT DECLARED, on the primary literature
 #  (PMID_claude_generated 12504581, Both 2002, the right-end transcription
 #  map of OAdV-7):
 #    RH0  exists in BAdV-4 and is related to the F-box family, but the paper
 #         says it was overlooked precisely because "it lacks an initiating
 #         methionine". An ORF with no start codon is not an independent
 #         protein.
 #    RH3  "the short ORF RH3 is not conserved in BAdV-4" and "No ORF
 #         RH3-specific transcript was identified". The paper set out to ask
 #         "whether each of the ORFs could actually code for a protein" and
 #         its answer lists the products as RH1, 2, 4, 5 and 6. RH3 is not
 #         among them.
 #  That both failed to build was the data agreeing with the paper, not a
 #  pipeline problem. Curator's call, 3 Nov 2026: get rid of them.
 "RH1": r"^(rh1|rh1 protein|right-?hand orf 1)$",
 "RH6": r"^(rh6|rh6 protein|right-?hand orf 6)$",
 "RH2": r"^(rh2|rh2 protein|right-?hand orf 2)$",
 "RH4": r"^(rh4|rh4 protein|right-?hand orf 4)$",
 "RH5": r"^(rh5|rh5 protein|right-?hand orf 5)$",
 "E4_1": r"^(e4\.1|e4\.1 protein)$",
 "E4_2": r"^(e4\.2|e4\.2 protein)$",
 "E4_3": r"^(e4\.3|e4\.3 protein)$",
 #  "core protein 1" and "core protein 2" are the atadenovirus deposits'
 #  names for what are pVII and pX elsewhere. NOT bound to those keys on
 #  the name -- homology decides, same as everything else here.
})

#  Mastadenovirus carries TWO fiber genes in the enteric species. HAdV-F40
#  and F41 encode a long and a short fiber from adjacent genes; they are
#  paralogues, not alternative starts, and they assemble onto different
#  penton vertices. Pass 1 left them in the accessory pool as UNC2 (median
#  555 aa, labelled "long fiber"/"fiber 1") and UNC5 (median 368 aa, "fiber
#  2"/"short fiber"), which is the clustering telling us the module was
#  missing two declared features.
MAST["FIBER1"] = AVI["FIBER1"]
MAST["FIBER2"] = AVI["FIBER2"]

#  E4 ORF5 / "E4 17 kDa protein". Pass 1 found it as UNC1, 59 sequences at a
#  median of 151 aa with 19 "Orf5" and 15 "E4 17 kDa protein" labels -- a
#  real E4 product the first rule set did not declare.
MAST["E4ORF5"] = r"^(e4 ?orf ?5|e4 orf5|e4 17 ?kda|e4 17 ?kda protein|orf5)$"

PSSM = {"Mastadenovirus": MAST, "Aviadenovirus": AVI,
        "Siadenovirus": SIA, "Barthadenovirus": BARTHA}

#  ---- PARALOGUE SETS: homologous features consolidation must not merge ----
#  consolidate() moves a sequence to whichever feature holds the majority of
#  its homology cluster, which is right for a mislabelled deposit and wrong
#  for two genuinely distinct genes that are homologous to each other.
#
#  The long and short fibers are ~35-50% identical over the shaft, so they
#  cluster together at --min-seq-id 0.3 and the majority rule would empty
#  one of them. The E1A 289R and 243R forms are worse: 243R IS 289R minus an
#  internal 46-residue block, so they are ~100% identical outside it and no
#  homology measure can ever separate them -- exactly the VP1/VP2 case the
#  parvo build had to carve out. A move BETWEEN two members of a set is
#  refused; a move INTO or OUT OF the set from elsewhere is still allowed,
#  which is what catches a foreign protein sitting in a fiber collection.
#  ---- NAMING A CLUSTER THE LITERATURE HAS ALREADY DESCRIBED ----
#  An accessory cluster is a homology group the build discovered. UNCn says
#  only "group number n", and a module full of UNCs is a module that calls
#  proteins without saying what they are.
#
#  Where the literature describes a gene AND the cluster's own dominant
#  deposited label names it, the cluster is named for that gene instead of
#  UNCn. This does NOT reintroduce the U-number trap: the trap is binding
#  sequences together BECAUSE they share an ORF number. Here the clustering
#  already decided the group by homology, independently of the label, and
#  the label is only being read off afterwards to name what was found. The
#  evidence for each name is in the module's LITERATURE_INVENTORY.md.
#
#  A name is listed here only when a paper describes the gene. Where the
#  literature is silent the number is KEPT -- a functional string the
#  literature does not support is worse than a position label.
CLUSTER_NAME = {
 "Siadenovirus": {
   #  Two of the FIVE genus-defining genes were sitting in UNC1 and UNC2.
   #  TAdV-3 NC_001958: ORF7 24509..25165 = 218 aa, ORF8 complement(
   #  25201..25698) = 165 aa; the clusters' median lengths are 225 and 174.
   "orf7": "ORF7",
   "orf8": "ORF8",
 },
 "Aviadenovirus": {
   #  CELO reannotation, PMID_claude_generated 14604445, Table 1.
   "orf2":  "REPLIKE",      # parvovirus NS1/Rep homologue, inactive helicase
   "orf12": "ORF12",        # same family as ORF2/ORF13, function unresolved
   "orf13": "ORF13",
   "orf16": "ART",          # mono-ADP-ribosyltransferase, Pfam PF01129
   "orf22": "ORF22",        # binds pRb and activates E2F
   "orf9":  "ORF9", "orf10": "ORF10", "orf11": "ORF11",
   "orf14": "ORF14", "orf17": "ORF17", "orf20": "ORF20",
 },
}

#  Lengths of each gene as the reference genome's DATABASE RECORD annotates
#  it. These are NOT literature boundaries and must never be used as one: a
#  GenBank /product= or join() is a submitter's claim with at most a curation
#  pass over it, and carries the same kind of error BV-BRC does. NC_004037
#  calls the LH3 locus "E1B 55K" and atadenoviruses have no E1B at all.
#
#  What they ARE good for is a CONSISTENCY check that a collection cannot
#  perform on itself: a profile trained on two mixed genes has a perfectly
#  self-consistent length distribution, and only an outside number exposes
#  it. Used that way by the paralogue audit, which additionally requires
#  homology to agree before it moves anything, and by check_coords.py.
#
#  Where the literature establishes a boundary, the literature wins and the
#  disagreement is written down. See each module's LITERATURE_INVENTORY.md.
REFLEN = {}
if os.path.exists(os.path.join(W, "dbcoords.json")):
    for _m, _f in json.load(open(os.path.join(W, "dbcoords.json"))).items():
        if isinstance(_f, dict):
            REFLEN[_m] = {k: v for k, v in _f.items()
                          if not k.startswith("_") and isinstance(v, int)}

# --------------------------------------------------------------------------
#  PRECURSOR vs MATURE PEPTIDE: they must not share a collection
# --------------------------------------------------------------------------
#  A cleaved protein is deposited BOTH ways. Barthadenovirus genome 114430
#  deposits p32K at 286 aa (precursor) and at 274 aa (mature) and both landed
#  in P32K, so one profile would have been trained on two forms. That is the
#  Parvoviridae VP1/VP2 defect exactly, and there it caused the real VP2
#  clusters to be dropped as substrings of VP1 and VP2 to ship as a VP1
#  profile, called at VP1's coordinates in 281 of 295 genomes.
#
#  THE TEST IS STRUCTURAL, NOT NOMINAL, and it needs two signals:
#    1. the shorter sequence is an EXACT C-terminal substring of a longer
#       deposit FROM THE SAME GENOME -- a cleaved product always is, and
#       nothing else reliably is; and
#    2. the shorter sequence does NOT start with Met -- a cleaved N-terminus
#       cannot, because the cut is made by the proteinase, not the ribosome.
#
#  Both are required. (1) alone also matches a submitter annotating a
#  downstream start codon, which still begins with Met. (2) alone matches a
#  partial deposit.
#
#  Only features whose cleavage is ESTABLISHED IN THE LITERATURE are listed.
#  Barthadenovirus p32K: measured by N-terminal sequencing of the mature
#  virion protein, AYRKRKKKFQRPKPR, cleaved at the type I site LAGG|AYR
#  (PMID_claude_generated 8659111, Vrati 1996, Table 1). Nothing is split on
#  a site this project inferred.
MATURE = {
 "Barthadenovirus": {"P32K": "P32K_MAT"},
}

PARALOG = {
 "Mastadenovirus":  [{"FIBER", "FIBER1", "FIBER2"}, {"E1A", "E1A_243R"}],
 "Aviadenovirus":   [{"FIBER", "FIBER1", "FIBER2"}],
 "Siadenovirus":    [],
 #  The atadenovirus left-hand and right-hand genes are FAMILIES of related
 #  paralogues, not one gene under several names. Without this, consolidation
 #  measured RH3 and RH2 as homologous to RH4 and emptied both into it --
 #  RH3 shipped with ZERO sequences while 4 real RH3 deposits sat inside RH4,
 #  and RH2 lost 7 of its 22. The same risk applies to LH1/LH2/LH3 and to
 #  E4.1/E4.2/E4.3, which are likewise neighbouring paralogues.
 "Barthadenovirus": [{"RH1", "RH2", "RH4", "RH5", "RH6"},
                     {"LH1", "LH2", "LH3"},
                     {"E4_1", "E4_2", "E4_3"}],
}

# --------------------------------------------------------------------------
#  THE POOL: masses and ORF numbers. Neither binds.
# --------------------------------------------------------------------------
#  A string matching POOL is NOT unassigned-and-forgotten; it is routed to
#  homology rescue and then to accessory clustering. It is kept out of the
#  name rules because the string carries no homology claim.
POOL = re.compile(
    r"^("
    r"orf ?[0-9]+[a-z]?( protein)?|"              # ORF19, ORF14A, ORF1b ...
    r"orf ?[0-9]+[a-z]?/[0-9a-z]+|"               # ORF1a/b
    r"orf ?[0-9]+ [a-z0-9. ]+|"                   # "ORF13 protease", "ORF11 pVI"
    r"[0-9]+(\.[0-9]+)? ?kda? ?(protein|glycoprotein|polypeptide)?|"
    r"[0-9]+(\.[0-9]+)?k( protein)?|"
    r"e3 [0-9]+(\.[0-9]+)?k( protein)?|"          # E3 10.1K, E3 19.4K, E3 31.6K
    r"e3 [0-9]+(\.[0-9]+)? ?kda protein|"
    r"e1a [0-9]+(\.[0-9]+)? ?kda protein|"
    r"putative protein [a-z]-?[0-9]+|"            # "putative protein A-106"
    r"hypothetical protein( [0-9]+)?|"
    r"predicted protein|unknown protein|unknown|"
    r"unnamed protein product|putative protein|undefined (cds|mat_peptide)|"
    r"conserved hypothetical protein|uncharacterized protein|"
    #  A BARE REGION LABEL is not a fusion and must not be discarded. "E3",
    #  "E4", "L2" and "early region 3" name the transcription unit a real
    #  protein came out of, and pass 3 sent 909 Mastadenovirus features to
    #  NOT_MODELLED on those strings alone. Homology can place them; the
    #  string cannot. Same for "splicing factor", which is 302 features and
    #  is L4 33K's ACTIVITY, not a locus of its own -- rescue puts it in 33K
    #  if it really is 33K and nowhere if it is not.
    r"e[1-4][ab]?[0-9]?|l[1-5]|"
    r"(early|late) region [0-9][ab]?|"
    r"splicing factor|mrna splicing factor"
    r")$", re.I)

#  Fusions and multi-protein labels: recorded, never binned. A fusion string
#  names two loci and cannot train either profile.
NOT_MODELLED = re.compile(
    r"^(e1a/e1b|e1b/e3|pvii/px|hexon/penton|"
    r"protein sequence is in conflict.*|"
    r"splicing factor|"                 # names an activity of 33K, not a locus
    r"early region [0-9]|late region [0-9]|"
    r"l[1-5]|e[1-4][ab]?|"              # a bare region label, not a protein
    r"polyprotein)$", re.I)

#  The source saying the sequence is not whole. Discarded before any rule is
#  tried -- no length window substitutes for this.
DISC = re.compile(r"truncated|,\s*partial|\bpartial\b|fragment|incomplete|"
                  r",\s*(n|c)-terminal", re.I)


#  ---- string normalisation, applied before any rule is tried --------------
#  Pass 1 declared E3_14_7K and then left 79 features of "E3 14.7K protein"
#  and "E3 14.7 protein" in the accessory pool, where they formed UNC7. The
#  rule said `e3 14\.7k` and the deposit said the same thing with a suffix.
#  Enumerating "( protein)?" on 60 alternations invites exactly that class of
#  miss, so the suffixes and hedges come off the STRING once, here, instead.
#
#  Nothing semantic is removed. "putative" is a claim about identity, not
#  about completeness, and stripping it only lets the rule match -- the
#  deposit still has to clear the same length window and the same homology
#  checks as any other member. Partial/truncated strings are discarded by
#  DISC long before this runs.
#  "hypothetical" is deliberately NOT in this list. "hypothetical protein" is
#  the canonical unnamed string and must survive intact to reach POOL -- it
#  is 1,115 Aviadenovirus features on its own. Stripping the hedge would
#  leave the bare word "protein" and POOL would stop recognising it.
_HEDGE = re.compile(r"^(putative|probable|predicted|possible)\s+(?=\S)", re.I)
_SUFF = re.compile(r"\s+(protein|polypeptide|gene product|cds)$", re.I)


def normalise(s):
    t = " ".join(s.split()).lower()
    t = t.replace("_", " ").replace("‐", "-").replace("–", "-")
    prev = None
    while t != prev:
        prev = t
        t = _HEDGE.sub("", t)
        t = _SUFF.sub("", t)
    return t.strip()


def expand(rx):
    """Add the NORMALISED form of every literal alternative to a rule.

    The rules are written the way the deposits read -- "maturation protein",
    "shutoff protein", "adenovirus death protein" -- but normalise() takes
    the suffix off before the rule ever sees the string, so those
    alternatives can never fire. Rather than hand-edit 40 of them and leave
    the next author the same trap, the expansion is done here, once, at
    compile time: every literal alternative gains its normalised form
    alongside the original.

    Only alternatives that are pure literals are expanded. Anything carrying
    regex metacharacters is left exactly as written.
    """
    body = rx
    if body.startswith("^(") and body.endswith(")$"):
        body = body[2:-2]
    alts = body.split("|")
    out = list(alts)
    for alt in alts:
        if not re.fullmatch(r"[a-z0-9 .\-/']+", alt.replace("\\", "")):
            continue
        n = normalise(alt.replace("\\", ""))
        if n and n != alt and n not in out:
            out.append(re.escape(n).replace("\\ ", " "))
    return "^(" + "|".join(out) + ")$"


def check_rules():
    """Every rule alternative must be able to match a NORMALISED string.

    normalise() strips a trailing " protein" and the hedges, so a rule
    alternative written as "i-leader protein" can never fire -- the string
    reaches it as "i-leader". That is not a hypothetical: I_LEADER shipped
    empty out of pass 3 for exactly this reason while 234 deposits of
    "i-leader protein" sat in the pool.

    A dead alternative is invisible: the feature simply comes out smaller,
    or empty, and nothing errors. So it is checked rather than remembered.
    """
    bad = []
    for mod, rules in PSSM.items():
        for key, rx in rules.items():
            body = rx
            if body.startswith("^(") and body.endswith(")$"):
                body = body[2:-2]
            for alt in body.split("|"):
                lit = alt.replace("\\", "")
                if not re.fullmatch(r"[a-z0-9 .\-/']+", lit):
                    continue           # has regex metacharacters; skip
                n = normalise(lit)
                #  Only an ORPHAN matters: an alternative whose normalised
                #  form matches nothing in its own rule. "hexon protein"
                #  normalises to "hexon", which the same rule also lists,
                #  so it is redundant and harmless. "shutoff protein"
                #  normalises to "shutoff", which the rule does NOT list,
                #  so the feature can never bind -- and SHUTOFF did ship
                #  empty out of pass 3.
                if n != lit and not re.match(expand(rx), n, re.I):
                    bad.append((mod, key, lit, n))
    return bad


# --------------------------------------------------------------------------
def rescue(pool, hits, min_id=80.0, min_qcov=0.85):
    """Map each pooled sequence to a named feature it is homologous to.

    Coverage is of the QUERY, so a protein matching a collection member over
    only part of itself is rejected. The two numbers must agree with
    filter_unchar.py, or a protein is adopted into a feature AND left in the
    grab-bag, and both profiles then fire on one locus.
    """
    if not pool or not hits:
        return {}
    t = tempfile.mkdtemp(prefix="adeno.rescue.")
    try:
        db = os.path.join(t, "named.faa")
        n = 0
        with open(db, "w") as fh:
            seen = set()
            for k, rows in hits.items():
                for _g, s, _a, _gen in rows:
                    if (k, s) in seen:
                        continue
                    seen.add((k, s))
                    fh.write(">%s|%d\n%s\n" % (k, n, s)); n += 1
        if not n:
            return {}
        subprocess.run(["makeblastdb", "-in", db, "-dbtype", "prot",
                        "-out", os.path.join(t, "db")], capture_output=True)
        q = os.path.join(t, "q.faa")
        idx = {}
        with open(q, "w") as fh:
            for i, s in enumerate(pool):
                idx[str(i)] = s
                fh.write(">%d\n%s\n" % (i, s))
        r = subprocess.run(
            ["blastp", "-query", q, "-db", os.path.join(t, "db"),
             "-outfmt", "6 qseqid sseqid pident length qlen bitscore",
             "-max_target_seqs", "5", "-evalue", "1e-5",
             "-num_threads", "4"], capture_output=True, text=True)
        best = {}
        for line in r.stdout.splitlines():
            f = line.split("\t")
            if len(f) < 6:
                continue
            qid, sid, pid, alen, qlen, bits = f[0], f[1], float(f[2]), int(f[3]), int(f[4]), float(f[5])
            if pid < min_id or alen < min_qcov * qlen:
                continue
            if qid not in best or bits > best[qid][1]:
                best[qid] = (sid.rsplit("|", 1)[0], bits)
        return {idx[q]: k for q, (k, _b) in best.items()}
    finally:
        shutil.rmtree(t, ignore_errors=True)


def cluster(seqs, min_id=0.3, cov=0.7):
    """mmseqs easy-cluster. Returns [[seq, ...], ...], largest first."""
    if len(seqs) < 2:
        return [[s] for s in seqs]
    t = tempfile.mkdtemp(prefix="adeno.clu.")
    try:
        with open(t + "/q.faa", "w") as fh:
            for i, s in enumerate(seqs):
                fh.write(">%d\n%s\n" % (i, s))
        subprocess.run(["mmseqs", "easy-cluster", t + "/q.faa", t + "/c", t + "/tmp",
                        "--min-seq-id", str(min_id), "-c", str(cov),
                        "--cov-mode", "0", "--threads", "4"], capture_output=True)
        grp = collections.defaultdict(list)
        p = t + "/c_cluster.tsv"
        if not os.path.exists(p):
            return [[s] for s in seqs]
        for l in open(p):
            r, m = l.split()
            grp[r].append(seqs[int(m)])
        return sorted(grp.values(), key=len, reverse=True)
    finally:
        shutil.rmtree(t, ignore_errors=True)


def consolidate(hits, mod):
    """One homology cluster, one feature.

    THE DEPOSITED NAME IS A CLAIM. A sequence sits in whichever collection
    its string landed in, so a mislabelled deposit puts it under the wrong
    key and its PSSM then out-scores the correct one at the OTHER protein's
    locus. Majority of the cluster wins.

    EXCEPT between two members of the same PARALOG set -- see above. Those
    are distinct genes that are homologous by descent, so the clustering
    cannot arbitrate and the majority rule would simply empty the smaller
    one.
    """
    sets = PARALOG.get(mod, [])

    def paired(a, b):
        return any(a in s and b in s for s in sets)

    keys = [k for k in hits if hits[k]]
    if len(keys) < 2:
        return collections.Counter()
    rows = [(k,) + tuple(e) for k in keys for e in hits[k]]
    t = tempfile.mkdtemp(prefix="adeno.cons.")
    try:
        with open(t + "/q.faa", "w") as fh:
            for i, r in enumerate(rows):
                fh.write(">%d\n%s\n" % (i, r[2]))
        subprocess.run(["mmseqs", "easy-cluster", t + "/q.faa", t + "/c", t + "/tmp",
                        "--min-seq-id", "0.3", "-c", "0.7", "--cov-mode", "0",
                        "--threads", "4"], capture_output=True)
        grp = collections.defaultdict(list)
        p = t + "/c_cluster.tsv"
        if not os.path.exists(p):
            return collections.Counter()
        for l in open(p):
            r, m = l.split()
            grp[r].append(int(m))
    finally:
        shutil.rmtree(t, ignore_errors=True)
    moved = collections.Counter()
    for members in grp.values():
        if len(members) < 3:
            continue
        c = collections.Counter(rows[i][0] for i in members)
        if len(c) < 2:
            continue
        win = c.most_common(1)[0][0]
        for i in members:
            k = rows[i][0]
            if k == win or paired(k, win):
                continue
            hits[k].remove(rows[i][1:])
            hits[win].append(rows[i][1:])
            moved["%s -> %s" % (k, win)] += 1
    return moved


# --------------------------------------------------------------------------
def main():
    orphans = check_rules()
    if orphans:
        for mod, key, lit, n in orphans:
            sys.stderr.write("DEAD RULE  %s %s: %r reaches the rule as %r\n"
                             % (mod, key, lit, n))
        sys.exit("%d rule alternative(s) can never match. Fix them first --\n"
                 "a dead alternative is invisible: the feature just comes out\n"
                 "smaller, or empty, and nothing errors." % len(orphans))
    only = sys.argv[sys.argv.index("--only") + 1] if "--only" in sys.argv else None
    nowrite = "--no-write" in sys.argv

    md5 = [l.strip() for l in open(os.path.join(W, "%s.uniq.md5" % T))]
    ann = [l.rstrip("\n").split("\t") for l in open(os.path.join(W, "%s.uniq.id_ann" % T), errors="replace")]
    m2a = {m: " ".join((a[1] if len(a) > 1 else "").split()) for m, a in zip(md5, ann)}
    seq = {}
    for l in open(os.path.join(W, "%s.uniq.seq" % T), errors="replace"):
        f = l.rstrip("\n").split("\t")
        if len(f) >= 2:
            seq[f[0]] = f[1]

    #  Genomes dropped at intake by the SKILL §1 unclassified-sp. rule.
    dropped = set()
    dp = os.path.join(W, "intake", "DROPPED_sp.tsv")
    if os.path.exists(dp):
        dropped = {l.split("\t")[0] for l in open(dp)}

    g2gen, g2name = {}, {}
    for l in open(os.path.join(W, "%s.id_name_gs" % T), errors="replace"):
        f = l.rstrip("\n").split("\t")
        if len(f) >= 2:
            g2gen[f[0]] = (f[3] if len(f) >= 4 and f[3] else "(unassigned)")
            g2name[f[0]] = f[1]

    MISFILED = {}
    mp = os.path.join(W, "misfiled.json")
    if os.path.exists(mp):
        for g, d in json.load(open(mp)).items():
            MISFILED[g] = d["to"]
        print("reassigned by homology: %d genome(s)\n" % len(MISFILED))
    else:
        print("WARNING: no misfiled.json -- module assignment is taking the\n"
              "         source genus at face value. Run find_misfiled.py.\n")

    g2mod = {g: m for m, gs in GENERA.items() for g in gs}

    def modof(g):
        return MISFILED.get(g) or g2mod.get(g2gen.get(g, ""))

    syn_rows = []
    summary = []
    for mod in GENERA:
        if only and mod != only:
            continue
        rules = [(k, re.compile(expand(v), re.I)) for k, v in PSSM[mod].items()]
        hits = collections.defaultdict(list)
        unass = collections.Counter(); nmc = collections.Counter()
        disc = collections.Counter(); pooled = collections.Counter()
        pool = {}

        for l in open(os.path.join(W, "%s.id_md5" % T), errors="replace"):
            f = l.rstrip("\n").split("\t")
            if len(f) < 2 or not f[1].strip():
                continue
            g = ".".join(f[0].split("|")[-1].split(".")[:2])
            if g in dropped or modof(g) != mod:
                continue
            s = m2a.get(f[1], ""); sq = seq.get(f[1], "")
            if not s or not sq:
                continue
            gen = g2gen.get(g, "") or ""
            if DISC.search(s):
                disc[(s, gen)] += 1
                continue
            sl = normalise(s)
            key = None
            for k, rx in rules:
                if rx.match(sl):
                    key = k
                    break
            if key is not None:
                hits[key].append((g, sq, s, gen))
                continue
            if NOT_MODELLED.match(sl):
                nmc[(s, gen)] += 1
                continue
            #  Mass-named, ORF-numbered, or simply unmatched: all three go
            #  to the pool. The pool is where homology, not the string,
            #  decides.
            (pooled if POOL.match(sl) else unass)[(s, gen)] += 1
            if sq in pool:
                pool[sq][3] += 1
            else:
                pool[sq] = [g, s, sq, 1]

        n_pool_before = len(pool)
        resc = collections.Counter()
        for sq, key in rescue(list(pool), hits).items():
            g, a, _s, nf = pool.pop(sq)
            gen = g2gen.get(g, "") or ""
            hits[key].append((g, sq, a, gen))
            resc["%s -> %s" % (a[:30], key)] += nf
            for ctr in (pooled, unass):
                if ctr.get((a, gen)):
                    ctr[(a, gen)] -= nf
                    if ctr[(a, gen)] <= 0:
                        del ctr[(a, gen)]

        moved = consolidate(hits, mod)

        #  ---- FIBER: the explicit labels are POSITIONAL. Trust them. ----
        #  Most adenoviruses have one fiber gene. The enteric mastadeno-
        #  viruses (HAdV-F40/41, PMID_claude_generated 7975882) and many
        #  fowl adenoviruses (FAdV-1 CELO, PMID_claude_generated 17698649 /
        #  18561970) have TWO, side by side.
        #
        #  A FIRST ATTEMPT ROUTED THEM BY LENGTH -- longest is fiber-1 --
        #  AND IT WAS WRONG, in a way only measurement showed. Of the 215
        #  Aviadenovirus genomes depositing both, fiber-1 is the longer in
        #  just 39: Fowl aviadenovirus 4 deposits f1=431 and f2=479 in 99
        #  genomes, and only CELO (f1=793, f2=410) matches the intuition.
        #  "fiber-1" and "fiber-2" name the FIRST and SECOND gene along the
        #  genome, not the long and the short one. The length rule swapped
        #  181 deposits each way and would have trained both profiles on a
        #  mixture of the two genes.
        #
        #  "long fibre" and "short fibre" DO mean what they say -- long is
        #  longer in 247 of the 252 mastadenovirus genomes depositing both.
        #  Those bind by name already.
        #
        #  So an explicitly labelled deposit is left exactly where its label
        #  put it, and only a GENERIC "fiber" is routed -- by homology to
        #  the two labelled collections, which are 35-50% identical over the
        #  shaft and so separate cleanly. A genome carrying only one fiber
        #  deposit has neither a long nor a short one; it just has the
        #  fiber, and it goes to FIBER.
        fib = collections.Counter()
        if "FIBER1" in PSSM[mod]:
            bygen = collections.defaultdict(list)
            for k in ("FIBER", "FIBER1", "FIBER2"):
                for e in hits.get(k, []):
                    bygen[e[0]].append((k, e))
            generic = [(g, e) for g, rows in bygen.items() if len(rows) > 1
                       for k, e in rows if k == "FIBER"]
            if generic and hits.get("FIBER1") and hits.get("FIBER2"):
                sub = {"FIBER1": hits["FIBER1"], "FIBER2": hits["FIBER2"]}
                got = rescue({e[1]: None for _g, e in generic}, sub,
                             min_id=55.0, min_qcov=0.60)
                for g, e in generic:
                    k = got.get(e[1])
                    if not k:
                        continue
                    hits["FIBER"].remove(e)
                    hits[k].append(e)
                    fib["generic -> %s" % k] += 1
            #  A genome with ONE fiber-family deposit must not sit in
            #  FIBER1/FIBER2 -- there is no second gene to be first or
            #  second of.
            for g, rows in bygen.items():
                if len(rows) != 1:
                    continue
                k, e = rows[0]
                if k == "FIBER" or e not in hits[k]:
                    continue
                hits[k].remove(e)
                hits["FIBER"].append(e)
                fib["%s -> FIBER (single)" % k] += 1

        #  ---- PARALOGUE AUDIT: a name rule can file into the wrong sibling --
        #  Consolidation is forbidden between members of a PARALOG set, which
        #  stops homology merging two real genes. It does nothing about the
        #  opposite error: a deposit whose NAME puts it in the wrong sibling.
        #
        #  MEASURED, Barthadenovirus: genome 2913659.5 deposits a 382 aa
        #  protein called "LH1". OAdV-7 LH1 is 120 aa and LH3 is 382 aa, and
        #  the sequence hits LH3 at 41-47% identity over 368 of its 382
        #  residues. It is LH3. The name rule bound it on the string.
        #
        #  A move needs BOTH signals to agree, because either alone is weak:
        #  homology inside a paralogue family is exactly what cannot separate
        #  the members, and length alone cannot tell a short paralogue from a
        #  truncated one. So: the best cross-collection hit must beat the best
        #  within-collection hit, AND the length must sit near the sibling's
        #  reference and away from its own. Silence here is not success --
        #  every move is printed.
        par = collections.Counter()
        refs = REFLEN.get(mod, {})
        for group in PARALOG.get(mod, []):
            grp = [k for k in group if k in refs and hits.get(k)]
            if len(grp) < 2:
                continue
            for k in list(grp):
                for e in list(hits[k]):
                    L = len(e[1])
                    own = refs[k]
                    near = [j for j in grp
                            if j != k and abs(L - refs[j]) < 0.20 * refs[j]]
                    if not near or abs(L - own) <= 0.20 * own:
                        continue
                    #  The query is STILL A MEMBER of hits[k], so leaving it
                    #  in the subject set makes it match itself at 100% and
                    #  no move can ever happen -- which is why the first run
                    #  of this audit reported nothing at all while the 382 aa
                    #  "LH1" sat in plain sight. Drop the query sequence from
                    #  every collection before asking.
                    sub = {j: [r for r in hits[j] if r[1] != e[1]]
                           for j in grp}
                    sub = {j: r for j, r in sub.items() if r}
                    got = rescue({e[1]: None}, sub, min_id=35.0, min_qcov=0.70)
                    j = got.get(e[1])
                    if j and j != k and j in near:
                        hits[k].remove(e)
                        hits[j].append(e)
                        par["%s -> %s (%daa, ref %d vs %d)"
                            % (k, j, L, own, refs[j])] += 1

        #  ---- split precursor from mature peptide (see MATURE above) ----
        mat = collections.Counter()
        for src, dst in MATURE.get(mod, {}).items():
            rows = hits.get(src) or []
            bygen = collections.defaultdict(list)
            for e in rows:
                bygen[e[0]].append(e)
            for g, rs in bygen.items():
                if len(rs) < 2:
                    continue
                for a in rs:
                    for b in rs:
                        if a is b or len(b[1]) <= len(a[1]):
                            continue
                        #  a is shorter; cleaved product, or a downstream
                        #  start codon? Met says which.
                        if b[1].endswith(a[1]) and not a[1].startswith("M"):
                            if a in hits[src]:
                                hits[src].remove(a)
                                hits[dst].append(a)
                                mat["%s -> %s (%d of %d aa, %d removed)"
                                    % (src, dst, len(a[1]), len(b[1]),
                                       len(b[1]) - len(a[1]))] += 1
                            break

        #  ---- accessory discovery on whatever homology could not place ----
        #  The floor is TWO numbers, not one. Pass 1 used features alone and
        #  declared 152 UNC features in Mastadenovirus, 28 of which were a
        #  SINGLE sequence that happened to recur in five genomes. A profile
        #  cannot be trained from one sequence, and a feature declared with
        #  nothing behind it is worse than an omission: the annotator then
        #  reports it missing on every genome. So a cluster must carry
        #  ACC_SEQ_FLOOR distinct sequences AND ACC_FLOOR features.
        clusters = cluster([s for s in pool])

        #  ---- CLUSTER ADOPTION: a cluster may BE a feature already named ----
        #  The per-sequence rescue above asks each unplaced sequence on its own
        #  at >=80% identity, which a divergent genus never reaches. A whole
        #  cluster is far better evidence than any one of its members: if most
        #  of a coherent group hits one declared feature, the group is that
        #  feature, whatever the deposits called it.
        #
        #  MEASURED, Barthadenovirus: 13 sequences deposited as "E1B 55K",
        #  "E1B-2" and "E1B protein, large T-antigen" clustered together at
        #  370-393 aa and every one of them hits LH3 at 42-77% identity over
        #  the full length. They ARE LH3 -- atadenoviruses have no E1B at all,
        #  and the capsid structure paper (PMID_claude_generated 33789897)
        #  says mastadenovirus E1B-55K was exapted FROM the atadenovirus LH3,
        #  which is exactly why a submitter reaching for a familiar name
        #  reached for that one. Without this pass LH3 shipped a profile built
        #  from 4 sequences while 13 more sat in an unnamed cluster.
        #
        #  Also caught here: 4 Barthadenovirus DNA polymerases, two of them
        #  spelled "DNA polyermase", which no name rule will ever match.
        #
        #  The threshold is deliberately far below the per-sequence rescue,
        #  because the agreement requirement does the work instead: a simple
        #  majority of the cluster must land on the SAME feature.
        #  ONE blastp for the whole pool, then votes are tallied per cluster.
        #  Calling rescue() per cluster rebuilt the named-sequence BLAST
        #  database every time -- 365 rebuilds against 81,000 Mastadenovirus
        #  sequences -- and the pass did not finish in 25 minutes. The
        #  information is identical either way.
        adopt = collections.Counter()
        kept = []
        votes = rescue({s: None for s in pool}, hits,
                       min_id=40.0, min_qcov=0.75) if pool else {}
        for cl in clusters:
            got = {s: votes[s] for s in cl if s in votes}
            if len(got) * 2 <= len(cl):
                kept.append(cl)
                continue
            vote = collections.Counter(got.values())
            key, n = vote.most_common(1)[0]
            if n * 2 <= len(cl):
                kept.append(cl)
                continue
            lab = collections.Counter(pool[s][1] for s in cl).most_common(2)
            adopt["%s <- %s" % (key, "/".join(a for a, _ in lab))] += len(cl)
            for s in cl:
                g, a, _s, nfeat = pool.pop(s)
                gen = g2gen.get(g, "") or ""
                hits[key].append((g, s, a, gen))
                for ctr in (pooled, unass):
                    if ctr.get((a, gen)):
                        ctr[(a, gen)] -= nfeat
                        if ctr[(a, gen)] <= 0:
                            del ctr[(a, gen)]
        clusters = kept

        unc = {}
        n_unc = 0
        thin = []
        named_cl = {}
        for cl in clusters:
            nf = sum(pool[s][3] for s in cl)
            if len(cl) < ACC_SEQ_FLOOR or nf < ACC_FLOOR:
                thin.append((cl, nf))
                continue
            lab = collections.Counter(pool[s][1] for s in cl)
            #  Name the cluster for the gene the literature describes, where
            #  its own dominant deposit label names one. See CLUSTER_NAME.
            name = None
            #  Ignore the contentless labels when asking what a cluster is
            #  called. Siadenovirus ORF7 and ORF8 are two of the five
            #  genus-defining genes and BOTH failed to be named, because
            #  "hypothetical protein" outnumbered "ORF7" 7 to 4 inside the
            #  cluster. "hypothetical protein" is not a competing claim about
            #  identity; it is the absence of one.
            #  Vote on the NORMALISED label, not the raw string. "ORF7" x4 and
            #  "ORF7 protein" x2 are the same claim, and counted apart they
            #  split the vote and neither reaches a majority -- which is why
            #  ORF7 failed to be named while ORF8 succeeded.
            inf = collections.Counter()
            for a, n in lab.items():
                key = re.sub(r"^(protein|hypothetical protein)\s+", "",
                             normalise(a.lower().strip()))
                if re.match(r"^(hypothetical|unknown|unnamed|undefined|"
                            r"uncharacteri[sz]ed|putative)\b", key) \
                        or key in ("protein", "cds", "orf", ""):
                    continue
                inf[key] += n
            if inf:
                key, ntop = inf.most_common(1)[0]
                cand = (CLUSTER_NAME.get(mod, {}).get(key)
                        if ntop * 2 > sum(inf.values()) else None)
                top = key
                if cand and cand not in hits:
                    name = cand
                    named_cl[cand] = top
            if name is None:
                n_unc += 1
                name = "UNC%d" % n_unc
            unc[name] = (cl, nf, lab.most_common(3))
            for s in cl:
                g, a, _s, nfeat = pool.pop(s)
                hits[name].append((g, s, a, g2gen.get(g, "") or ""))
                #  The tracking file is a list of UNFINISHED work. A string
                #  whose sequence has now been adopted into a declared
                #  feature is finished, so it comes off both counters --
                #  pass 1 left ORF20 sitting in UNASSIGNED_TRACKING at 324
                #  features while UNC2 was built out of 72 of its sequences.
                gen = g2gen.get(g, "") or ""
                for ctr in (pooled, unass):
                    if ctr.get((a, gen)):
                        ctr[(a, gen)] -= nfeat
                        if ctr[(a, gen)] <= 0:
                            del ctr[(a, gen)]

        od = os.path.join(W, mod, "collections")
        if not nowrite:
            os.makedirs(od, exist_ok=True)
            #  REMOVE COLLECTIONS THIS RUN NO LONGER PRODUCES. Writing is not
            #  the same as replacing: pass 1 declared 152 Mastadenovirus UNC
            #  clusters and pass 11 declares 55, but UNC56..UNC152 were still
            #  sitting on disk from the earlier run, holding 1-2 sequences
            #  each -- below the floor, belonging to a clustering that no
            #  longer exists, and invisible in every summary this script
            #  prints. The next step in the workflow globs *.fasta and would
            #  have built profiles from all of them.
            #
            #  Only .fasta files are touched, and only in this module's own
            #  collections directory, which this script owns outright. Nothing
            #  recursive, nothing outside `od`, and synonyms.tsv is left alone.
            stale = sorted(f for f in os.listdir(od)
                           if f.endswith(".fasta") and f[:-6] not in hits)
            for f in stale:
                os.remove(os.path.join(od, f))
            if stale:
                print("  removed %d stale collection(s) from a previous run: %s"
                      % (len(stale), ", ".join(s[:-6] for s in stale[:12])
                         + (" ..." if len(stale) > 12 else "")))
            for k, rows in hits.items():
                with open(os.path.join(od, "%s.fasta" % k), "w") as fh:
                    seen = set()
                    for g, s, a, gen in rows:
                        if s in seen:
                            continue
                        seen.add(s)
                        fh.write(">%s %s\n%s\n" % (g, a, s))
            #  synonyms.tsv is the join between the deposited vocabulary and
            #  the feature each string ended up in. check_literature_gaps.py
            #  scores a probe list against THIS file, so a module without one
            #  cannot be literature-checked at all -- which is exactly how a
            #  protein the literature describes stays invisible. The schema
            #  is fixed by that script: count, annotation, genus, bins_to.
            syn = collections.Counter()
            for k, rows in hits.items():
                for _g, _s, a, gen in rows:
                    syn[(a, gen, k)] += 1
            for src, dest in ((unass, "UNASSIGNED"), (pooled, "UNASSIGNED"),
                              (nmc, "NOT_MODELLED"), (disc, "DISCARDED")):
                for (a, gen), n in src.items():
                    syn[(a, gen, dest)] += n
            with open(os.path.join(od, "synonyms.tsv"), "w") as fh:
                fh.write("count\tannotation\tgenus\tbins_to\n")
                for (a, gen, k), n in syn.most_common():
                    fh.write("%d\t%s\t%s\t%s\n" % (n, a, gen, k))
            with open(os.path.join(W, mod, "UNASSIGNED_TRACKING.tsv"), "w") as fh:
                fh.write("n_features\tgenus\tannotation\n")
                for (a, gen), n in (unass + pooled).most_common():
                    fh.write("%d\t%s\t%s\n" % (n, gen, a))
            with open(os.path.join(W, mod, "NOT_MODELLED.tsv"), "w") as fh:
                fh.write("n_features\tgenus\tannotation\n")
                for (a, gen), n in nmc.most_common():
                    fh.write("%d\t%s\t%s\n" % (n, gen, a))
            with open(os.path.join(W, mod, "DISCARDED_partial.tsv"), "w") as fh:
                fh.write("n_features\tgenus\tannotation\n")
                for (a, gen), n in disc.most_common():
                    fh.write("%d\t%s\t%s\n" % (n, gen, a))
            with open(os.path.join(W, mod, "ACCESSORY_CLUSTERS.tsv"), "w") as fh:
                fh.write("feature\tn_seqs\tn_features\tmedian_aa\ttop_labels\n")
                for name, (cl, nf, lab) in unc.items():
                    L = sorted(len(s) for s in cl)
                    fh.write("%s\t%d\t%d\t%d\t%s\n" %
                             (name, len(cl), nf, L[len(L) // 2],
                              "; ".join("%s x%d" % (a, n) for a, n in lab)))

        ngen = len({r[0] for rows in hits.values() for r in rows})
        print("######## %s  %d genomes" % (mod, ngen))
        print("  features bound     %d" % sum(len(v) for v in hits.values()))
        print("  declared features  %d  (%d named + %d UNC)"
              % (len(hits), len(hits) - n_unc, n_unc))
        print("  pooled by mass/ORF %d features, %d unique seqs"
              % (sum(pooled.values()), n_pool_before))
        print("  rescued by homology %d -> named features" % sum(resc.values()))
        print("  moved by consolidation %d" % sum(moved.values()))
        if fib:
            print("  fiber routed by genome %d  (%s)"
                  % (sum(fib.values()),
                     ", ".join("%s %d" % (k, n) for k, n in fib.most_common())))
        if mat:
            print("  precursor/mature split %d  (%s)"
                  % (sum(mat.values()), "; ".join(sorted(mat))))
        if par:
            print("  paralogue audit moved %d  (%s)"
                  % (sum(par.values()), "; ".join(sorted(par))))
        if adopt:
            print("  adopted whole clusters %d seqs  (%s)"
                  % (sum(adopt.values()),
                     ", ".join("%s %d" % (k, n) for k, n in adopt.most_common(8))))
        if named_cl:
            print("  clusters named from literature %d  (%s)"
                  % (len(named_cl),
                     ", ".join("%s<-%s" % (k, v) for k, v in sorted(named_cl.items()))))
        print("  declared as UNC    %d clusters, %d seqs"
              % (n_unc, sum(len(v[0]) for v in unc.values())))
        print("  below the floor    %d clusters, %d seqs, %d features"
              % (len(thin), sum(len(c) for c, _n in thin),
                 sum(n for _c, n in thin)))
        print("  still unplaced     %d unique seqs, %d features"
              % (len(pool), sum(pool[s][3] for s in pool)))
        print("  discarded partial  %d features" % sum(disc.values()))
        top = sorted(((len(v), k) for k, v in hits.items()), reverse=True)
        print("  collections: " + ", ".join("%s(%d)" % (k, n) for n, k in top))
        print()
        summary.append((mod, ngen, len(hits), n_unc, len(pool)))

    print("=" * 70)
    for mod, ngen, nf, nu, left in summary:
        print("%-18s %5d genomes  %3d features (%d UNC)  %d seqs unplaced"
              % (mod, ngen, nf, nu, left))


if __name__ == "__main__":
    main()
