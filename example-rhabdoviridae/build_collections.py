#!/usr/bin/env python3
"""
Build per-protein collections for the Rhabdoviridae LowVan module from the
BV-BRC dump in this directory.

Inputs (all already here):
    Rhabdo.id_name_gs   genome_id \t genome_name \t family \t genus
    Rhabdo.id_md5       feature_id \t md5
    Rhabdo.uniq.md5     md5                      (one per unique sequence)
    Rhabdo.uniq.seq     md5 \t sequence
    Rhabdo.uniq.id_ann  feature_id \t annotation

uniq.seq and uniq.id_ann are joined BY KEY -- md5 and feature id respectively,
the latter resolved through id_md5. They were once positionally aligned and
this loader relied on it, but that alignment is only a side effect of SOLR
returning rows in input order and it did not survive the re-export.

Outputs:
    synonyms.tsv            every annotation string x genus, with counts.
                            This is the file to read first -- it is the real
                            BV-BRC vocabulary the module is replacing, and it
                            is what the binning rules below are derived from.
    collections/<Module>/<FEATURE>.fasta
    collections/report.tsv          per-collection counts and length ranges
    collections/UNBINNED.tsv        annotation x genus a rule did not match
    collections/MISSING.tsv         genomes lacking each core protein
"""

import collections
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile

# --------------------------------------------------------------------------
# N-TERMINAL CURATION OVERRIDES
# --------------------------------------------------------------------------
# feature_id -> residues to strip from the N-terminus, with the reason.
#
# These are start-codon misassignments in BV-BRC, not biology. Applied here
# rather than by hand-editing collections/, so a rebuild does not silently
# revert them.
#
# Dichorhavirus G, three Orchid fleck dichorhavirus records (3 Sep 2026).
# BV-BRC calls these 578-580 aa where every other dichorhavirus G is 523-552.
# The extra residues are an upstream in-frame Met, on three lines of evidence:
#
#   1. Selection. Comparing the two long forms that differ (152177.57 vs .61),
#      the extension is 30.6% identical while the body downstream of the
#      internal Met is 88.0% identical. The extension is not under selection.
#   2. Alignment. All three carry an in-frame Met at exactly the position
#      where 14 other dichorhavirus G proteins begin (residue 39, 39, 37).
#   3. Hydrophobicity. The annotated N-terminus peaks at Kyte-Doolittle +1.56
#      (w=9) and is Tyr/His/Ile-rich with prolines; from the internal Met it
#      peaks at +3.28 with an LMLFIISG core, matching the other 14. Measured,
#      not assumed from RABV: 96% of the 90 plant Betarhabdovirinae G proteins
#      carry the same N-terminal h-region (median KD 2.88) at the same
#      strength as the animal set (Alpharhabdovirinae median 2.91).
#
# NOTE the cleavage POSITION remains inference -- no plant rhabdovirus has
# published signal-peptide coordinates, as Rhabdoviridae.notes records. What
# is established here is that the h-region exists and the extension sits
# upstream of it, which is all the trim rests on.
NTERM_TRIM = {
    "fig|152177.57.CDS.5": 38,
    "fig|152177.59.CDS.5": 38,
    "fig|152177.61.CDS.5": 36,
}

# --------------------------------------------------------------------------
# PER-FEATURE BINNING OVERRIDES
# --------------------------------------------------------------------------
# feature_id -> (correct key or None to drop, reason).
#
# For source annotations that are simply wrong about which protein this is.
# No regex can fix these: the offending string is a legitimate name for the
# feature it was binned into, it is just attached to the wrong sequence. They
# are listed individually so a rebuild does not silently reintroduce them, and
# so each one carries the evidence that settled it.
#
# Found by scripts/qc_cross_feature.py in the lowvan-module skill, which flags
# any cluster near-identical to a different feature's members.
MISBINNED = {
    # 298 aa, annotated "matrix", in a genome whose peg.1/peg.3 are N (451 aa)
    # and M (202 aa). Lyssavirus M is 202 aa and P is ~298: this is the
    # phosphoprotein. It is 99.7% identical to fig|38768.6.CDS.2, which the
    # same dump annotates "phosphoprotein".
    "fig|57482.650.peg.2": ("P", "298 aa 'matrix' in the P slot; real M is peg.3 at 202 aa"),

    # Nine rabies genomes whose only feature, CDS.1, is annotated
    # "nucleoprotein" but is 524 aa and 98-99% identical to members of the G
    # collection. Rabies N is 450 aa. Their PSSM out-scored the real N profile
    # at the glycoprotein locus, so the annotator emitted "Nucleocapsid
    # protein" on G's coordinates and never called N at all.
    #
    # This was first fixed by retiring the cluster by hand. The re-export
    # rebuilt the tree and the cluster came straight back under a new number,
    # which is why it belongs here instead: a curation decision that lives only
    # in a directory does not survive a rebuild.
    "fig|11292.22889.CDS.1": ("G", "524 aa 'nucleoprotein'; rabies N is 450, this is the glycoprotein"),
    "fig|11292.22893.CDS.1": ("G", "524 aa 'nucleoprotein'; rabies N is 450, this is the glycoprotein"),
    "fig|11292.22901.CDS.1": ("G", "524 aa 'nucleoprotein'; rabies N is 450, this is the glycoprotein"),
    "fig|11292.22907.CDS.1": ("G", "524 aa 'nucleoprotein'; rabies N is 450, this is the glycoprotein"),
    "fig|11292.22909.CDS.1": ("G", "524 aa 'nucleoprotein'; rabies N is 450, this is the glycoprotein"),
    "fig|11292.22915.CDS.1": ("G", "524 aa 'nucleoprotein'; rabies N is 450, this is the glycoprotein"),
    "fig|11292.22916.CDS.1": ("G", "524 aa 'nucleoprotein'; rabies N is 450, this is the glycoprotein"),
    "fig|11292.22924.CDS.1": ("G", "524 aa 'nucleoprotein'; rabies N is 450, this is the glycoprotein"),
    "fig|11292.22930.CDS.1": ("G", "524 aa 'nucleoprotein'; rabies N is 450, this is the glycoprotein"),

    # 513 aa "matrix protein" in a genome whose peg.1 is also "matrix protein"
    # at 279 aa and whose polymerase is peg.3. Bole Tick Virus 2 is N-P-M-G-L;
    # 513 aa in that slot is the glycoprotein.
    "fig|1608041.29.peg.2": ("G", "513 aa 'matrix protein'; the real M is peg.1 at 279 aa"),

    # 235 aa "polymerase protein" -- rabies L is 2,127 aa. 98% identical to a
    # 239 aa member of the P collection, so it is a phosphoprotein.
    "fig|11292.30343.peg.1": ("P", "235 aa 'polymerase protein'; rabies L is 2127 aa"),

    # "putative L3/G1" binds to L on the L3, but the polymerase in this genome
    # is peg.10 at 2,140 aa. It is 98% identical to a member of the G
    # collection: the submitter's G1, not an L.
    "fig|1481139.13.peg.8": ("G", "'putative L3/G1'; the real L is peg.10 at 2140 aa"),
}

AMBIG = set("BJOUXZ*")
MIN_AA = 30

# --------------------------------------------------------------------------
# genus -> module. Routing is by architecture, not taxonomy.
# --------------------------------------------------------------------------
ALPHA = "Alpharhabdovirinae"
BETA = "Betarhabdovirinae"
GYMNO = "Alphagymnorhavirus"

ALPHA_GENERA = [
    # Alpharhabdovirinae
    "Lyssavirus", "Vesiculovirus", "Ephemerovirus", "Sigmavirus", "Hapavirus",
    "Tibrovirus", "Ledantevirus", "Sripuvirus", "Almendravirus", "Curiovirus",
    "Tupavirus", "Sunrhavirus", "Barhavirus", "Sawgrhavirus", "Arurhavirus",
    "Perhabdovirus", "Sprivivirus", "Ohlsrhavirus", "Lostrhavirus",
    "Caligrhavirus", "Alphapaprhavirus", "Alphanemrhavirus",
    "Alphathriprhavirus", "Amplylivirus", "Merhavirus",
    "Alpharicinrhavirus", "Siniperhavirus",
    # Deltarhabdovirinae -- all carry the canonical five
    "Primrhavirus", "Stangrhavirus", "Betaricinrhavirus", "Betapaprhavirus",
    "Betanemrhavirus", "Alphahymrhavirus", "Gammahymrhavirus",
    "Gammaricinrhavirus",
    # Platrhaviruses. All three clear the 5-genome cutoff in this dump, so
    # they are routed here. NB Betaplatrhavirus commonly lacks G -- the JSON
    # therefore must not carry copy_num on G if these stay in.
    "Alphaplatrhavirus", "Betaplatrhavirus", "Gammaplatrhavirus",
    # Below the cutoff: core five folded in, no dedicated accessory PSSMs.
    "Scophrhavirus", "Cetarhavirus", "Alphacrustrhavirus",
    "Alphadrosrhavirus", "Mousrhavirus", "Betahymrhavirus", "Zarhavirus",
    "Replylivirus", "Uniorhavirus", "Betathriprhavirus", "Margarhavirus"]
BETA_GENERA = [
    "Alphacytorhabdovirus", "Betacytorhabdovirus", "Gammacytorhabdovirus",
    "Alphanucleorhabdovirus", "Betanucleorhabdovirus",
    "Gammanucleorhabdovirus", "Deltanucleorhabdovirus"]

GENUS_MODULE = {}
for _g in ALPHA_GENERA:
    GENUS_MODULE[_g.lower()] = ALPHA
for _g in BETA_GENERA:
    GENUS_MODULE[_g.lower()] = BETA
for _g in ("Novirhabdovirus", "Dichorhavirus", "Trirhavirus"):
    GENUS_MODULE[_g.lower()] = _g

# Genera deliberately NOT modelled. Their sequences are recorded in
# collections/NOT_MODELLED.tsv rather than binned, so the decision stays
# visible and is trivially reversible -- put the genus back in the routing
# above and restore its TAXA block in build_rhabdoviridae_json.py.
#
# Varicosavirus, dropped 3 Sep 2026 per JJD: ICTV shows accessory genes that
# cannot be supported with PSSMs from what BV-BRC holds. 262 genomes yielded
# 73 in-range proteins over six collections, 53 of them the coat protein, and
# 196 of the 260 genomes with any feature carry exactly one. L had three
# sequences, all Lettuce big-vein associated virus at 2040 aa -- under the
# -m 5 clustering minimum -- and ORF2/ORF4/ORF5/ORF6 were never buildable.
NOT_MODELLED = {
    "varicosavirus": "accessory genes not supportable with BV-BRC sequences",
    # Dropped 10 Sep 2026 after the re-export, which gave the genus its first
    # sequences and made the problem measurable rather than hypothetical.
    # 12 genomes, 4 species, 37 proteins. The species' proteins are only
    # 30-36% identical to one another -- a matrix protein from Alnus and one
    # from Chrysanthemum trirhavirus share 36% over 152 aa. MMseqs groups
    # nothing at the documented 0.8, nor anywhere down to 0.3 with coverage
    # relaxed to 0.5, for G, M or P. Five of nine features would build and
    # four would not, so the module could never call a whole genome.
    "trirhavirus": "4 species only 30-36% identical; G/M/P unclusterable at any threshold",
    # Dropped 11 Sep 2026 per JJD, covering both gymnorhavirus genera, which
    # shared one module. 16 genomes (14 Alpha, 2 Beta) yielded two collections
    # -- L with 7 sequences and N with 8 -- and neither forms a cluster at the
    # -mi 0.6 floor, so the module had SIX declared features and ZERO profiles.
    # Four of those six (P2, P3, P4, P5) never had a collection at all: they
    # were positional labels invented for one protein family, the U-number trap
    # in plain sight.
    #
    # What settled it was running a real genome rather than reading the table.
    # NC_139364.1 (Cupressus virus 1) matched the module's own rep contig and
    # came back with ONE feature, a nucleocapsid, and that came from a stale
    # profile left in the installed repo by an install step that never prunes.
    # An empty module is worse than no module: the rep contig captures the
    # genome and dead-ends it, where routing would otherwise have offered it to
    # Alpharhabdovirinae and its 225 profiles.
    #
    # Reversible the same way Trirhavirus is -- restore these two GENUS_MODULE
    # lines and the TAXA block in build_rhabdoviridae_json.py. The blocker is
    # sequence depth, not taxonomy: 8 nucleocapsids cannot make a profile.
    "alphagymnorhavirus": "6 features, 0 profiles; L=7 and N=8 will not cluster at -mi 0.6",
    "betagymnorhavirus": "shared the Alphagymnorhavirus module; 2 genomes",
}
#  Retired with the module, kept for the revival path:
#    GENUS_MODULE["alphagymnorhavirus"] = GYMNO
#    GENUS_MODULE["betagymnorhavirus"] = GYMNO
# Pre-split NCBI names still live in this dump.
GENUS_MODULE["nucleorhabdovirus"] = "LEGACY_Nucleorhabdovirus"
GENUS_MODULE["cytorhabdovirus"] = "LEGACY_Cytorhabdovirus"
# 1728 genomes have NO genus assigned ("Rhabdoviridae sp.", "Farmington",
# "Apis rhabdovirus 1", ...). They are held separately rather than guessed
# into a module: assigning them needs a BLASTn call against the rep contigs,
# which is a curation step, not a regex.
UNASSIGNED = "UNASSIGNED_no_genus"

BELOW_CUTOFF = {
    "margarhavirus": 2, "mousrhavirus": 3, "zarhavirus": 2,
    "replylivirus": 2, "uniorhavirus": 2, "scophrhavirus": 4,
    "cetarhavirus": 4, "betathriprhavirus": 2, "betagymnorhavirus": 2,
    "alphacrustrhavirus": 4, "alphadrosrhavirus": 4, "betahymrhavirus": 3,
}

# --------------------------------------------------------------------------
# Core rules. Derived from the actual synonyms census in this dump, not
# guessed -- 327 distinct strings, dominated by a handful of forms.
# (feature_key, gene_or_None, product_regex)
# --------------------------------------------------------------------------
# Note the deliberate typo tolerance. The real BV-BRC vocabulary in this dump
# contains "Nucleicapsid" (127x), "nucleotprotein" (27x), "gylcoprotein"
# (24x), "phosophoprotein" (10x), "rNA dependent RNA polymerase" (38x). Those
# are not edge cases, they are hundreds of proteins, so the patterns below
# match them rather than leaving them for a human to sweep up.
CORE = [
    # M1 / M2 FIRST, before anything else. M1 is the phosphoprotein and M2 is
    # the matrix protein. ICTV: novirhabdovirus P "formerly designated M1",
    # M "formerly designated M2", old gene order 3'-N-M1-M2-G-NV-L-5'.
    # The RefSeq reference record states it outright, so the mapping is not an
    # inference. NC_001542 (Rabies lyssavirus, Pasteur virus) annotates
    #
    #     1514..2407   /product="phosphoprotein M1"   297 aa
    #     2496..3104   /product="M2 protein"          202 aa
    #
    # -- "phosphoprotein M1" names both in one string, and M2 sits at the
    # matrix position and the matrix length. Verified against the record, not
    # recalled: an earlier version of this comment cited UniProt from memory.
    #
    # Note which rule actually fires for rabies: that product contains the
    # word "phosphoprotein", so it binds on the phosphoprotein pattern below,
    # not here. The bare ^M1\b rule catches only records annotated "M1
    # protein" with no other clue, and every one of the 20 such proteins in
    # this dump is Novirhabdovirus, 222-230 aa, blasting 20/20 to P.
    #
    # These MUST precede the generic M pattern. "M1 matrix protein" and
    # "M1 maxtrix protein" both contain the word matrix, so without this
    # ordering they bind to M -- i.e. silently to the wrong protein. Matching
    # on ^M1\b / ^M2\b rather than the exact string catches every variant in
    # the dump: "M1 protein", "M1 matrix protein", "M1 maxtrix protein".
    # These bin to the ordinary P and M collections and therefore carry the
    # conventional S1 Table annotations, not anything bespoke.
    #  Also the spelled-out forms. JJD caught "matrix protein 1" and "matrix
    #  protein 2" both sitting in Novirhabdovirus M: the ^M1\b pattern only
    #  anticipated "M1 protein" / "M1 matrix protein", so the written-out
    #  version fell straight through to the generic ^matrix rule below.
    #  Measured on Hirame novirhabdovirus, each blasted against the
    #  Novirhabdovirus collections with itself excluded:
    #
    #    "matrix protein 1"  227 aa  95.6% to P, NO hit to M   -> phosphoprotein
    #    "matrix protein 2"  193 aa  92.7% to M, NO hit to P   -> matrix
    #
    #  which is the same M1=P / M2=M mapping as the abbreviations, and the
    #  reason both must sit ahead of the ^matrix rule.
    ("P", None, r"^M1\b|^matrix protein[- ]?1\b"),
    ("M", None, r"^M2\b|^matrix protein[- ]?2\b"),
    #  The L synonym list is long because submitters have many words for the
    #  polymerase, and every one that is missing costs real proteins. The
    #  homology rescue surfaced 48 in this dump that no rule matched: "large
    #  polymerase protein" x15, "viral polymerase" x14, "rna-directed rna
    #  polymerase" x8 (DIRECTED, not dependent -- one word), "polyprotein" x6,
    #  "replicase" x5. Every one had a clear annotation and was lost to regex
    #  coverage alone.
    #
    #  Prefer the submitter's own word to a homology score wherever the word is
    #  unambiguous. "polyprotein" and "replicase" are broad terms, but nothing
    #  in a rhabdovirus except L runs 1800-2300 aa, and EXPECTED sends anything
    #  outside that range to L.outliers.fasta rather than into the collection.
    #  The length gate is what makes a loose word safe here.
    ("L", None, r"RNA[- ]?dependent RNA[- ?]?polymerase|^rNA dependent"
               r"|RNA[- ]?directed RNA[- ?]?polymerase"
               r"|^L (protein|polymerase)|^polymerase( protein)?$"
               r"|^RNA polymerase$|^large (structural )?protein$"
               r"|^L polymerase|polymerase L|RdRp|polymerase activity module"
               r"|^RNA-dependent RNA-polymerase$|^putative L|^RNA replicase"
               r"|^protein L$|^L$"
               r"|^large polymerase( protein)?$|^viral polymerase$"
               r"|^(putative )?replicase$|^polyprotein$"
               r"|^(putative )?viral RNA polymerase$"),
    # "outer coat protein" is the ENVELOPE glycoprotein, not the nucleocapsid.
    # It is Sigmavirus's word for G (176 proteins in this dump). It used to be
    # in the N rule below, which sits ahead of G, so all 176 were binned as N;
    # the cross-feature QC caught it when G cluster 16's master turned out to
    # be 99.8% identical to a member of the N collection. Plain "coat protein"
    # really is the nucleocapsid (Varicosavirus, Lyssavirus) and stays with N.
    ("G", None, r"^outer coat protein$"),
    ("N", None, r"nucle[oi][ct]?[- ]?(protein|capsid)|nucleocapsid"
               r"|^N( protein)?$|^coat protein$"
               r"|^putative N\b|^nucleoprotein"),
    ("P", None, r"phos[po]*phoprotein|polymerase[- ]associated"
               r"|^P( protein)?$|^PP3( protein)?$|^NS( protein)?$"
               r"|^putative P[- ]protein$|^phosphoprotein"),
    ("M", None, r"^matrix|matrix protein|^M( protein)?$|^putative M\b"),
    # "gylcoprotein" (24x) and "glcyoprotein" transpose different letters, so
    # they are listed explicitly rather than fudged with a character class.
    ("G", None, r"glycoprotein|gylcoprotein|glcyoprotein|glycprotein"
               r"|^G[- ]protein( precursor)?$|^G$|^spike protein$"
               r"|^putative G[- ]protein$")]

# Catch-all for the unnamed accessory ORFs. 152 in Hapavirus, 38 in
# Curiovirus, 27 in Tibrovirus... these are real proteins with no name, and
# S1 Table's "Uncharacterized lineage-specific protein" (empty symbol) is
# exactly the annotation for them. Evaluated LAST, after every named rule.
UNCHAR_RULE = ("UNCHAR", r"^hypothetical protein$|^unknown|^uncharacteri[sz]ed"
                         r"|^putative uncharacterized|^unnamed|^ORF\d*$"
                         r"|^protein of unknown function")

# Accessory rules keyed on GENUS, evaluated before CORE.
GENUS_RULES = {
    "vesiculovirus": [
        ("VESI_CPRIME", r"^C[-'] ?prime|C prime|^C'"),
        ("VESI_C", r"^C protein$|^C$")],
    "ephemerovirus": [
        # "non-structural glycoprotein" and "non-structural transmembrane
        # glycoprotein" are GNS too. The old pattern only caught the
        # unhyphenated spelling, so 21 BEFV GNS proteins matched the
        # generic G rule instead -- caught by qc_cross_feature.py when the
        # EPHEM_GNS cluster came back 95% identical to a G member.
        #  "non-structural transmembrane protein", without "glyco", is the same
        #  feature -- two Porcine ephemerovirus entries use it. The rule is
        #  genus-gated to Ephemerovirus, so the looser wording cannot reach
        #  anything else.
        ("EPHEM_GNS", r"\bGNS\b|non-?\s?structural\s+transmembrane\s+(?:glyco)?protein"
                      r"|non-?\s?structural\s+glycoprotein"),
        ("EPHEM_ALPHA1", r"alpha[- ]?1\b|^a1 "),
        ("EPHEM_ALPHA2", r"alpha[- ]?2\b|^a2 "),
        # BEFV has a third alpha ORF (~5.7 kDa) overlapping alpha2
        # [Walker 1997, PMID 9191923]; 8 records carry it in this dump.
        ("EPHEM_ALPHA3", r"alpha[- ]?3\b|^a3 "),
        ("EPHEM_BETA", r"^beta protein$|^beta$"),
        ("EPHEM_GAMMA", r"^gamma protein$|^gamma$")],
    # M1/M2 are handled family-wide at the top of CORE, so nothing genus
    # -specific is needed here.
    "novirhabdovirus": [("NV", r"non[- ]?virion|^NV\b")],
    "sigmavirus": [("SIGMA_X", r"^PP3 protein$|^PP3$|^X protein$")],
        "hapavirus": [("HAPA_VIROPORIN", r"viroporin"), ("HAPA_U3", r"^U3\b"),
                  ("HAPA_U2", r"^U2\b"), ("HAPA_U1", r"^U1\b")],
                    "curiovirus": [("CURIO_VIROPORIN", r"viroporin"),
                   ("CURIO_PMIP", r"^U1\b")],
        "sunrhavirus": [("SUNRHA_SH", r"small hydrophobic|^SH\b|^U1\b")],
        "sawgrhavirus": [("SAWGR_GY", r"^Gy\b")],
                            #  Retired 11 Sep 2026 with the module. P2-P5 were positional labels with
#  no collection behind any of them; see NOT_MODELLED above.
#    "alphagymnorhavirus": [("P5", r"^P5\b"), ("P4", r"^P4\b"),
#                           ("P3", r"^P3\b"), ("P2", r"^P2\b")],
}
#GENUS_RULES["betagymnorhavirus"] = GENUS_RULES["alphagymnorhavirus"]

MOV = ("MOV", r"movement protein|^sc4\b|^4b\b|30K|^P3( protein)?$"
              r"|^ORF ?3$|^3$")
for _g in ("alphacytorhabdovirus", "betacytorhabdovirus",
           "gammacytorhabdovirus", "alphanucleorhabdovirus",
           "betanucleorhabdovirus", "deltanucleorhabdovirus",
           "dichorhavirus", "cytorhabdovirus", "nucleorhabdovirus"):
    GENUS_RULES.setdefault(_g, []).append(MOV)
GENUS_RULES["alphacytorhabdovirus"] += [("ACYTO_VIROPORIN", r"viroporin")]
GENUS_RULES["betacytorhabdovirus"] += [("BCYTO_P5", r"^P5\b"),
                                       ("BCYTO_P4", r"^P4\b")]
GENUS_RULES["alphanucleorhabdovirus"] += [
    ("ANUCLEO_X", r"^X protein$")]
GENUS_RULES["betanucleorhabdovirus"] += [("BNUCLEO_U", r"^U protein$")]
# Gammanucleorhabdovirus P3/P4 are explicitly NOT movement proteins, so the
# MOV rule is deliberately absent for this genus.
GENUS_RULES["gammanucleorhabdovirus"] = []
GENUS_RULES["varicosavirus"] = [
    ("ORF6", r"^(protein|ORF) ?6$"), ("ORF5", r"^(protein|ORF) ?5$"),
    ("ORF4", r"^(protein|ORF) ?4$"),
    ("MOV", r"^(protein|ORF) ?3$|movement protein|30K"),
    ("ORF2", r"^(protein|ORF) ?2$")]
# Segment 3's P6/P7/P8 are the ordinary phosphoprotein, matrix and
# glycoprotein -- the numbers are the segment layout, not protein names -- so
# they bin to the conventional core keys and the JSON emits the conventional
# annotations. The positional numbering survives in the regex, which is where
# it belongs, and in the S2 accessory keys where no functional name is known.
# Trirhavirus rules retired with the module (see NOT_MODELLED). Kept here,
# commented, because the genus is real and a deeper sequence set may yet
# support it -- the blocker is diversity, not taxonomy.
#   ("G", r"^P8\b"), ("M", r"^P7\b"), ("P", r"^P6\b"),
#   ("ORF5", r"silencing suppressor|^P5\b"), ("ORF4", r"^P4\b"),
#   ("ORF3", r"^P3\b"), ("ORF2", r"^P2\b")

CORE_KEYS = {"N", "P", "M", "G", "L"}

def compile_rules(genus):
    out = []
    for key, prx in GENUS_RULES.get(genus, []):
        out.append((key, re.compile(prx, re.I)))
    for key, _g, prx in CORE:
        out.append((key, re.compile(prx, re.I)))
    out.append((UNCHAR_RULE[0], re.compile(UNCHAR_RULE[1], re.I)))
    return out

COMPILED = {g: compile_rules(g) for g in
            set(list(GENUS_RULES) + list(GENUS_MODULE))}
CORE_ONLY = ([(k, re.compile(p, re.I)) for k, _g, p in CORE]
             + [(UNCHAR_RULE[0], re.compile(UNCHAR_RULE[1], re.I))])

# 1728 genomes in this dump have no genus assigned, and they are not all one
# kind of thing -- they include ephemeroviruses, novirhabdoviruses and plant
# viruses. Giving them core rules only stranded their accessory proteins, so
# they get an explicit rule set of the accessories that are unambiguous
# WITHOUT knowing the genus. Deliberately excluded: everything keyed on a
# U-number, because a U-number cannot be interpreted without the genus.
UNASSIGNED_RULES = [
    ("NV", r"non[- ]?virion|^NV\b"),
    # "non-structural glycoprotein" and "non-structural transmembrane
        # glycoprotein" are GNS too. The old pattern only caught the
        # unhyphenated spelling, so 21 BEFV GNS proteins matched the
        # generic G rule instead -- caught by qc_cross_feature.py when the
        # EPHEM_GNS cluster came back 95% identical to a G member.
        #  "non-structural transmembrane protein", without "glyco", is the same
        #  feature -- two Porcine ephemerovirus entries use it. The rule is
        #  genus-gated to Ephemerovirus, so the looser wording cannot reach
        #  anything else.
        ("EPHEM_GNS", r"\bGNS\b|non-?\s?structural\s+transmembrane\s+(?:glyco)?protein"
                      r"|non-?\s?structural\s+glycoprotein"),
    ("EPHEM_ALPHA1", r"alpha[- ]?1\b"),
    ("EPHEM_ALPHA2", r"alpha[- ]?2\b"),
    ("EPHEM_ALPHA3", r"alpha[- ]?3\b"),
    ("EPHEM_BETA", r"^beta protein$|^beta$"),
    ("EPHEM_GAMMA", r"^gamma protein$|^gamma$"),
    ("MOV", r"movement protein|^sc4\b|^4b\b|30K")]
COMPILED[""] = ([(k, re.compile(p, re.I)) for k, p in UNASSIGNED_RULES]
                + CORE_ONLY)

# --------------------------------------------------------------------------
# Expected full-length ranges, for the "look at them for length" step.
# A sequence outside its range is not discarded -- it is written to
# <FEATURE>.outliers.fasta so it stays visible. Most outliers are partial
# gene submissions (this dump is full of single-gene lyssavirus records), but
# some are genuine misannotations: a 2128 aa "glycoprotein" is an L protein.
# Ranges are generous on purpose; tighten after the first alignment pass.
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# Length rescue for the unnamed.
#
# The binning rules read annotation text, so a protein BV-BRC labelled only
# "hypothetical protein" lands in UNCHAR no matter what it is. Nine of them are
# L polymerases: Nishimuro ledantevirus at 2111 aa, two Metorhabdovirus at 2108
# and 1528, Sphaeridiorhabdovirus at 2109, and so on -- and none of those nine
# genomes had any L binned at all, so the collection was simply missing them.
#
# Length settles it without reading anything. Excluding L, the longest protein
# in any Rhabdoviridae collection is 814 aa; L runs 1804-2295. Nothing sits
# between, so a 1400 aa floor cannot take anything but an L.
#
# The rescue routes to "L" rather than writing L.fasta directly, so EXPECTED
# below still applies: full-length ones join the collection, the four partials
# (1452-1637 aa) go to L.outliers.fasta where partials belong. That is the
# whole point of routing instead of appending -- a rescued sequence gets the
# same length scrutiny as one that was named correctly.
# --------------------------------------------------------------------------
UNCHAR_L_RESCUE = 1400

# --------------------------------------------------------------------------
# The uncharacterized bag must not be one collection.
#
# "hypothetical protein", "uncharacterized protein", a blank -- these strings
# say the submitter did not know, and pooling on them builds one profile out of
# proteins that have nothing to do with each other. Alpharhabdovirinae UNCHAR
# was 444 sequences recalling 31% of itself, which is what a model of nothing
# scores.
#
# So the pool is split by homology into UNC1, UNC2 ... UNCn: each a real
# homology group, each its own feature key with its own alignment directory and
# its own profiles. Per JJD, they all carry the SAME annotation string, the S1
# Table form for an uncharacterized lineage-specific protein, because that is
# what is known about them. The point of separate keys is that the day one of
# them gets a function from the literature, it is renamed in the JSON on its own
# and the others are untouched -- impossible while they share a collection.
#
# Grouped at the same 80% line used to decide that two proteins are the same
# protein. A group needs UNC_MIN members to become a feature; below that there
# is no profile to build, and the sequences are left out.
#
# The numbering is by group size, descending, with the representative sequence
# as tiebreak so a rebuild on unchanged data reproduces it. It is NOT stable
# across a data update: add sequences and UNC3 may become UNC2. The number is a
# handle, never a claim -- treat it exactly as cautiously as a U-number.
# --------------------------------------------------------------------------
UNC_PIDENT = 0.80
UNC_COV = 0.60
UNC_MIN = 3

def _drop_named_duplicates(module, pool, named):
    """Remove pool members that already ARE a named feature of this module.

    Exact-sequence dedup does not catch them: the same protein from two
    genomes differs by a residue or two. blastp does. Anything at or above the
    adoption threshold is that feature, and leaving it in the uncharacterized
    pool would put two profiles on one locus.
    """
    if not named or not pool or not shutil.which("blastp"):
        return pool
    items = list(pool.items())
    tmp = tempfile.mkdtemp(prefix="uncdup.")
    try:
        db = os.path.join(tmp, "named.faa")
        with open(db, "w") as fh:
            for (_m, k), d in sorted(named.items()):
                for i, s in enumerate(d):
                    fh.write(">%s@%d\n%s\n" % (k, i, s))
        subprocess.run(["makeblastdb", "-in", db, "-dbtype", "prot",
                        "-out", os.path.join(tmp, "db")], capture_output=True)
        q = os.path.join(tmp, "q.faa")
        with open(q, "w") as fh:
            for i, (s, _h) in enumerate(items):
                fh.write(">q%d\n%s\n" % (i, s))
        r = subprocess.run(
            ["blastp", "-query", q, "-db", os.path.join(tmp, "db"), "-outfmt",
             "6 qseqid pident length qlen", "-evalue", "1e-5",
             "-max_target_seqs", "3", "-num_threads", "4"],
            capture_output=True, text=True)
        drop = set()
        for line in r.stdout.splitlines():
            qid, pid, alen, qlen = line.split("\t")
            if float(pid) >= RESCUE_PIDENT and int(alen) / max(int(qlen), 1) >= RESCUE_COV:
                drop.add(int(qid[1:]))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if drop:
        print("%s: %d uncharacterized protein(s) are already a named feature, "
              "removed before grouping" % (module, len(drop)))
    return {s: h for i, (s, h) in enumerate(items) if i not in drop}

def reroute_length_outliers(bins):
    """Re-home proteins whose text rule put them in a feature their length denies.

    A text rule can be right about the words and wrong about the protein.
    Eight Alpharhabdovirinae proteins annotated "polymerase-associated
    protein" are 2114-2160 aa: the phrase contains "polymerase-associated",
    which is a phosphoprotein synonym, so they bound to P -- and then EXPECTED
    parked them in P.outliers.fasta because no phosphoprotein is 2000 aa. They
    blast to L at 54-98% identity over their whole length. They are
    polymerases, and they were being lost from L rather than merely misfiled.

    So: for anything outside its own feature's expected range, ask homology
    where it belongs. Move it only if another feature of the same module
    claims it at the adoption threshold AND its length fits that feature's
    range. Otherwise leave it where it is, to end up in <FEAT>.outliers as
    before -- the outlier file is for genuine partials, not for proteins we
    could have identified.
    """
    if not shutil.which("blastp"):
        return []
    moved = []
    for module in sorted({m for m, _k in bins}):
        keys = [k for (m, k) in bins if m == module and not k.startswith("UNC")]
        stray = []
        for k in keys:
            lo, hi = EXPECTED.get(k, (0, 10 ** 6))
            for s, h in bins[(module, k)].items():
                if not (lo <= len(s) <= hi):
                    stray.append((k, s, h))
        if not stray:
            continue
        tmp = tempfile.mkdtemp(prefix="reroute.")
        try:
            db = os.path.join(tmp, "named.faa")
            with open(db, "w") as fh:
                for k in keys:
                    lo, hi = EXPECTED.get(k, (0, 10 ** 6))
                    for i, s in enumerate(bins[(module, k)]):
                        if lo <= len(s) <= hi:      # only in-range members vote
                            fh.write(">%s@%d\n%s\n" % (k, i, s))
            subprocess.run(["makeblastdb", "-in", db, "-dbtype", "prot",
                            "-out", os.path.join(tmp, "db")], capture_output=True)
            q = os.path.join(tmp, "q.faa")
            with open(q, "w") as fh:
                for i, (_k, s, _h) in enumerate(stray):
                    fh.write(">s%d\n%s\n" % (i, s))
            r = subprocess.run(
                ["blastp", "-query", q, "-db", os.path.join(tmp, "db"), "-outfmt",
                 "6 qseqid sseqid pident length qlen bitscore", "-evalue", "1e-5",
                 "-max_target_seqs", "5", "-num_threads", "4"],
                capture_output=True, text=True)
            best = {}
            for line in r.stdout.splitlines():
                qid, sid, pid, alen, qlen, bits = line.split("\t")
                i = int(qid[1:])
                cand = (sid.split("@")[0], float(pid), int(alen) / max(int(qlen), 1),
                        float(bits))
                if i not in best or cand[3] > best[i][3]:
                    best[i] = cand
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        for i, (k, s, h) in enumerate(stray):
            if i not in best:
                continue
            tgt, pid, cov, _b = best[i]
            if tgt == k or pid < RESCUE_PIDENT or cov < RESCUE_COV:
                continue
            lo, hi = EXPECTED.get(tgt, (0, 10 ** 6))
            if not (lo <= len(s) <= hi):
                continue
            del bins[(module, k)][s]
            bins[(module, tgt)].setdefault(s, h)
            moved.append((module, k, tgt, len(s), pid, h.split()[0]))
    return moved

# --------------------------------------------------------------------------
# Proteins recovered from genomes, not from the protein export.
#
# BV-BRC holds a genome without holding all of its proteins. Twenty-one
# Alpharhabdovirinae phosphoproteins are in no protein record anywhere in the
# dump, yet the gene is plainly there: ORF-calling the interval between the
# called N and M finds it, full length, matching the P collection by blastp.
#
# Anything in Extracted/<Module>.<KEY>.faa is merged into that collection. The
# directory is the durable form of "we went and got it" -- hand-adding the
# sequences to collections/ would be erased by the next rebuild, which is the
# same trap as every other curation-by-file in this pipeline.
#
# The sequences still go through EXPECTED afterwards, so a recovered protein of
# the wrong length lands in <FEATURE>.outliers.fasta like any other.
# --------------------------------------------------------------------------
def merge_extracted(bins, here):
    import glob as _glob
    added = []
    for fa in sorted(_glob.glob(os.path.join(here, "Extracted", "*.faa"))):
        base = os.path.basename(fa)[:-4]
        if base.count(".") < 1:
            continue
        module, key = base.split(".", 1)
        n = 0
        hdr, seq = None, []
        def flush():
            nonlocal n
            if hdr and seq:
                s = "".join(seq)
                if s not in bins[(module, key)]:
                    bins[(module, key)][s] = hdr
                    n += 1
        for line in open(fa, errors="replace"):
            if line.startswith(";") or not line.strip():
                continue
            if line.startswith(">"):
                flush(); hdr, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
        flush()
        if n:
            added.append((module, key, n))
    return added


def split_unchar(bins):
    """Replace each module's UNCHAR pool with homology groups UNC1..UNCn."""
    import hashlib
    out = []
    for (module, key) in [k for k in list(bins) if k[1] == "UNCHAR"]:
        pool = bins[(module, key)]
        if not pool:
            del bins[(module, key)]
            continue
        if not shutil.which("mmseqs"):
            print("\nWARNING: mmseqs not on PATH -- %s/UNCHAR left as one pool"
                  % module)
            continue
        #  Drop pool members that are already a named feature BEFORE grouping.
        #  Order matters: filtering afterwards empties groups that were built
        #  around those sequences and leaves gaps in the numbering, which is
        #  how UNC5, UNC8 and UNC9 came out with zero members on the first
        #  attempt. Same 80% line as adoption.
        named = {(m, k): d for (m, k), d in bins.items()
                 if m == module and not k.startswith("UNC")}
        pool = _drop_named_duplicates(module, pool, named)
        if len(pool) < UNC_MIN:
            del bins[(module, key)]
            print("%s: UNCHAR pool reduced to %d after removing proteins already "
                  "covered by a named feature; no groups" % (module, len(pool)))
            continue
        items = sorted(pool.items(), key=lambda kv: kv[1])
        tmp = tempfile.mkdtemp(prefix="unc.")
        try:
            fa = os.path.join(tmp, "u.faa")
            with open(fa, "w") as fh:
                for i, (s, _h) in enumerate(items):
                    fh.write(">u%d\n%s\n" % (i, s))
            subprocess.run(
                ["mmseqs", "easy-cluster", fa, os.path.join(tmp, "c"),
                 os.path.join(tmp, "t"), "--min-seq-id", str(UNC_PIDENT),
                 "-c", str(UNC_COV), "--cov-mode", "0"],
                capture_output=True)
            groups = collections.defaultdict(list)
            tsv = os.path.join(tmp, "c_cluster.tsv")
            if os.path.exists(tsv):
                for line in open(tsv):
                    rep, mem = line.split()[:2]
                    groups[rep].append(int(mem[1:]))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        if not groups:
            continue
        del bins[(module, key)]
        ordered = sorted(groups.values(),
                         key=lambda g: (-len(g),
                                        hashlib.md5(items[min(g)][0].encode()).hexdigest()))
        n = 0
        kept = dropped = 0
        for g in ordered:
            if len(g) < UNC_MIN:
                dropped += len(g)
                continue
            n += 1
            k = "UNC%d" % n
            for i in g:
                s, h = items[i]
                bins[(module, k)][s] = h
            kept += len(g)
            out.append((module, k, len(g)))
        print("%s: UNCHAR pool of %d split into %d homology group(s) "
              "(%d sequences kept, %d in groups below %d and left out)"
              % (module, len(items), n, kept, dropped, UNC_MIN))
    return out

# --------------------------------------------------------------------------
# Homology rescue for proteins the source never named.
#
# Every binning rule above reads annotation text, so a protein BV-BRC left
# blank, or labelled positionally as "3 protein", matches nothing and is lost.
# That is not a rare edge: 158 proteins in this dump carry a blank or
# positional annotation, and 150 of them are recognisably a named feature --
# 145 at 90% identity or better. Eighty of those are glycoproteins.
#
# So the rescue is by homology, not by text. Blast each unnamed protein against
# the named collections OF ITS OWN MODULE and adopt the best hit.
#
# One threshold, and what happens on each side of it, per JJD:
#
#   >= RESCUE_PIDENT over >= RESCUE_COV of the query   adopted into the feature.
#        At 80% identity across 60% of the length within one module, the
#        protein is that feature.
#
#   < RESCUE_PIDENT                                    LEFT OFF entirely.
#        Not adopted, and not swept into UNCHAR either. Below 80% it is its
#        own lineage-specific protein, and without a literature-based function
#        to say what that protein is, there is nothing to build. Rice yellow
#        stunt gene 3 is here at 69% to a movement protein: recognisably
#        movement-protein-like and not the same protein, so calling it MOV
#        would assert a homology that 69% does not support, and filing it
#        under "uncharacterized" would only dilute that bag with proteins we
#        declined to identify. They are listed in RESCUE_REVIEW.tsv so the
#        decision is visible, and left out of the collections.
#
# Note what this rules OUT. A bad annotation -- "P6", "hypothetical protein" --
# must never be gathered by regex, because the string is a position or a
# shrug, not a protein. Grouping on it manufactures a feature out of proteins
# that have nothing in common. If those proteins are to be used at all, only
# homology may group them, which is exactly what UNCHAR plus MMseqs does.
#
# Nor does every protein need a home. A collection is training data, not an
# inventory: the annotator blasts each finished profile back against the
# incoming genome, so a protein left out here can still be called by the
# profile its relatives built.
#
# Routing is to the feature key, never straight to the fasta, so EXPECTED
# still applies afterwards and a rescued protein of the wrong length goes to
# <FEATURE>.outliers.fasta exactly as a correctly-named one would.
# --------------------------------------------------------------------------
RESCUE_PIDENT = 80.0
RESCUE_COV = 0.60

def homology_rescue(bins, unnamed):
    """Adopt unnamed proteins into the named feature they are homologous to.

    Returns (adopted, review). Mutates bins. Needs blastp; without it the
    collections are still valid, just missing the rescues, so this warns
    rather than failing.
    """
    adopted, lineage = [], []
    if not unnamed:
        return adopted, lineage
    if not shutil.which("blastp"):
        print("\nWARNING: blastp not on PATH -- %d unnamed protein(s) not "
              "rescued. Collections are usable but incomplete." % len(unnamed))
        return adopted, lineage

    by_mod = collections.defaultdict(list)
    for rec in unnamed:
        by_mod[rec[0]].append(rec)

    for module, recs in sorted(by_mod.items()):
        named = {(m, k): d for (m, k), d in bins.items()
                 if m == module and k != "UNCHAR"}
        if not named:
            continue
        tmp = tempfile.mkdtemp(prefix="rescue.")
        try:
            db = os.path.join(tmp, "named.faa")
            with open(db, "w") as fh:
                for (m, k), d in sorted(named.items()):
                    for i, (s, _h) in enumerate(d.items()):
                        fh.write(">%s@%d\n%s\n" % (k, i, s))
            subprocess.run(["makeblastdb", "-in", db, "-dbtype", "prot",
                            "-out", os.path.join(tmp, "db")], capture_output=True)
            q = os.path.join(tmp, "q.faa")
            with open(q, "w") as fh:
                for i, rec in enumerate(recs):
                    fh.write(">q%d\n%s\n" % (i, rec[1]))
            r = subprocess.run(
                ["blastp", "-query", q, "-db", os.path.join(tmp, "db"), "-outfmt",
                 "6 qseqid sseqid pident length qlen bitscore", "-evalue", "1e-5",
                 "-max_target_seqs", "5", "-num_threads", "4"],
                capture_output=True, text=True)
            best = {}
            for line in r.stdout.splitlines():
                qid, sid, pid, alen, qlen, bits = line.split("\t")
                i = int(qid[1:])
                cand = (float(pid), int(alen) / max(int(qlen), 1), float(bits),
                        sid.split("@")[0])
                if i not in best or cand[2] > best[i][2]:
                    best[i] = cand
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

        for i, rec in enumerate(recs):
            if i not in best:
                continue
            pid, cov, _bits, key = best[i]
            _mod, s, fid, name, gen, ann, _gid = rec
            if pid >= RESCUE_PIDENT and cov >= RESCUE_COV:
                hdr = "%s [%s] %s | %s" % (fid, name, gen,
                                           ann or "(unannotated, rescued to %s)" % key)
                bins[(module, key)].setdefault(s, hdr)
                adopted.append((module, key, fid, name, gen, len(s), pid, cov, ann))
            else:
                #  Related but not the same protein, and nothing in the
                #  literature says what it is. Left out of the collections
                #  entirely -- recorded, not used.
                lineage.append((module, key, fid, name, gen, len(s), pid, cov, ann))
    return adopted, lineage

EXPECTED = {
    "N": (330, 560), "P": (190, 360), "M": (150, 380),
    "G": (400, 700), "L": (1800, 2300), "MOV": (250, 380),
    "NV": (95, 140),
    "EPHEM_GNS": (500, 640), "EPHEM_ALPHA1": (75, 125),
    "EPHEM_ALPHA2": (85, 140), "EPHEM_BETA": (100, 170),
    "EPHEM_GAMMA": (90, 130), "EPHEM_DELTA": (95, 120),
    "SIGMA_X": (215, 360), "VESI_C": (45, 75), "VESI_CPRIME": (55, 80),
}

def classify(genus, ann):
    """Map an annotation string to a feature key.

    Matched TWICE, deliberately in this order:

      1. the annotation exactly as BV-BRC supplies it
      2. a normalised form, with hyphens and underscores turned into spaces
         and runs of whitespace collapsed

    The second pass is strictly additive -- anything that matched before still
    matches on pass 1, so no existing binning changes. It exists because the
    rules were written against the spellings present in one dump and BV-BRC
    punctuates the same product string inconsistently: "L-protein" missed
    `^L (protein|polymerase)` on a single hyphen and sat in UNBINNED, which is
    how Varicosavirus lost a third of its L sequences. Family-wide this
    recovers L-protein (Varicosavirus), L-polymerase (Hapavirus),
    RNA-dependent-RNA-polymerase (Alphanucleorhabdovirus), RNA-dependent
    polymerase (Dichorhavirus) and RNA-dependent RNA polyermaser
    (Betanucleorhabdovirus).
    """
    if not ann:
        return None
    rules = COMPILED.get(genus, CORE_ONLY)
    for key, prx in rules:
        if prx.search(ann):
            return key
    norm = re.sub(r"[-_]", " ", ann)
    norm = re.sub(r"\s+", " ", norm).strip()
    if norm != ann:
        for key, prx in rules:
            if prx.search(norm):
                return key
    return None

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    os.chdir(here)

    # genome_id -> (name, genus)
    genome = {}
    with open("Rhabdo.id_name_gs") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                # (name, genus, family) -- family is needed to say what
                # taxonomic group an unassigned protein came from.
                genome[p[0]] = (p[1], p[3].strip(), p[2].strip())
    print(f"genomes: {len(genome)}")

    # md5 -> sequence  (uniq.seq)
    seq_of = {}
    with open("Rhabdo.uniq.seq") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                seq_of[p[0]] = p[1].strip().upper().rstrip("*")
    print(f"unique sequences: {len(seq_of)}")

    # md5 -> representative feature id + product.
    #
    # Joined BY KEY, not by line number. The three uniq.* files used to be
    # positionally aligned and the loader relied on that, but the alignment is
    # a side effect of SOLR returning rows in input order and it does not
    # survive a re-export: uniq.seq comes back md5-sorted while uniq.id_ann is
    # in whatever order the feature query produced, and a few md5s can be
    # missing a sequence entirely. Both files carry an explicit key, so use it.
    id2md5 = {}
    with open("Rhabdo.id_md5") as fh:
        for line in fh:
            p_ = line.rstrip("\n").split("\t")
            if len(p_) >= 2 and p_[1].strip():
                id2md5[p_[0]] = p_[1].strip()
    ann_of = {}
    rep_id_of = {}
    for line in open("Rhabdo.uniq.id_ann"):
        row = line.rstrip("\n").split("\t")
        if not row or not row[0].strip():
            continue
        m = id2md5.get(row[0].strip())
        if not m:
            continue
        rep_id_of[m] = row[0].strip()
        ann_of[m] = row[1].strip() if len(row) > 1 else ""
    # the working set: every unique sequence we have both a product and a
    # sequence for. Replaces the old positional `md5s` list.
    md5s = [m for m in ann_of if m in seq_of]
    _no_seq = [m for m in ann_of if m not in seq_of]
    print(f"annotated unique sequences: {len(ann_of)}"
          + (f"  ({len(_no_seq)} with no sequence, skipped)" if _no_seq else ""))

    # md5 -> set of genome_ids that carry it; and genome -> set of md5s
    md5_genomes = collections.defaultdict(set)
    genome_md5s = collections.defaultdict(set)
    n_feat = 0
    with open("Rhabdo.id_md5") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or not p[1]:
                continue
            fid, m = p[0], p[1]
            n_feat += 1
            #  BUG FIX. This used to strip only a ".CDS.<n>" suffix and
            #  otherwise take the whole feature id as the genome id. BV-BRC
            #  also uses .peg., .mRNA., .rna., .misc_feature., .gap. and
            #  others, so every non-CDS feature got a genome id like
            #  "1034378.4.peg.2" that is absent from Rhabdo.id_name_gs. The
            #  genus lookup then failed and the feature was treated as having
            #  no genus. That affected 36389 .peg. features alone, whose real
            #  genomes resolve to Lyssavirus (21517), Vesiculovirus (2755)
            #  and so on. It also broke the core-protein tally, because each
            #  feature of such a genome got its OWN distinct bogus id, so
            #  found[gid] could never accumulate five keys for it.
            #  A BV-BRC genome id is always <taxon_id>.<version>, so anchor
            #  on that rather than trying to enumerate feature-type suffixes
            #  (.CDS. .peg. .mRNA. .rna. .misc_feature. .gap. .3'UTR. ...).
            _body = fid.split("|", 1)[-1]
            _m = re.match(r"(\d+\.\d+)(?:\.|$)", _body)
            gid = _m.group(1) if _m else _body
            md5_genomes[m].add(gid)
            genome_md5s[gid].add(m)
    print(f"features: {n_feat}")

    # ---- synonyms census: annotation x genus -------------------------------
    syn = collections.Counter()
    for m in md5s:
        a = ann_of.get(m, "")
        for gid in md5_genomes.get(m, ()):
            gen = genome.get(gid, ("", "", ""))[1] or "(no genus)"
            syn[(a, gen)] += 1
    with open("synonyms.tsv", "w") as fh:
        fh.write("count\tannotation\tgenus\tbins_to\n")
        for (a, gen), c in syn.most_common():
            key = classify(gen.lower(), a) or ""
            fh.write(f"{c}\t{a}\t{gen}\t{key}\n")
    print(f"synonyms.tsv: {len(syn)} annotation x genus pairs")

    # ---- bin the unique sequences -----------------------------------------
    bins = collections.defaultdict(dict)
    misbinned_hits = []
    l_rescued = []
    unnamed = []
    unbinned = collections.Counter()
    dropped = collections.Counter()
    no_ann = []
    unassigned = []   # features from genomes with no genus in BV-BRC
    not_modelled = [] # features from genera deliberately not modelled
    trimmed = []      # N-terminal curation overrides actually applied
    # genome -> set of core features found, for the MISSING report
    found = collections.defaultdict(set)

    for m in md5s:
        s = seq_of.get(m, "")
        _cut = NTERM_TRIM.get(rep_id_of.get(m, ""))
        if _cut:
            s = s[_cut:]
            trimmed.append((rep_id_of.get(m, ""), _cut, len(s)))
        a = ann_of.get(m, "")
        if not s:
            dropped["no sequence"] += 1
            continue
        if len(s) < MIN_AA:
            dropped[f"shorter than {MIN_AA} aa"] += 1
            continue
        if set(s) & AMBIG:
            dropped["ambiguous residues"] += 1
            continue

        for gid in md5_genomes.get(m, ()):
            name, gen, fam = genome.get(gid, ("", "", ""))
            gl = gen.lower()

            #  Genomes with no genus in BV-BRC no longer form a collection of
            #  their own. They used to go to an UNASSIGNED_no_genus module,
            #  but that module is not in the JSON, and most of what landed in
            #  it was a second copy of a sequence already training a real
            #  module (the loop below emits one record per genome carrying an
            #  md5, so a sequence shared with a named-genus genome appears in
            #  both). Per JJD: duplicates are fine to drop, because the shared
            #  md5 is still covered by the real module at annotation time.
            #  What must not be lost is the record of the ones that are NOT
            #  covered anywhere, and which taxonomic group they came from --
            #  those are tracked here and written to UNASSIGNED_TRACKING.tsv.
            if not gen:
                unassigned.append((classify(gl, a) or "(unbinned)", len(s),
                                   rep_id_of.get(m, ""), a, fam, name, gid, s))
                continue

            if gl in NOT_MODELLED:
                not_modelled.append((gen, classify(gl, a) or "(unbinned)", len(s),
                                     rep_id_of.get(m, ""), a, name, gid,
                                     NOT_MODELLED[gl]))
                continue

            module = GENUS_MODULE.get(gl)
            if module is None:
                dropped[f"genus not routed: {gen}"] += 1
                continue
            key = classify(gl, a)
            _ov = MISBINNED.get(rep_id_of.get(m, ""))
            if _ov is not None:
                key = _ov[0]
                misbinned_hits.append((rep_id_of.get(m, ""), a, key, _ov[1]))
                if key is None:
                    continue
            if key == "UNCHAR" and len(s) >= UNCHAR_L_RESCUE:
                l_rescued.append((rep_id_of.get(m, ""), name, gen, len(s)))
                key = "L"
            if not key:
                unbinned[(a, gen)] += 1
                if not a:
                    no_ann.append((rep_id_of.get(m, ""), name, gen, len(s)))
                #  no rule named it -- homology gets a turn after the loop
                unnamed.append((module, s, rep_id_of.get(m, ""), name, gen, a, gid))
                continue
            hdr = f"{rep_id_of.get(m,'')} [{name}] {gen} | {a}"
            bins[(module, key)].setdefault(s, hdr)
            if key in CORE_KEYS:
                found[gid].add(key)

    # ---- homology rescue, before anything is written ------------------------
    _adopted, _review = homology_rescue(bins, unnamed)   # _review = lineage-specific
    for _m, _k, _f, _n, _g, _l, _p, _c, _a in sorted(_adopted, key=lambda r: (r[0], r[1])):
        print("rescued by homology: %-22s -> %-10s %5d aa  %.1f%% id  %s [%s]"
              % (_f, "%s/%s" % (_m, _k), _l, _p, _n[:30], _a or "unannotated"))
    if _adopted:
        print("%d unnamed protein(s) adopted at >=%.0f%% identity over >=%.0f%% "
              "of length" % (len(_adopted), RESCUE_PIDENT, 100 * RESCUE_COV))
    if _review:
        print("%d unnamed protein(s) below %.0f%% identity LEFT OFF "
              "(lineage-specific, no literature function); listed in "
              "collections/RESCUE_REVIEW.tsv" % (len(_review), RESCUE_PIDENT))

    # ---- merge proteins recovered directly from genomes ---------------------
    for _m, _k, _n in merge_extracted(bins, here):
        print("merged %d recovered protein(s) into %s/%s from Extracted/" % (_n, _m, _k))

    # ---- re-home length outliers the text rules misplaced -------------------
    _moved = reroute_length_outliers(bins)
    for _m, _from, _to, _L, _p, _f in sorted(_moved):
        print("re-homed by length+homology: %-22s %-6s -> %-6s %5d aa  %.1f%% id  %s"
              % (_m, _from, _to, _L, _p, _f))
    if _moved:
        print("%d protein(s) moved to the feature their length and homology agree on"
              % len(_moved))

    # ---- split the uncharacterized pool into homology groups ---------------
    _unc = split_unchar(bins)
    if _unc:
        print("%d uncharacterized homology group(s) across %d module(s)"
              % (len(_unc), len({u[0] for u in _unc})))

    # ---- write collections ------------------------------------------------
    out = "collections"
    os.makedirs(out, exist_ok=True)
    #  Clear EVERY collection fasta from a previous run before writing. Not
    #  just the UNC groups: retiring a binning rule leaves its collection on
    #  disk otherwise, and BCYTO_P6.fasta, TUPA_SH.fasta and the rest sat
    #  there with 22, 9 and 6 sequences after their rules were removed. A
    #  collection nothing produces any more is the same defect as a profile
    #  nothing produces any more -- it reads as current and is not.
    for _d in glob.glob(os.path.join(out, "*")):
        if not os.path.isdir(_d):
            continue
        for _f in glob.glob(os.path.join(_d, "*.fasta")):
            os.remove(_f)

    def write_fa(path, items):
        with open(path, "w") as fh:
            for s, hdr in sorted(items, key=lambda kv: kv[1]):
                fh.write(f">{hdr}\n")
                for i in range(0, len(s), 60):
                    fh.write(s[i:i + 60] + "\n")

    # ---- keep UNCHAR disjoint from the named features ---------------------
    # UNCHAR is a grab-bag of proteins no rule named. If one of them is also
    # near-identical to a member of a named collection in the same module, both
    # profiles fire on the same locus and the annotator emits two features on
    # identical coordinates -- always a defect. Rabies survives this because
    # its UNCHAR members happen to be distinct; Novirhabdovirus P and
    # Dichorhavirus N did not, and each produced a duplicate call.
    #
    # The grab-bag is the side that yields: a protein covered by a named
    # feature does not belong in it.
    _named = collections.defaultdict(set)
    for (module, key), d in bins.items():
        if key != "UNCHAR":
            _named[module] |= set(d)
    _dropped_unchar = 0
    for (module, key) in list(bins):
        if key != "UNCHAR":
            continue
        keep = {}
        for s, hdr in bins[(module, key)].items():
            if s in _named[module]:
                _dropped_unchar += 1
                continue
            keep[s] = hdr
        bins[(module, key)] = keep
    if _dropped_unchar:
        print(f"UNCHAR: dropped {_dropped_unchar} sequence(s) already covered "
              f"by a named feature in the same module")

    rows = []
    for (module, key), d in sorted(bins.items()):
        dd = os.path.join(out, module)
        os.makedirs(dd, exist_ok=True)
        lo, hi = EXPECTED.get(key, (0, 10 ** 6))
        keep = [(s, h) for s, h in d.items() if lo <= len(s) <= hi]
        out_l = [(s, h) for s, h in d.items() if len(s) < lo]
        out_h = [(s, h) for s, h in d.items() if len(s) > hi]
        write_fa(os.path.join(dd, f"{key}.fasta"), keep)
        if out_l or out_h:
            write_fa(os.path.join(dd, f"{key}.outliers.fasta"), out_l + out_h)
        lens = sorted(len(s) for s, _h in keep) or [0]
        rows.append((module, key, len(keep), len(out_l), len(out_h),
                     lens[0], lens[-1], lens[len(lens) // 2],
                     "PSSM-ok" if len(keep) >= 3 else "too few"))

    with open(os.path.join(out, "RESCUE_REVIEW.tsv"), "w") as fh:
        fh.write("# Below %.0f%% identity to their nearest named feature, so NOT\n"
                 "# adopted into it and NOT added to any collection. Each is a\n"
                 "# lineage-specific protein with no literature-based function to\n"
                 "# identify it. Recorded here so the omission is visible; the\n"
                 "# annotator may still call some of them when a related profile\n"
                 "# is blasted back against the genome.\n" % RESCUE_PIDENT)
        fh.write("module\tnearest_feature\tfeature_id\tgenome\tgenus\taa\t"
                 "pident\tqcov\tsource_annotation\n")
        for r in sorted(_review, key=lambda x: (-x[6], x[0])):
            fh.write("%s\t%s\t%s\t%s\t%s\t%d\t%.1f\t%.2f\t%s\n"
                     % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8]))

    with open(os.path.join(out, "report.tsv"), "w") as fh:
        fh.write("module\tfeature\tn_in_range\tn_short\tn_long"
                 "\tmin_aa\tmax_aa\tmedian_aa\tstatus\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    #  Which real module/feature collections contain each sequence, so an
    #  unassigned protein can be marked as already covered or genuinely lost.
    seq_modules = collections.defaultdict(set)
    for (module, key), d in bins.items():
        for s_ in d:
            seq_modules[s_].add(f"{module}/{key}")

    with open(os.path.join(out, "UNASSIGNED_TRACKING.tsv"), "w") as fh:
        fh.write("feature\tlength\trep_id\tannotation\tfamily\tgenome_name"
                 "\tgenome_id\tcovered_by\n")
        n_lost = 0
        for key, ln, rid, a, fam, name, gid, s_ in sorted(
                unassigned, key=lambda r: (r[0], -r[1], r[5])):
            cov = ",".join(sorted(seq_modules.get(s_, ())))
            if not cov:
                n_lost += 1
            fh.write(f"{key}\t{ln}\t{rid}\t{a}\t{fam}\t{name}\t{gid}\t{cov}\n")
    print(f"UNASSIGNED_TRACKING.tsv: {len(unassigned)} features from genomes "
          f"with no genus; {n_lost} not covered by any real module")

    if trimmed:
        print(f"N-terminal trims applied: {len(trimmed)}")
        for fid, cut, newlen in sorted(trimmed):
            print(f"  {fid}  -{cut} aa  -> {newlen} aa")
    for _fid, _ann, _key, _why in sorted(set(misbinned_hits)):
        print(f"MISBINNED override: {_fid} '{_ann}' -> {_key}  ({_why})")
    for _fid, _nm, _gen, _L in sorted(l_rescued, key=lambda r: -r[3]):
        print(f"L rescued by length: {_fid} {_L} aa  [{_nm}] {_gen}")
    if l_rescued:
        print(f"{len(l_rescued)} unnamed protein(s) >= {UNCHAR_L_RESCUE} aa routed to L")
    _mb_missed = sorted(set(MISBINNED) - {h[0] for h in misbinned_hits})
    if _mb_missed:
        print("WARNING: MISBINNED entries that matched no sequence (stale ids?): "
              + ", ".join(_mb_missed))

    missed = sorted(set(NTERM_TRIM) - {t[0] for t in trimmed})
    if missed:
        print("WARNING: NTERM_TRIM entries that matched no sequence (stale ids?): "
              + ", ".join(missed))

    with open(os.path.join(out, "NOT_MODELLED.tsv"), "w") as fh:
        fh.write("genus\tfeature\tlength\trep_id\tannotation\tgenome_name"
                 "\tgenome_id\treason\n")
        for row in sorted(not_modelled, key=lambda r: (r[0], r[1], -r[2])):
            fh.write("\t".join(str(x) for x in row) + "\n")
    if not_modelled:
        _g = sorted({r[0] for r in not_modelled})
        print(f"NOT_MODELLED.tsv: {len(not_modelled)} features from "
              f"{len(_g)} deliberately unmodelled genus/genera ({', '.join(_g)})")

    with open(os.path.join(out, "UNBINNED.tsv"), "w") as fh:
        fh.write("count\tannotation\tgenus\n")
        for (a, gen), c in unbinned.most_common():
            fh.write(f"{c}\t{a}\t{gen}\n")

    # Proteins with NO annotation string at all. These cannot be binned by
    # name; they need assignment by length + genomic position, which means
    # running the annotator, not a regex. Listed with lengths so the obvious
    # ones (a 2100 aa protein is L) can be swept manually.
    with open(os.path.join(out, "NO_ANNOTATION.tsv"), "w") as fh:
        fh.write("# Unique sequences carrying an empty product string.\n")
        fh.write("rep_feature_id\tgenome_name\tgenus\taa_length"
                 "\tlikely_by_length\n")
        for fid, name, gen, ln in sorted(no_ann, key=lambda r: -r[3]):
            guess = ("L" if ln >= 1700 else "G" if 400 <= ln <= 700
                     else "N" if 330 <= ln < 400 else "M" if 150 <= ln < 330
                     else "P" if 190 <= ln < 360 else "")
            fh.write(f"{fid}\t{name}\t{gen}\t{ln}\t{guess}\n")

    # ---- MISSING report: genomes lacking each core protein ---------------
    # Only meaningful for genomes that have at least one called feature.
    # Restricted to NEAR-COMPLETE genomes (>=4 of the 5 core proteins found).
    # Without that restriction the report is useless: this dump is mostly
    # single-gene submissions -- 34,777 of the 48,074 "genomes" are
    # Lyssavirus records and most carry one gene -- so almost everything is
    # "missing" something. A genome with 4 of 5 is the interesting case: it
    # is a real genome with one gene the rules failed to call.
    with open(os.path.join(out, "MISSING.tsv"), "w") as fh:
        fh.write("# Near-complete genomes (>=4 of 5 core proteins called)\n"
                 "# that are still missing one. These are the cases worth\n"
                 "# looking at -- either a missing PSSM or a bad annotation.\n")
        fh.write("genome_id\tgenome_name\tgenus\tmissing\tn_found\n")
        n_missing = 0
        for gid, keys in sorted(found.items()):
            miss = sorted(CORE_KEYS - keys)
            if miss and len(keys) >= 4:
                name, gen, _fam = genome.get(gid, ("", "", ""))
                fh.write(f"{gid}\t{name}\t{gen}\t{','.join(miss)}\t"
                         f"{len(keys)}\n")
                n_missing += 1

    complete = sum(1 for k in found.values() if len(k) == 5)
    print(f"\ncollections: {len(rows)}")
    print(f"in-range sequences:  {sum(r[2] for r in rows)}")
    print(f"length outliers:     {sum(r[3] + r[4] for r in rows)} "
          f"({sum(r[3] for r in rows)} short, {sum(r[4] for r in rows)} long)")
    print(f"unbinned: {sum(unbinned.values())} occurrences, "
          f"{len(unbinned)} distinct annotation x genus")
    print(f"genomes with all 5 core proteins: {complete}")
    print(f"near-complete genomes missing exactly one: {n_missing}")
    thin = [r for r in rows if r[8] == "too few"]
    if thin:
        print(f"\n{len(thin)} collections have <3 in-range sequences and "
              f"cannot support a PSSM:")
        for r in thin:
            print(f"  {r[0]}/{r[1]}: {r[2]}")
    for reason, n in dropped.most_common(6):
        print(f"  dropped, {reason}: {n}")

def _run_unchar_filter():
    """Re-apply the blastp UNCHAR filter, because regenerating undoes it.

    The disjointness step above only catches exact sequence matches. The real
    overlap is near-identical proteins from different genomes, which needs
    blastp -- filter_unchar.py. That script rewrites UNCHAR.fasta in place, so
    its result lives only in a file, and every run of this script silently
    reverts it: Alpharhabdovirinae UNCHAR went 453 -> 515 mid-rebuild with
    nothing to indicate the named-feature duplicates were back.

    A curation step that has to be remembered is a curation step that will be
    forgotten. Invoke it here so regenerating the collections cannot leave them
    half-curated.
    """
    import shutil
    import subprocess
    script = os.path.expanduser(
        "~/.claude/skills/lowvan-module/scripts/filter_unchar.py")
    if not os.path.exists(script):
        print("\nWARNING: filter_unchar.py not found -- UNCHAR still contains "
              "proteins already covered by named features. Run it before "
              "building PSSMs.")
        return
    if not shutil.which("blastp"):
        print("\nWARNING: blastp not on PATH -- skipping the UNCHAR filter. "
              "UNCHAR still overlaps the named features.")
        return
    print("\n--- filter_unchar.py (blastp) ---")
    subprocess.run([sys.executable, script, "--workdir", os.getcwd(), "--write"])

if __name__ == "__main__":
    main()
    _run_unchar_filter()
