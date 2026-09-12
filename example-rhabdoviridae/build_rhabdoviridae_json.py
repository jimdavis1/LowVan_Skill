#!/usr/bin/env python3
"""
Generate the Rhabdoviridae taxon blocks for LowVan's Viral_PSSM.json.

REVISION 2. Three changes from the first draft, all per JJD:

 1. Taxa with fewer than 5 sequenced genomes in BV-BRC are not modelled at all.
    Counts were measured from the BV-BRC taxonomy collection (the `genomes`
    field) -- see GENOME_COUNTS below, and the notes for which are measured
    vs assumed.

 2. Single-segment genus-level groups are COLLAPSED to subfamily-level
    modules. Minor proteins from the wrong genus simply fail to match when
    the bit_cutoff is set correctly -- the same arrangement already used in
    Betacoronavirus, where nobeco-specific PSSMs are blasted against a
    SARS-CoV-2 genome and find nothing. Only the multi-segmented genera stay
    separate, because they need their own `segments` schema.
    34 modules -> 6.

 3. Nomenclature follows S1 Table. In particular: taxon-specific information
    is NOT purged. It is carried in the gene symbol as a lineage prefix, the
    way S1 does it for Sarbeco_ORF6 / Nobeco_NS7a / Merbeco_ORF4b /
    G_COV_ORF4a / D_COV_NS7b. Collapsing to subfamily modules makes this
    mandatory rather than merely tidy: Ephemero_U1, Tibro_U1 and Hapa_U1 are
    three different proteins that now live in one module.

Other S1 conventions adopted here that the first draft got wrong:
  * "Uncharacterized lineage-specific <NAME> protein", or
    "Uncharacterized lineage-specific protein" with an EMPTY symbol when the
    protein has no name at all (cf. Ferlavirus, Fimoviridae, Gammacoronavirus).
  * Signal peptides ARE annotated: "Signal peptide of <parent>" / SP /
    mat_peptide (cf. Orthoebolavirus, Orthopneumovirus, Arenaviridae).
    So rhabdovirus G gets one, and a mature chain.
  * "C-prime protein" / "C-Prime", not "C' protein" / "C'".
  * Plant movement protein is "Movement protein" / "Mov", matching Fimoviridae.
  * Segments are named "Segment 1", "Segment 2", ... matching Fimoviridae.
  * Descriptor precedes the name: "Helicase NSP13 protein" style.

Run:  python3 build_rhabdoviridae_json.py > Rhabdoviridae_Viral_PSSM.json
"""

import glob
import json
import math
import os
import re
import sys

SINGLE = "Single RNA Segment"

# --------------------------------------------------------------------------
# BV-BRC genome counts, measured from
#   /api/taxonomy/?lineage_ids=<subfamily>&taxon_rank=genus
# on 2026-08-26. Genera below the cutoff are dropped.
# --------------------------------------------------------------------------
MIN_GENOMES = 5

GENOME_COUNTS = {
    # --- Alpharhabdovirinae, measured
    "Alphapaprhavirus": 269, "Perhabdovirus": 200, "Ledantevirus": 106,
    "Sunrhavirus": 57, "Almendravirus": 49, "Tupavirus": 34, "Tibrovirus": 31,
    "Sripuvirus": 26, "Curiovirus": 15, "Arurhavirus": 12,
    "Alphanemrhavirus": 10, "Caligrhavirus": 10, "Sawgrhavirus": 10,
    "Lostrhavirus": 9, "Barhavirus": 8, "Alphathriprhavirus": 7,
    "Amplylivirus": 5,
    "Scophrhavirus": 4, "Cetarhavirus": 4, "Mousrhavirus": 3,
    "Zarhavirus": 2, "Replylivirus": 2, "Uniorhavirus": 2,
    "Betathriprhavirus": 2,
    # --- Gammarhabdovirinae, measured
    "Novirhabdovirus": 3090, "Margarhavirus": 2,
    # --- Betarhabdovirinae, measured
    "Deltanucleorhabdovirus": 866, "Alphanucleorhabdovirus": 356,
    "Alphacytorhabdovirus": 313, "Varicosavirus": 262, "Dichorhavirus": 226,
    "Betacytorhabdovirus": 216, "Betanucleorhabdovirus": 129,
    "Gammacytorhabdovirus": 21, "Alphagymnorhavirus": 14, "Trirhavirus": 12,
    "Gammanucleorhabdovirus": 8, "Betagymnorhavirus": 2,
    # --- Deltarhabdovirinae, measured (all 11)
    "Primrhavirus": 112, "Stangrhavirus": 50, "Betaricinrhavirus": 14,
    "Betapaprhavirus": 13, "Betanemrhavirus": 8, "Alphahymrhavirus": 7,
    "Gammahymrhavirus": 6, "Gammaricinrhavirus": 6,
    "Alphacrustrhavirus": 4, "Alphadrosrhavirus": 4, "Betahymrhavirus": 3,
    # --- ASSUMED >= 5, not measured (the 25-row API cap cut the listing off
    #     before these). All are long-established, heavily sequenced genera;
    #     Lyssavirus and Vesiculovirus alone run to thousands of genomes.
    "Lyssavirus": None, "Vesiculovirus": None, "Ephemerovirus": None,
    "Sigmavirus": None, "Hapavirus": None, "Sprivivirus": None,
    "Merhavirus": None, "Ohlsrhavirus": None,
}

DROPPED = sorted(g for g, n in GENOME_COUNTS.items()
                 if n is not None and n < MIN_GENOMES)


def bounds(lo, hi):
    return int(math.floor(lo * 0.9)), int(math.ceil(hi * 1.1))


def bit_for(min_len):
    """Scale the tBLASTn cutoff to protein size.

    This matters more now that the modules are collapsed: the bit_cutoff is
    what stops an Ephemerovirus-specific PSSM from calling something in a
    lyssavirus genome. Small accessory proteins cannot reach 100, so they get
    proportionally lower cutoffs -- and those are the entries most in need of
    empirical tuning once the module is run over real genomes.
    """
    if min_len >= 300:
        return 100
    if min_len >= 180:
        return 80
    if min_len >= 120:
        return 60
    if min_len >= 80:
        return 40
    return 30


# --------------------------------------------------------------------------
# DESIGNATION PROVENANCE
# --------------------------------------------------------------------------
# A "designation" is a short positional tag inside an annotation string: the U1
# in "Uncharacterized lineage-specific U1 protein", the Gx, the Px, the alpha1.
#
# A designation goes into the annotation string ONLY when the source data
# already uses it for that genus. Everything below was verified against
# synonyms.tsv: BV-BRC calls Hapavirus's accessory ORFs U1/U2/U3, so HAPA_U1
# keeps its U1; nothing calls Arurhavirus's ORFs U1 or U2, so ARURHA_U1 does
# not get one.
#
# The feature KEY and the gene_symbol still carry the tag -- they are internal
# identifiers and have to be unique. It is the published annotation string that
# must not invent a designation, because a reader has no way to tell an
# invented U1 from a real one, and the number implies a correspondence between
# genera that does not exist. The cost is that many features share the string
# "Uncharacterized lineage-specific protein". That is correct: they are
# uncharacterized lineage-specific proteins, and the Taxon column says which
# lineage.
#
# A designation earns its way onto this list by appearing in the source data or
# in a DLIT that names it. Add to it when that evidence turns up.
BACKED_DESIGNATION = {
    # BV-BRC uses U1/U2/U3 for these genera (Walker 2015, PMID 25679389)
    "HAPA_U1", "HAPA_U2", "HAPA_U3",
    # ephemerovirus alpha/beta/gamma/delta ORFs and GNS are named in the data
    "EPHEM_ALPHA1", "EPHEM_ALPHA2", "EPHEM_ALPHA3", "EPHEM_BETA",
    "EPHEM_GAMMA",
    "EPHEM_GNS",
    # plant rhabdovirus P-numbers that the submitters use "BCYTO_P4", "BCYTO_P5",
    #  Keys retired 11 Sep 2026 (too few sequences to build) were removed
    #  from this set with their declarations: ANUCLEO_P6, BCYTO_P6, EPHEM_DELTA, SRIPU_GX, SRIPU_MX, SRIPU_U1, TIBRO_U1, TIBRO_U2, TIBRO_U3.
    #  A designation is only "backed" while the feature it names exists.

}

_DESIG = re.compile(
    r"\s+\b(U\d?|G[xy]|[NMP]x|P\d|alpha\d?|beta|gamma|delta)\b(?=\s+protein\b)")


UNCHAR_ANNO = "Uncharacterized lineage-specific protein"

# The family-wide controlled vocabulary: strings that are either already in
# S1 Table or describe a protein whose identity is not in question. These are
# exempt from the DLIT requirement below.
def _controlled():
    g = globals()
    return {g[n] for n in g if n.startswith("ANNO_") and isinstance(g[n], str)} | {UNCHAR_ANNO}

# Strings that describe a protein's function rather than merely naming it. A
# claim like this needs a citation a human has actually read.
_FUNCTIONAL = re.compile(
    r"suppressor|antagonist|polymerase|protease|helicase|channel|viroporin|"
    r"movement|matrix|nucleocapsid|phospho|glycoprotein|peptide|factor|"
    r"paralog|virion|hydrophobic", re.I)


def enforce_designations(taxa):
    CONTROLLED = _controlled()
    """Remove unsupported designations and unsupported functional claims.

    Two rules, both about not asserting more than the evidence carries:

    1. A designation (U1, Gx, Px) goes in the annotation only for features in
       BACKED_DESIGNATION, i.e. where the source data already uses it.

    2. A functional claim needs a DLIT a human has read. If every citation on
       a feature is in CLAUDE_PMIDS -- proposed by a model, resolving in PubMed
       and topically plausible, but unread -- the claim is not established, and
       the annotation falls back to the uncharacterized form. The citation
       stays on the feature so it can be checked and the name restored.
    """
    edits = []
    for module, block in taxa.items():
        for key, ent in block.get("features", {}).items():
            old = ent.get("anno", "")
            new = old
            if key not in BACKED_DESIGNATION:
                new = _DESIG.sub("", new)
            pm = ent.get("PMID") or []
            # Established vocabulary needs no per-feature citation. Calling a
            # 2,100 aa rhabdovirus ORF "RNA-dependent RNA polymerase" does not
            # rest on any particular paper, and neither does a designation the
            # source data already uses. Only a NOVEL functional claim has to be
            # carried by a citation someone has read.
            exempt = new in CONTROLLED or key in BACKED_DESIGNATION
            #  Now that PMID holds only curated citations, "every citation is
            #  model-proposed" is simply "there are proposed ones and no
            #  curated ones". Reading pm alone would silently stop firing.
            proposed = ent.get("PMID_claude_generated") or []
            if (proposed and not pm and not exempt
                    and _FUNCTIONAL.search(new) and new != UNCHAR_ANNO):
                new = UNCHAR_ANNO
            if new != old:
                ent["anno"] = new
                edits.append((module, key, old, new))
    return edits


def feat(anno, symbol, lo, hi, *, ftype="CDS", segment=SINGLE, quality=True,
         kmers=None, pmid=None, cov=0.65, up=1, down=1, special=None,
         partner=None, bit=None, cleave=None):
    """cleave -- which termini of a mature peptide are protease cleavage sites.

    "n", "c" or "nc". A cleaved terminus is not a start or a stop codon, so the
    annotator must not go looking for one there: upstream_ext / downstream_ext
    are forced to 0 on that side.

    Getting this wrong is not a subtle error. G_MAT was left at the default
    upstream_ext=1, so the annotator scanned upstream from the mature
    N-terminus, found the precursor's own start codon, and emitted the mature
    peptide on exactly the same coordinates as the precursor -- rabies G and
    G_mature both starting at 3317 instead of the mature form starting 19
    codons in.
    """
    if cleave:
        if "n" in cleave:
            up = 0
        if "c" in cleave:
            down = 0
    mn, mx = bounds(lo, hi)
    d = {"anno": anno, "gene_symbol": symbol, "feature_type": ftype,
         "segment": segment, "min_len": mn, "max_len": mx}
    if special:
        d["special"] = special
    else:
        d["bit_cutoff"] = bit if bit is not None else bit_for(mn)
        if bit is not None:
            #  a hand-tuned cutoff -- the signal-peptide profiles in particular
            #  were set by measurement, high enough to reject one-off matches
            #  and low enough to still detect the peptide. derive_bounds() must
            #  not recompute those from min_len. Stripped before serialization.
            d["_bit_explicit"] = True
        d["coverage_cutoff"] = cov
        d["upstream_ext"] = up
        d["downstream_ext"] = down
        d["kmers"] = kmers if kmers is not None else (1 if quality else 0)
    if quality:
        d["copy_num"] = 1
    if pmid:
        ids = [str(p) for p in pmid]
        #  The two fields are DISJOINT, per JJD. A model-proposed citation
        #  appears only in PMID_claude_generated, never in PMID.
        #
        #  Listing it in both was the wrong default. PMID is where a reader
        #  looks for the evidence behind a feature, and a proposal that
        #  resolves in PubMed and reads plausibly is indistinguishable there
        #  from one a curator checked. Putting it in both fields means the
        #  flag has to be noticed to do its job, and a flag that has to be
        #  noticed is not a safeguard.
        #
        #  So a feature whose every citation is model-proposed ends up with NO
        #  PMID field at all, which is the honest statement: nothing here has
        #  been read. Five features are in that position. When a curator
        #  confirms one, remove it from CLAUDE_PMIDS and it moves to PMID by
        #  itself on the next regeneration.
        proposed = [p for p in ids if p in CLAUDE_PMIDS]
        curated = [p for p in ids if p not in CLAUDE_PMIDS]
        if curated:
            d["PMID"] = curated
        if proposed:
            d["PMID_claude_generated"] = proposed
    if partner:
        d["non_pssm_partner"] = partner
    return d


def seg(lo, hi):
    mn, mx = bounds(lo, hi)
    return {"min_len": mn, "max_len": mx, "replicon_geometry": "linear"}


def unchar(name=None, symbol="", lo=50, hi=250, **kw):
    """S1 Table's uncharacterized-protein forms."""
    anno = ("Uncharacterized lineage-specific %s protein" % name if name
            else "Uncharacterized lineage-specific protein")
    return feat(anno, symbol, lo, hi, quality=False, kmers=0, **kw)


# --------------------------------------------------------------------------
# Controlled vocabulary, family-wide. One string per function.
# --------------------------------------------------------------------------
ANNO_N = "Nucleocapsid protein"
ANNO_P = "Phosphoprotein"
ANNO_M = "Matrix protein"
ANNO_L = "RNA-dependent RNA polymerase"
ANNO_G = "Envelope glycoprotein precursor"      # precursor: SP is cleaved
#  G_MAT: ONLY the N-terminus is cleaved.
#
#  These carried down=0 as well, which says "do not scan to a stop, the
#  C-terminus is a cut site". It is not. Signal peptidase removes the leader
#  and nothing else; the mature chain runs to the precursor's own stop codon.
#  RABV is the arithmetic: a 524 aa precursor is 19 aa of signal plus 505 aa of
#  mature chain, so the two C-termini are the same residue.
#
#  Measured cost of getting it wrong, over 1108 annotated genomes: 45% of
#  mature-G calls ended somewhere other than the precursor's C-terminus, median
#  28 codons short and up to 321, because the end was wherever tblastn's
#  alignment happened to stop. With down=1 the scan runs to the stop codon and
#  the end is determined by the sequence instead.
#
#  Worth stating generally, because it costs more elsewhere: set cleave only on
#  termini that are genuinely protease cut sites. A polyprotein with several
#  internal products has two such termini per product, and every one left to
#  drift compounds along the chain.
ANNO_G_MAT = "Mature envelope glycoprotein"
ANNO_G_SP = "Signal peptide of G"
ANNO_MOV = "Movement protein"                    # matches Fimoviridae

# --------------------------------------------------------------------------
# DLIT PROVENANCE
# --------------------------------------------------------------------------
# PMIDs proposed by a language model rather than chosen by a human curator.
#
# Every one below was checked to resolve in PubMed and to be topically correct
# for the feature it is attached to. That is NOT the same as a curator having
# read the paper and confirmed it establishes the function or coordinates the
# feature claims, which is what a DLIT is supposed to mean. Until someone does
# that, these are proposals.
#
# They are emitted into the JSON as `PMID_claude_generated` alongside the
# ordinary `PMID` list, so a reader can tell curated citations from proposed
# ones at a glance and so they can be re-checked or replaced. Delete an entry
# from this dict once a human has verified it; the flag disappears on the next
# regeneration.
CLAUDE_PMIDS = {
    "1413521":  "BEFV genome contains two related glycoprotein genes (EPHEM_GNS)",
    "4038520":  "IHNV mRNA species (Novirhabdovirus NV)",
    "8337841":  "Adelaide River virus polycistronic glycoprotein genes (EPHEM_GNS)",
    "8607260":  "VSV replicates without C proteins (VESI_C)",
    "8683214":  "identification of the fish rhabdovirus NV protein (NV)",
    "9191923":  "GNS-L intergenic region organisation (EPHEM alpha/beta/gamma)",
    "22305623": "Kotonkan and Obodhiang ephemeroviruses (EPHEM_DELTA)",
    "24257609": "BEFV alpha1 has viroporin-like properties (EPHEM_ALPHA1)",
    "26700068": "cytorhabdovirus P3 encodes a 30K movement protein (MOV)",
    "28276468": "novirhabdovirus NV recruits PPM1Bb (NV)",
    "38140643": "tri-segmented rhabdoviruses described (Trirhavirus)",
    "38819305": "plant rhabdovirus viroporin P9 (ACYTO_VIROPORIN); replaces the "
                "DOI 10.1093/plcell/koae162 that was previously in this field",
}


PMID_FAMILY = ["35723908"]        # ICTV Rhabdoviridae profile
PMID_ACC = ["25679389"]           # Walker 2015 accessory-ORF survey
#  ALPHARHABDOVIRINAE ONLY. 6897030 is the direct sequencing of rabies virion
#  G that establishes the 19 aa signal + 505 aa mature chain, and the
#  Alpharhabdovirinae clusters carrying the KFP[ILM]YTIP mature N-terminus are
#  the same measurement. It was previously attached to G_SP and G_MAT in all
#  four modules, which asserted evidence that does not exist: the
#  Betarhabdovirinae, Novirhabdovirus and Dichorhavirus cleavage sites come
#  from SignalP 6.0, not from this paper. A rabies citation on a plant
#  rhabdovirus signal peptide is the U-number problem in another costume -- a
#  real reference projected past what it covers. Those three now carry no PMID;
#  their evidence is the per-cluster method, probability and conservation
#  recorded in each feature's PROVENANCE.
PMID_G_SP = ["6897030"]           # RABV G precursor / mature N-terminus

TAXA = {}

# ==========================================================================
# 1. Alpharhabdovirinae -- the general single-segment animal-rhabdovirus
#    module. Despite the name it is routed by ARCHITECTURE, not taxonomy:
#    everything with the canonical N-P-M-G-L plan plus positionally-numbered
#    accessories lands here, including all of Deltarhabdovirinae and the
#    Gammarhabdovirinae genus Margarhavirus.
#
#    Alpharhabdovirinae (>=5 genomes): Lyssavirus, Vesiculovirus,
#      Ephemerovirus, Sigmavirus, Hapavirus, Tibrovirus, Ledantevirus,
#      Sripuvirus, Almendravirus, Curiovirus, Tupavirus, Sunrhavirus,
#      Barhavirus, Sawgrhavirus, Arurhavirus, Perhabdovirus, Sprivivirus,
#      Ohlsrhavirus, Lostrhavirus, Caligrhavirus, Alphapaprhavirus,
#      Alphanemrhavirus, Alphathriprhavirus, Amplylivirus, Merhavirus.
#    Deltarhabdovirinae (>=5): Primrhavirus, Stangrhavirus,
#      Betaricinrhavirus, Betapaprhavirus, Betanemrhavirus, Alphahymrhavirus,
#      Gammahymrhavirus, Gammaricinrhavirus. All eight carry the canonical
#      five; every extra is a numbered unknown.
#    Below the cutoff, core five folded in, no dedicated accessory PSSMs:
#      Scophrhavirus (4, Px), Cetarhavirus (4), Alphacrustrhavirus (4, Px),
#      Alphadrosrhavirus (4, U1a/U1b), Mousrhavirus (3), Betahymrhavirus
#      (3, U1/U1x), Zarhavirus (2), Replylivirus (2), Uniorhavirus (2),
#      Betathriprhavirus (2, U1), Margarhavirus (2, U1/U2).
#
#    NOT folded in, and why:
#      * The three platrhavirus genera. Betaplatrhavirus "commonly lacks the
#        G protein", which would force copy_num off G for all 33 genera above
#        to accommodate three unassigned metagenomic genera. ICTV also
#        contradicts itself there (body text says four genes, its own figure
#        legend says five), and at a 9.9 kb lower bound a missing internal G
#        is equally consistent with a truncated assembly. Deferred.
#      * The two gymnorhavirus genera -- genuinely different proteins, see
#        module 7.
#
#    Known routing caveat: Uniorhavirus has TWO glycoprotein paralogs
#    (ICTV G1/G2; the primary paper calls them G and GNS), which overlap in
#    frame by 40 nt. With n=2 there is no basis for a second G PSSM, so the
#    G PSSM may call both and trip copy_num on G. Flagged, not solved.
# ==========================================================================
_alpha = {
    "N": feat(ANNO_N, "N", 415, 480, pmid=PMID_FAMILY),
    "P": feat(ANNO_P, "P", 250, 330, pmid=PMID_FAMILY),
    "M": feat(ANNO_M, "M", 190, 245, pmid=PMID_FAMILY),
    "G": feat(ANNO_G, "G", 490, 640, pmid=PMID_FAMILY),
    "L": feat(ANNO_L, "L", 2020, 2150, pmid=PMID_FAMILY),
    # G is a class I membrane glycoprotein with a cleaved signal peptide.
    # RABV: 524 aa precursor, 19 aa SP, 505 aa mature chain.
    "G_SP": feat(ANNO_G_SP, "SP", 16, 30, ftype="mat_peptide",
                 quality=False, kmers=0, cleave="c", pmid=PMID_G_SP),
    "G_MAT": feat(ANNO_G_MAT, "G_mature", 470, 620, ftype="mat_peptide",
                  quality=False, kmers=0, cleave="n", pmid=PMID_G_SP),

    # --- Vesiculovirus: overlapping alternative frame inside P, leaky
    #     scanning. ICTV contradicts itself on which label is longer; primary
    #     data (Spiropoulou & Nichol 1993) gives the minor form +10 residues,
    #     so C = 55 and C-prime = 65.
    "VESI_C": feat("C protein", "Vesi_C", 55, 55, quality=False, kmers=0,
                   pmid=["8388490", "8610460", "8607260"]),
    "VESI_CPRIME": feat("C-prime protein", "Vesi_C-Prime", 65, 67,
                        quality=False, kmers=0, pmid=["8388490", "8610460"]),

    # --- Ephemerovirus. alpha2 is expressed by termination-reinitiation via a
    #     TURBS: a purely translational mechanism, so it is an ordinary
    #     templated ORF and needs no special module.
    "EPHEM_GNS": feat("Nonstructural envelope glycoprotein GNS", "GNS",
                      534, 609, quality=False, kmers=0,
                      pmid=["1413521", "8337841"]),
    "EPHEM_ALPHA1": feat("Class Ia viroporin alpha1 protein",
                         "Ephemero_alpha1", 88, 108, quality=False, kmers=0,
                         pmid=["24257609", "9191923"]),
    "EPHEM_ALPHA2": unchar("alpha2", "Ephemero_alpha2", 92, 116,
                           pmid=["25679389", "9191923"]),
    "EPHEM_BETA": unchar("beta", "Ephemero_beta", 146, 157, pmid=["9191923"]),
    "EPHEM_GAMMA": unchar("gamma", "Ephemero_gamma", 100, 115,
                          pmid=["9191923"]),
    # BEFV's third alpha ORF. The binning rule has always been here; the
    # re-export is what finally put sequences behind it (13 of them). BV-BRC
    # calls them "alpha 3 protein", so the designation is source-backed.
    "EPHEM_ALPHA3": unchar("alpha3", "Ephemero_alpha3", 45, 70,
                           pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"EPHEM_DELTA": unchar("delta", "Ephemero_delta", 109, 109,
    #                      pmid=["22305623"]),

    # --- Sigmavirus X: between P and M, no cross-species homology.
    #     BV-BRC currently carries this as product "PP3", gene "X".
    "SIGMA_X": unchar("X", "Sigma_X", 224, 321, pmid=PMID_ACC),

    # --- Hapavirus PMIPs: the one sequence-established accessory family in
    #     the subfamily. NGAV alone has a GNS-like glycoprotein.
    "HAPA_U1": unchar("U1", "Hapa_U1", 100, 250, pmid=PMID_ACC),
    "HAPA_U2": unchar("U2", "Hapa_U2", 100, 250, pmid=PMID_ACC),
    "HAPA_U3": unchar("U3", "Hapa_U3", 100, 250, pmid=PMID_ACC),

    # --- Tibrovirus. NB the ICTV genome paragraph miscalls the G-L gene "U1";
    #     the figure legend and distinguishing features both say U3.
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"TIBRO_U1": unchar("U1", "Tibro_U1", 170, 216, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"TIBRO_U2": unchar("U2", "Tibro_U2", 100, 200, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"TIBRO_U3": feat("Class Ia viroporin U3 protein", "Tibro_U3", 60, 120,
    #                 quality=False, kmers=0, pmid=PMID_ACC),

    # --- Sripuvirus. Mx (stop-start/TURBS) and Gx (leaky scanning) in all
    #     members; U1 arose by duplication of the P gene in CHOV and SMV.
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"SRIPU_U1": feat("Phosphoprotein paralog U1 protein", "Sripu_U1",
    #                 150, 320, quality=False, kmers=0, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"SRIPU_MX": unchar("Mx", "Sripu_Mx", 50, 150, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"SRIPU_GX": unchar("Gx", "Sripu_Gx", 50, 150, pmid=PMID_ACC),

    # --- Viroporins. Positionally scattered and too divergent to establish
    #     orthology (Walker 2015), so one entry per lineage, never shared.
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"ALMEN_U1": feat("Class Ia viroporin U1 protein", "Almen_U1", 51, 80,
    #                 quality=False, kmers=0, pmid=PMID_ACC),
    "CURIO_VIROPORIN": feat("Class Ia viroporin protein", "Curio_viroporin",
                            55, 130, quality=False, kmers=0, pmid=PMID_ACC),
    "HAPA_VIROPORIN": feat("Class Ia viroporin protein", "Hapa_viroporin",
                           55, 130, quality=False, kmers=0, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"BARHA_GX": feat("Class Ia viroporin Gx protein", "Barha_Gx", 55, 110,
    #                 quality=False, kmers=0, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"ARURHA_GX": feat("Class Ia viroporin Gx protein", "Arurha_Gx", 55, 110,
    #                  quality=False, kmers=0, pmid=PMID_ACC),

    # --- Curiovirus P-M intergenic protein. The viroporin slot carries a
    #     different U-number per virus (U3 in CURV/IRIRV, U4 in RBUV, U2 in
    #     ITAV), so it is keyed on function above, not on the number.
    "CURIO_PMIP": unchar("U1", "Curio_U1", 90, 250, pmid=PMID_ACC),

    # --- Tupavirus / Sunrhavirus small hydrophobic protein. SH proteins share
    #     hydropathy but not orthology-level identity outside close relatives.
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"TUPA_SH": feat("Small hydrophobic protein", "Tupa_SH", 77, 93,
    #                quality=False, kmers=0, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"TUPA_PX": unchar("Px", "Tupa_Px", 40, 130, pmid=PMID_ACC),
    "SUNRHA_SH": feat("Small hydrophobic protein", "Sunrha_SH", 77, 93,
                      quality=False, kmers=0, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"SUNRHA_PX": unchar("Px", "Sunrha_Px", 40, 130, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"SUNRHA_U2": unchar("U2", "Sunrha_U2", 19, 200, pmid=PMID_ACC),

    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"LEDAN_U1": unchar("U1", "Ledan_U1", 60, 200, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"CALIG_U1": unchar("U1", "Calig_U1", 120, 160, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"BARHA_NX": unchar("Nx", "Barha_Nx", 60, 200, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"SAWGR_GX": unchar("Gx", "Sawgr_Gx", 50, 180, pmid=PMID_ACC),
    "SAWGR_GY": unchar("Gy", "Sawgr_Gy", 50, 180, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"ARURHA_U1": unchar("U1", "Arurha_U1", 55, 200, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"ARURHA_U2": unchar("U2", "Arurha_U2", 55, 200, pmid=PMID_ACC),

    # --- Merhavirus / CTRV. The ONLY rhabdovirus where naive translation of
    #     the genome gives a wrong protein: a 76 nt GU-AG spliceosomal intron
    #     at antigenome 8648-8723 interrupts L. Excision shifts +2 -> +3 and
    #     yields one 2123 aa ORF. Without this, L is emitted truncated at the
    #     in-frame stop at 8681.
    #     Symbol is prefixed to keep it distinct from the ordinary L PSSM in
    #     this same collapsed module. Note the ordinary L PSSM will also hit
    #     CTRV and stop at the in-frame stop at 8681, so the splice module
    #     needs to win: check that get_splice_variant_features.pl replaces
    #     rather than duplicates the L call on CTRV contigs.
    "MERHA_L_SPLICED": feat(ANNO_L, "Merha_L", 2123, 2123, quality=False,
                            special="splice", pmid=["21507977"]),

    # ---------------------------------------------------------------------
    # Deltarhabdovirinae. Every genus in that subfamily carries the canonical
    # five; the extras are all positionally-numbered unknowns, so they fold
    # in here rather than needing a module of their own. Dedicated accessory
    # features only where the genus clears the 5-genome cutoff -- below that
    # a PSSM cannot meet the >=3-sequences-per-cluster rule anyway, so those
    # genera contribute their core five and nothing else.
    # ---------------------------------------------------------------------
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"BPAP_U1": unchar("U1", "Bpap_U1", 60, 200, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"BNEM_U1": unchar("U1", "Bnem_U1", 60, 250, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"GRICIN_U1": unchar("U1", "Gricin_U1", 280, 330, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"STANG_U1": unchar("U1", "Stang_U1", 50, 70, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"BRICIN_NX": unchar("Nx", "Bricin_Nx", 50, 200, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"BRICIN_PX": unchar("Px", "Bricin_Px", 50, 200, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"ATHRIP_GX": feat("Class Ia viroporin Gx protein", "Athrip_Gx",
    #                  121, 124, quality=False, kmers=0, pmid=PMID_ACC),

    # --- Catch-all for the metagenomic and thinly-described genera that pass
    #     the genome cutoff but whose accessory ORFs ICTV leaves unnamed
    #     (Gammahymrhavirus, Alphapaprhavirus, Alphanemrhavirus,
    #     Alphaplatrhavirus, Ohlsrhavirus, Lostrhavirus, Primrhavirus alt
    #     ORFs). Empty symbol per S1.
}

TAXA["Alpharhabdovirinae"] = {
    "close_genomes": {
            #  Added 11 Sep 2026 after the whole-taxon coverage run. These six
            #  genera contributed sequences to this module's profiles and then
            #  could not route a single genome to them: the BLASTn step found
            #  no reference above 150 bits for any of Almendravirus, Alphaplatrhavirus, Primrhavirus, Betaplatrhavirus.
            #  That is a wiring gap, not a modelling gap -- the profiles were
            #  trained on their N, P, M, G and L.
            #
            #  Alpharhabdovirinae.10 (Alphaplatrhavirus) and .12 (Betaplatrhavirus) were
            #  added here on 11 Sep and REMOVED the same day. They routed 8 genomes
            #  and every one came back with L alone -- no N, P, M or G -- because the
            #  platrhaviruses contributed sequences to this module's collections but
            #  no profile of theirs survives clustering. Three further Alphaplatrhavirus
            #  references were tried and recovered 6 more genomes, also L-only.
            #
            #  A reference that routes a genome the profiles cannot describe makes the
            #  routed percentage look better and the annotation no better. Worse, the
            #  genome goes from cleanly unrouted to annotated-and-flagged-poor for four
            #  missing essential features. The platrhaviruses need profiles, not
            #  references; until they have them they belong outside the module.
            #
            #  Chosen on taxonomic coverage, NOT on yield. The greedy set cover
            #  over all 315 unrouted genomes is almost flat (5 references buy
            #  10%, 50 buy 39%, full coverage needs 242) because 200 of them
            #  are singletons. Adding references to chase that tail is not
            #  worth it. Making a genus routable at all is, and it keeps paying
            #  as new isolates are deposited.
            "Alpharhabdovirinae.9.dna": {
                "genome_ids": "LC270812", "genome_name": "Menghai rhabdovirus kunoichi"},
            "Alpharhabdovirinae.11.dna": {
                "genome_ids": "BK059423", "genome_name": "San Gabriel mononegavirus San Gabriel Valley"},
        # One rep-contig file per genus; all route to this single module.
        # See Rhabdoviridae-Rep-Contigs.tsv for accessions still to resolve.
        "Alpharhabdovirinae.1.dna": {
            "genome_ids": "NC_001542.1", "genome_name": "Rabies lyssavirus Pasteur virus"},
        "Alpharhabdovirinae.2.dna": {
            "genome_ids": "NC_001560.1",
            "genome_name": "Vesicular stomatitis Indiana virus"},
        "Alpharhabdovirinae.3.dna": {
            "genome_ids": "NC_002526.1",
            "genome_name": "Bovine ephemeral fever virus BB7721"},
        "Alpharhabdovirinae.4.dna": {
            "genome_ids": "NC_038279.1",
            "genome_name": "Drosophila ananassae sigmavirus Kilifi Kenya"},
        "Alpharhabdovirinae.5.dna": {
            "genome_ids": "NC_034447.1", "genome_name": "Hart Park virus"},
        "Alpharhabdovirinae.6.dna": {
            "genome_ids": "NC_020804.1", "genome_name": "Tibrogargan virus"},
        "Alpharhabdovirinae.7.dna": {
            "genome_ids": "NC_002803.1", "genome_name": "Spring viraemia of carp virus"},
        "Alpharhabdovirinae.8.dna": {
            "genome_ids": "NC_025384.1",
            "genome_name": "Culex tritaeniorhynchus rhabdovirus"},
    },
    "features": _alpha,
    # Wide because the module now spans 10.1-16.1 kb of genome sizes. This is
    # the cost of collapsing; tighten per-genus only if the quality module
    # turns out to pass too much.
    "segments": {SINGLE: seg(10100, 16100)},
}

# ==========================================================================
# 2. Novirhabdovirus (Gammarhabdovirinae) -- kept separate. 3090 genomes,
#    fish pathogens, and the NV gene. Margarhavirus (n=2) is dropped, so this
#    module is the whole of the subfamily worth modelling.
# ==========================================================================
_novi = {
    # No UNCHAR for Novirhabdovirus: all 13 of its unnamed proteins turned out
    # to be near-identical to members of M, NV or P, so the grab-bag was empty
    # once made disjoint -- and while it existed it produced a duplicate call
    # on the phosphoprotein locus of both test genomes.
    "N": feat(ANNO_N, "N", 391, 404, pmid=PMID_FAMILY),
    "P": feat(ANNO_P, "P", 222, 230, pmid=PMID_FAMILY),
    "M": feat(ANNO_M, "M", 195, 201, pmid=PMID_FAMILY),
    "G": feat(ANNO_G, "G", 507, 508, pmid=PMID_FAMILY),
    "L": feat(ANNO_L, "L", 1984, 1986, pmid=PMID_FAMILY),
    "G_SP": feat(ANNO_G_SP, "SP", 16, 30, ftype="mat_peptide",
                 quality=False, kmers=0, cleave="c"),
    "G_MAT": feat(ANNO_G_MAT, "G_mature", 480, 495, ftype="mat_peptide",
                  quality=False, kmers=0, cleave="n"),
    # NV has its own transcription unit between G and L and is universal in
    # the genus, but NV sequences show NO detectable similarity across the
    # four species, so a single genus-wide PSSM will fail. Build per-species
    # PSSMs in the Novirhabdovirus directory. kmers=0 for the same reason.
    "NV": feat("Non-virion protein", "NV", 111, 130, kmers=0,
               pmid=["4038520", "8683214", "28276468", "9010293"]),
}
TAXA["Novirhabdovirus"] = {
    "close_genomes": {
        "Novirhabdovirus.1.dna": {
            "genome_ids": "NC_001652.1",
            "genome_name": "Infectious hematopoietic necrosis virus"},
        "Novirhabdovirus.2.dna": {
            "genome_ids": "NC_000855.1",
            "genome_name": "Viral hemorrhagic septicemia virus isolate Fil3"},
    },
    "features": _novi,
    "segments": {SINGLE: seg(11100, 11500)},
}

# ==========================================================================
# 3. Betarhabdovirinae -- the single-segment plant genera whose gene set is
#    recognisably N-P-Mov-M-G-L: Alpha/Beta/Gammacytorhabdovirus and
#    Alpha/Beta/Gamma/Deltanucleorhabdovirus.
#
#    The gymnorhaviruses are deliberately NOT here -- see module 7. They have
#    no recognisable P, M, G or Mov, so folding them in would align their
#    P2/P3/P4 against Mov/M/G on genomic position alone.
# ==========================================================================
PMID_MOV = ["26700068", "8091658", "8178430", "40999317"]

_beta = {
    "N": feat(ANNO_N, "N", 397, 521, pmid=PMID_FAMILY),
    "P": feat(ANNO_P, "P", 237, 330, pmid=PMID_FAMILY),
    "L": feat(ANNO_L, "L", 1877, 2130, pmid=PMID_FAMILY),
    # Movement protein. Called P3, sc4, 4b, Y, ORF3 and "protein 3" in the
    # literature for the same protein; all collapse to this one string, with
    # the symbol matching Fimoviridae.
    "MOV": feat(ANNO_MOV, "Mov", 290, 340, quality=False, kmers=1,
                pmid=PMID_MOV),
    # G and M are absent in whole genera (no gammacytorhabdovirus has a G;
    # several also lack M), so neither can carry copy_num.
    "G": feat(ANNO_G, "G", 500, 669, quality=False, kmers=1,
              pmid=PMID_FAMILY),
    "M": feat(ANNO_M, "M", 170, 360, quality=False, kmers=1,
              pmid=PMID_FAMILY),
    # Plant-rhabdovirus G is also a class I membrane glycoprotein, but no
    # signal-peptide coordinates are published for any member. These two are
    # inferred by homology to RABV G, which is the paper's stated policy for
    # unvalidated cleavage sites -- but they are inference, not measurement.
    "G_SP": feat(ANNO_G_SP, "SP", 16, 35, ftype="mat_peptide",
                 quality=False, kmers=0, cleave="c"),
    "G_MAT": feat(ANNO_G_MAT, "G_mature", 480, 650, ftype="mat_peptide",
                  quality=False, kmers=0, cleave="n"),

    # RYSV P6 is virion-associated and an RNA silencing suppressor. BYSMV P6
    # is a jasmonate-signalling effector -- a different protein with the same
    # number, hence separate prefixed entries.
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"ANUCLEO_P6": feat("RNA silencing suppressor P6 protein", "ANucleo_P6",
    #                   150, 250, quality=False, kmers=0, pmid=["23634838"]),
    "ANUCLEO_X": unchar("X", "ANucleo_X", 50, 250, pmid=["20362316"]),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"ANUCLEO_P5": unchar("P5", "ANucleo_P5", 50, 250, pmid=PMID_ACC),
    "BNUCLEO_U": unchar("U", "BNucleo_U", 50, 250, pmid=PMID_ACC),

    # Gammanucleorhabdovirus P3 and P4 have NO homology to other rhabdovirus
    # proteins, "including the P3 movement proteins of other plant
    # rhabdoviruses" (ICTV, after Tsai et al. 2005). They must NOT be folded
    # into the Movement protein entry above.
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"GNUCLEO_P3": unchar("P3", "GNucleo_P3", 250, 350, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"GNUCLEO_P4": unchar("P4", "GNucleo_P4", 50, 250, pmid=PMID_ACC),

    "ACYTO_VIROPORIN": feat("Class Ia viroporin P9 protein", "ACyto_P9",
                            55, 130, quality=False, kmers=0,
                            pmid=["38819305"]),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"ACYTO_PX": unchar("Px", "ACyto_Px", 50, 200, pmid=PMID_ACC),
    "BCYTO_P4": unchar("P4", "BCyto_P4", 50, 250, pmid=PMID_ACC),
    "BCYTO_P5": unchar("P5", "BCyto_P5", 50, 250, pmid=PMID_ACC),
    #  RETIRED 11 Sep 2026 per JJD: positional label with too few
    #  sequences to build a profile. Not enough data is not a reason
    #  to ship a call. Sequences now flow to the UNC homology groups.
    #"BCYTO_P6": unchar("P6", "BCyto_P6", 50, 250, pmid=PMID_ACC),

}
TAXA["Betarhabdovirinae"] = {
    "close_genomes": {
            #  Added 11 Sep 2026 after the whole-taxon coverage run. These six
            #  genera contributed sequences to this module's profiles and then
            #  could not route a single genome to them: the BLASTn step found
            #  no reference above 150 bits for any of Gammacytorhabdovirus, Deltanucleorhabdovirus.
            #  That is a wiring gap, not a modelling gap -- the profiles were
            #  trained on their N, P, M, G and L.
            #
            #  Chosen on taxonomic coverage, NOT on yield. The greedy set cover
            #  over all 315 unrouted genomes is almost flat (5 references buy
            #  10%, 50 buy 39%, full coverage needs 242) because 200 of them
            #  are singletons. Adding references to chase that tail is not
            #  worth it. Making a genus routable at all is, and it keeps paying
            #  as new isolates are deposited.
            "Betarhabdovirinae.7.dna": {
                "genome_ids": "BK064345", "genome_name": "Argyranthemum gammacytorhabdovirus 1"},
            "Betarhabdovirinae.8.dna": {
                "genome_ids": "BK070508", "genome_name": "Artemisia deltanucleorhabdovirus 1"},
        "Betarhabdovirinae.1.dna": {
            "genome_ids": "NC_007642.1", "genome_name": "Lettuce necrotic yellows virus"},
        "Betarhabdovirinae.2.dna": {
            "genome_ids": "NC_001615.3", "genome_name": "Sonchus yellow net virus"},
        "Betarhabdovirinae.3.dna": {
            "genome_ids": "NC_016136.1", "genome_name": "Potato yellow dwarf virus"},
        "Betarhabdovirinae.4.dna": {
            "genome_ids": "NC_005974.1", "genome_name": "Maize fine streak virus"},
        "Betarhabdovirinae.5.dna": {
            "genome_ids": "NC_002251.1", "genome_name": "Northern cereal mosaic virus"},
        # Added Sep 2026: Rice yellow stunt virus returned ZERO BLASTn hits
        # against the five rep contigs above and was rejected before any PSSM
        # ran. Plant rhabdoviruses are too divergent at the nucleotide level for
        # one exemplar per genus.
        "Betarhabdovirinae.6.dna": {
            "genome_ids": "NC_003746.1", "genome_name": "Rice yellow stunt virus"},
    },
    "features": _beta,
    "segments": {SINGLE: seg(9900, 16700)},
}

# ==========================================================================
# 4-6. The segmented plant genera. These CANNOT be collapsed into
#      Betarhabdovirinae because each needs its own `segments` schema.
#
#      NB the polymerase segment index is INVERTED between the bisegmented
#      genera: Dichorhavirus puts L alone on Segment 2, Varicosavirus and
#      Trirhavirus put L alone on Segment 1. Any code assuming "L is always
#      on segment n" mis-annotates one of them.
# ==========================================================================
S1, S2, S3 = "Segment 1", "Segment 2", "Segment 3"

TAXA["Dichorhavirus"] = {
    "close_genomes": {
        "Dichorhavirus.1.dna": {
            "genome_ids": "NC_009608.1,NC_009609.1", "genome_name": "Orchid fleck virus"}},
    "features": {
        # No UNCHAR for Dichorhavirus either: 27 of its 33 unnamed proteins
        # were near-identical to N, M, MOV, G or P, and the 6 that remain do
        # not cluster into a profile at any setting.
        "N": feat(ANNO_N, "N", 450, 450, segment=S1, pmid=PMID_FAMILY),
        "P": feat(ANNO_P, "P", 237, 237, segment=S1, pmid=PMID_FAMILY),
        "MOV": feat(ANNO_MOV, "Mov", 335, 335, segment=S1, pmid=PMID_MOV),
        "M": feat(ANNO_M, "M", 183, 183, segment=S1, pmid=PMID_FAMILY),
        "G": feat(ANNO_G, "G", 542, 542, segment=S1, pmid=PMID_FAMILY),
        # G_SP was DROPPED for Dichorhavirus on 3 Sep 2026 and REINSTATED on
        # 10 Sep 2026, per JJD. Both decisions were right on the evidence at
        # the time, and the thing that changed was the clustering, not the data.
        #
        # The 3 Sep withdrawal: SignalP placed the cleavage at 20|21 (A|L),
        # invariant across all six sequences, so the coordinate was never in
        # doubt. The profile was. One PSSM built from all six scored 40.0 bits
        # on its best member and 23.4-25.0 on the other five, against a
        # bit_cutoff of 30 -- one of its own six training sequences cleared its
        # own cutoff. A feature the annotator can never call should not ship.
        #
        # What changed: N-terminal reclustering split G into 1_1, 1_2, 2_1,
        # 2_2, 3_1, so those six signal peptides now sit in homogeneous groups
        # instead of one mixed one. A 20-residue profile can only describe one
        # homogeneous set -- the same per-cluster principle established for
        # Alpharhabdovirinae. Re-derived and re-tested by the identical method
        # (psiblast, -comp_based_stats 0, against the same cutoff of 30):
        #
        #     1_1   n=7   38.7-42.1   7/7 clear 30
        #     1_2   n=4   43.3-45.6   4/4 clear 30
        #     2_2   n=3   46.7-46.7   3/3 clear 30
        #
        # Cluster 2_1 is absent by design: SignalP gave it probability 0.46,
        # below the 0.70 floor in scripts/apply_signalp.py, so it was refused.
        # The self-recall test is now enforced automatically by that script's
        # --self-cutoff, so a profile that cannot reach its cutoff is withheld
        # without anyone having to remember this episode.
        "G_SP": feat(ANNO_G_SP, "SP", 16, 30, ftype="mat_peptide",
                     segment=S1, quality=False, kmers=0, cleave="c"),
        "G_MAT": feat(ANNO_G_MAT, "G_mature", 505, 530, ftype="mat_peptide",
                      segment=S1, quality=False, kmers=0, cleave="n"),
        "L": feat(ANNO_L, "L", 1877, 1877, segment=S2, pmid=PMID_FAMILY),
    },
    "segments": {S1: seg(6400, 6700), S2: seg(5900, 6100)},
}

# Varicosavirus -- DROPPED as a module, 3 Sep 2026, per JJD.
#
# ICTV shows the genomes carry a set of accessory genes that cannot be
# supported with PSSMs from what BV-BRC holds. The dump has 262 genomes but
# only 73 in-range proteins across six collections, 53 of them the coat
# protein; 196 of the 260 genomes with any feature carry exactly one. L had
# three sequences, all Lettuce big-vein associated virus at 2040 aa, which
# is below the -m 5 clustering minimum, and the ORF2/ORF4/ORF5/ORF6 set was
# never buildable. Annotating the genus positionally off ICTV's hedged
# 2/3/4/5 -> P/Mov/M/G correspondence was not defensible without the
# sequences to back it, so the module is withdrawn rather than shipped thin.
#
# Its sequences are not discarded: build_collections.py records them in
# collections/NOT_MODELLED.tsv. Restoring the module means restoring this
# block and the Varicosavirus entry in build_collections.py's GENUS_MODULE.


# Trirhavirus: P6/P7/P8 on Segment 3 are the P/M/G equivalents. Segment 2
# "P3" is NOT a movement protein -- the genus has no identifiable MP, and
# only N and L are recognisably rhabdoviral.
# ICTV publishes NO length data of any kind for this genus, so every bound
# here is a PLACEHOLDER and must be replaced from real data before the
# quality module is trusted on it.
# --------------------------------------------------------------------------
# Trirhavirus -- DROPPED 10 Sep 2026, per JJD: "too small and too diverse".
#
# The re-export finally gave the genus sequences (37 proteins, 12 genomes,
# 4 species) and they showed why it cannot work: the species' proteins are
# 30-36% identical to each other. G, M and P produce no cluster of two at any
# setting tried, down to the -mi 0.6 floor with -mc 0.5. Five of nine built and
# four did not, so the module could never call a complete genome, which is the
# bar for shipping a taxon at all.
#
# Its proteins are not lost -- with the genus no longer routed they are
# recorded in collections/NOT_MODELLED.tsv with their taxonomic origin.
# --------------------------------------------------------------------------


# ==========================================================================
# 7. Alphagymnorhavirus -- the one straggler with genuinely different
#    proteins, so it gets its own module rather than being folded in.
#
#    Architecture is 3'-N-P2-P3-P4-L-5', with P5 between P4 and L in one
#    member. ICTV identifies ONLY N and L; there is no recognisable P, M, G
#    or movement protein anywhere in either gymnorhavirus genus, despite
#    these being plant viruses. The intervening proteins are named
#    positionally (P2, P3, P4, P5) and their functions are unknown.
#
#    ICTV twice notes that the intervening proteins have "distant similarity
#    to those of betagymnorhaviruses and varicosaviruses", and both genera
#    cluster phylogenetically with Varicosavirus -- so when these alignments
#    are curated, the varicosavirus ORF2/4/5 collections are the right thing
#    to compare against. They cannot simply join Varicosavirus, though,
#    because gymnorhavirus genomes are UNSEGMENTED and varicosaviruses are
#    bisegmented.
#
#    Alphagymnorhavirus has 14 BV-BRC genomes and clears the cutoff.
#    Betagymnorhavirus (n=2) has the same architecture and folds in here,
#    contributing sequence without getting features of its own.
# ==========================================================================
#    DROPPED 11 Sep 2026 per JJD, with Betagymnorhavirus which folded into it.
#    The module declared six features and built ZERO profiles. L had 7
#    sequences and N had 8, neither clustering at the -mi 0.6 floor, and
#    P2-P5 never had a collection at all -- positional labels with nothing
#    behind them.
#
#    Running a genome is what decided it. NC_139364.1, the module's own
#    reference, matched its rep contig and came back with a single
#    nucleocapsid, and even that came from a stale profile the installer had
#    failed to prune. A module with a rep contig and no profiles captures
#    genomes and dead-ends them; without it, routing offers them to
#    Alpharhabdovirinae and its 225 profiles instead.
#
#    Kept commented rather than deleted because the ICTV reasoning below is
#    sound and the blocker is sequence depth, not taxonomy. To revive:
#    uncomment this block and the two GENUS_MODULE lines in
#    build_collections.py. The varicosavirus comparison noted below is still
#    the right first step whenever these genera get more sequences.
# ==========================================================================
#TAXA["Alphagymnorhavirus"] = {
#    "close_genomes": {
#        "Alphagymnorhavirus.1.dna": {
#            "genome_ids": "NC_139364.1", "genome_name": "Cupressus virus 1"}},
#    "features": {
#        "N": feat(ANNO_N, "N", 400, 500, pmid=PMID_FAMILY),
#        "L": feat(ANNO_L, "L", 1900, 2130, pmid=PMID_FAMILY),
#        # Positionally numbered, functions unknown. Bounds are placeholders:
#        # no per-protein lengths are published for either genus.
#        "P2": unchar("P2", "Gymno_P2", 50, 400),
#        "P3": unchar("P3", "Gymno_P3", 50, 400),
#        "P4": unchar("P4", "Gymno_P4", 50, 400),
#        "P5": unchar("P5", "Gymno_P5", 50, 400),
#    },
#    "segments": {SINGLE: seg(10200, 12200)},
#}


def sort_deep(obj):
    if isinstance(obj, dict):
        return {k: sort_deep(obj[k]) for k in sorted(obj)}
    return obj


def declare_unchar_groups(taxa, workdir=os.path.dirname(os.path.abspath(__file__))):
    """Declare one feature per UNC<n> homology group that has a collection.

    The uncharacterized proteins are no longer one pool per module. Each is a
    homology group with its own collection, alignment directory and profiles,
    so each needs its own key here -- and they all carry the SAME annotation,
    the S1 Table form for an uncharacterized lineage-specific protein, because
    that is genuinely all that is known about any of them.

    Separate keys with one shared annotation is the whole point. The day UNC4
    turns out to be a viroporin, it is renamed in the JSON by itself; the other
    groups keep their annotation and their profiles. While they shared a key
    that was impossible, and the single profile they shared recalled 31% of its
    own collection because it was a model of nothing.

    Declared from what is on disk rather than written out by hand: the group
    count follows the data, and a rebuild that finds one more group gets one
    more feature without anyone editing this file.
    """
    added = []
    for module, block in taxa.items():
        feats = block.get("features")
        if feats is None:
            continue
        for fa in sorted(glob.glob(os.path.join(workdir, "collections", module,
                                                "UNC[0-9]*.fasta")),
                         key=lambda q: int(re.sub(r"\D", "", os.path.basename(q)) or 0)):
            key = os.path.basename(fa)[:-6]
            if key in feats:
                continue
            n = sum(1 for line in open(fa, errors="replace") if line.startswith(">"))
            if n < 3:
                continue
            lens = []
            cur = 0
            for line in open(fa, errors="replace"):
                if line.startswith(">"):
                    if cur:
                        lens.append(cur)
                    cur = 0
                else:
                    cur += len(line.strip())
            if cur:
                lens.append(cur)
            lens.sort()
            feats[key] = unchar(None, "", lens[0], lens[-1])
            added.append((module, key, n, lens[0], lens[-1]))
    return added


def derive_bounds(taxa, workdir=os.path.dirname(os.path.abspath(__file__))):
    """Re-derive min_len/max_len from the collections, at p1/p99 rather than min/max.

    The hand-set (lo, hi) in each feat() call came from the pre-re-export
    collections. Those collections then grew three to five fold, and the bounds
    did not follow: ten essential features ended up excluding their own
    reference sequences -- Alpharhabdovirinae P alone had 141 members shorter
    than its own min_len, so viral_genome_quality.pl would stamp "Feature is
    too short" on proteins the profile was built from.

    Percentiles, not min/max, because the collections still contain truncated
    proteins -- that is the defect the length flag exists to catch. One 191-aa
    fragment must not drag P's floor down far enough to stop flagging the next
    one. p1/p99 with the established +/-10% pad tolerates real variation and
    still flags genuine truncation: across the essential features it takes
    members-outside-their-own-bounds from 384 to 15.

    Below 20 sequences a percentile cannot distinguish an outlier from the
    distribution, so fall back to observed min/max; below 5, keep the hand-set
    value.
    """
    def lengths(path):
        out, cur = [], 0
        for line in open(path, errors="replace"):
            if line.startswith(">"):
                if cur:
                    out.append(cur)
                cur = 0
            else:
                cur += len(line.strip())
        if cur:
            out.append(cur)
        return sorted(out)

    def pct(v, q):
        i = (len(v) - 1) * q
        lo, hi = math.floor(i), math.ceil(i)
        return v[lo] if lo == hi else v[lo] + (v[hi] - v[lo]) * (i - lo)

    changed = []
    for module, block in taxa.items():
        for key, ent in block.get("features", {}).items():
            if ent.get("special") or "min_len" not in ent:
                continue
            fa = os.path.join(workdir, "collections", module, key + ".fasta")
            if not os.path.exists(fa):
                continue
            v = lengths(fa)
            if len(v) < 5:
                continue
            lo, hi = (pct(v, 0.01), pct(v, 0.99)) if len(v) >= 20 else (v[0], v[-1])
            mn, mx = int(math.floor(lo * 0.9)), int(math.ceil(hi * 1.1))
            if (mn, mx) != (ent["min_len"], ent["max_len"]):
                changed.append((module, key, len(v),
                                ent["min_len"], ent["max_len"], mn, mx))
                ent["min_len"], ent["max_len"] = mn, mx
            #  bit_cutoff deliberately does NOT follow min_len here. It is a
            #  specificity threshold, and min_len is now a 1st-percentile floor:
            #  letting one truncated protein at p1 drag the floor across a
            #  bit_for() boundary would loosen the cutoff for the whole feature.
            #  Alpharhabdovirinae P demonstrated it -- min_len 225 -> 174 crosses
            #  the 180 boundary and would have taken the cutoff 80 -> 60, making
            #  a core protein easier to call in the wrong genus in exchange for
            #  nothing. The cutoffs stay as feat() set them, from the hand-set
            #  bounds, and are tuned against measured cross-genus hits instead.
    return changed


def canonical(obj):
    """Serialise exactly as JSON::XS->new->pretty->canonical->encode would.

    One format, so a diff shows only real changes. 3-space indent and " : " are
    JSON::XS's `pretty`; the indent width is fixed at 3 in that module and is
    not configurable, so it is a fact rather than a preference. `canonical`
    sorts keys, which is the part that actually matters: Perl randomises hash
    order, so any Perl tool that reads and rewrites this file reorders every
    line differently on each run unless the keys are sorted.

    Integral floats are folded to ints because Perl has no separate float type
    for them and would write 1 where Python writes 1.0. That is the only place
    the two languages disagree on this schema.

    Same output as the skill's scripts/json_canon.py; kept inline so this
    generator has no dependency outside the working directory.
    """
    def norm(o):
        if isinstance(o, bool):
            return o
        if isinstance(o, float) and o.is_integer():
            return int(o)
        if isinstance(o, dict):
            #  keys beginning with "_" are generator bookkeeping, not schema
            return {k: norm(v) for k, v in o.items() if not k.startswith("_")}
        if isinstance(o, list):
            return [norm(v) for v in o]
        return o
    return json.dumps(norm(obj), indent=3, separators=(",", " : "),
                      sort_keys=True, ensure_ascii=False) + "\n"


if __name__ == "__main__":
    _unc = declare_unchar_groups(TAXA)
    for _m, _k, _n, _lo, _hi in _unc:
        sys.stderr.write("uncharacterized group declared: %-22s %-8s %4d seqs  %d-%d aa\n"
                         % (_m, _k, _n, _lo, _hi))
    sys.stderr.write("%d uncharacterized homology group(s) declared as features\n"
                     % len(_unc))
    _bounds = derive_bounds(TAXA)
    for _m, _k, _n, _omn, _omx, _nmn, _nmx in _bounds:
        sys.stderr.write("bounds re-derived: %-22s %-16s n=%-5d %d-%d -> %d-%d\n"
                         % (_m, _k, _n, _omn, _omx, _nmn, _nmx))
    sys.stderr.write("%d feature(s) had length bounds re-derived from the collections\n"
                     % len(_bounds))
    _edits = enforce_designations(TAXA)
    for _m, _k, _o, _n in _edits:
        sys.stderr.write("designation stripped: %-18s %-46s -> %s\n" % (_k, _o, _n))
    sys.stderr.write("%d annotation string(s) had an unsupported designation removed\n"
                     % len(_edits))
    sys.stdout.write(canonical(sort_deep(TAXA)))
    n_feat = sum(len(v["features"]) for v in TAXA.values())
    annos = {f["anno"] for v in TAXA.values() for f in v["features"].values()}
    sys.stderr.write(
        "modules: %d   features: %d   unique annotation strings: %d\n"
        % (len(TAXA), n_feat, len(annos)))
    sys.stderr.write("dropped, fewer than %d BV-BRC genomes (%d genera):\n  %s\n"
                     % (MIN_GENOMES, len(DROPPED), ", ".join(DROPPED)))
