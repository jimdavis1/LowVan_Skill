#!/usr/bin/env python3
"""
Build per-protein collections for the Togaviridae LowVan module from the
BV-BRC dump in this directory.

Inputs
    Toga.id_name_gs      genome_id \t genome_name \t family \t genus
    Toga.genome_meta.tsv genome_id \t taxon \t name \t length \t contigs \t status \t species \t subgenus \t acc
    Toga.feat_full.tsv   feature_id \t product \t md5 \t feature_type \t start \t end \t strand \t na_length
    Toga.id_md5          feature_id \t md5
    Toga.uniq.md5        md5                       (one per unique sequence)
    Toga.uniq.seq        md5 \t sequence
    Toga.uniq.id_ann     representative_feature_id \t modal annotation

uniq.seq / uniq.id_ann / uniq.md5 are line-aligned by construction here (see
build_dump.py), but this loader joins BY KEY anyway -- the alignment is a
property of one export, and a loader that depends on it breaks silently on the
next one.

Outputs
    synonyms.tsv                 every annotation string x species, with counts
    collections/Togaviridae/<FEATURE>.fasta
    collections/Togaviridae/<FEATURE>.outliers.fasta
    collections/report.tsv       per-collection counts and length ranges
    collections/UNASSIGNED_TRACKING.tsv   strings no rule matched
    collections/NOT_MODELLED.tsv          what was deliberately excluded, and why
    collections/MISSING.tsv      near-complete genomes lacking each core protein
"""

import collections
import os
import re
import shutil
import subprocess
import sys
import tempfile

MODULE = "Togaviridae"
AMBIG = set("BJOUXZ*")
MIN_AA = 30

# --------------------------------------------------------------------------
# How much assembly ambiguity a training sequence may carry.
#
# Dropping every sequence with an X costs 3,953 of 25,858 unique proteins here
# -- 15% of the dump -- and it does not cost them at random. X comes from
# read-assembled submissions, which is where the recent and the divergent
# sequences are; a hard filter quietly trades coverage of the genus for tidiness
# of the alignments. The distribution says the same: the median affected
# sequence is 1.15% X, i.e. a handful of positions in a 500 aa protein.
#
# So the line is drawn on proportion rather than presence. 1% admits 1,867 of
# the 3,953 (about five X in a 500 aa protein, which one PSSM column absorbs)
# and still rejects the 6%-and-up tail that would genuinely blur a profile.
# B/J/O/U/Z and internal stops are never admitted at any proportion -- unlike X
# they indicate a translation or a record that is wrong, not merely uncertain.
# Trailing and leading X are stripped first; they are padding, not data.
# --------------------------------------------------------------------------
MAX_AMBIG_FRAC = 0.01
HARD_AMBIG = set("BJOUZ*")

# ==========================================================================
# WHAT THIS MODULE COVERS
# ==========================================================================
# Togaviridae is monogeneric: ICTV 2024 places only Alphavirus in it, Rubivirus
# having moved to Matonaviridae. Every alphavirus has ONE genome organisation,
#
#   5'-[ nsP1 - nsP2 - nsP3 - nsP4 ]---(junction)---[ C - E3 - E2 - 6K - E1 ]-3'
#        non-structural polyprotein P1234            structural polyprotein
#
# so by the partitioning rule -- a module is a group of taxa sharing a gene
# layout -- Togaviridae is exactly one module. The BV-BRC genus column agrees:
# 19,147 of 19,180 genomes are "Alphavirus" and the remaining 33 are blank.
#
# Salmonid alphavirus (SAV/SPDV, 475 genomes) is the divergent end of the
# genus and was checked separately rather than assumed: it carries the same
# nine proteins in the same order, so it stays in this module and gets its own
# rep contig instead of its own module. Routing distance is a rep-contig
# problem, not a partitioning one.
# ==========================================================================

# --------------------------------------------------------------------------
# SPECIES DELIBERATELY NOT MODELLED
# --------------------------------------------------------------------------
# Keyed on the BV-BRC `species` field. Everything not named here is treated as
# an alphavirus, so a genuinely new alphavirus in a later export is included by
# default and only a deliberate act excludes one.
#
# Three groups, all found by profiling every species for genome length, record
# count and annotation vocabulary (species_profile.py):
#
#  1. SEGMENTED ORBIVIRUS-LIKE. Eighteen "viruses" that are 10 BV-BRC records
#     each of 3.4-4.0 kb, annotated VP1-VP7 plus NS1-NS4. That is a
#     Sedoreoviridae genome, not a togavirus -- these are Kemerovo / Great
#     Island group orbiviruses, historically parked under Togaviridae and still
#     carrying the family label here. Not one of them has a capsid, E1, E2, E3
#     or 6K annotation anywhere in the dump.
#
#  2. SEGMENTED BUNYAVIRUS-LIKE. Eight viruses that are 3 records each (L, M
#     and S, visible in the record names) annotated N / M / L protein and in
#     five cases NSs. That is a Peribunyaviridae genome. Same reasoning.
#
#  3. ALPHAVIRUS-*LIKE* SUPERFAMILY. Fungal and invertebrate viruses whose only
#     protein annotation is "RNA-dependent RNA polymerase", "polyprotein",
#     "replicase" or "RdRp" -- the alphavirus-like superfamily replicase
#     (methyltransferase-helicase-RdRp), which is a domain architecture shared
#     far outside Togaviridae. None carries the alphavirus structural cassette,
#     none exceeds 2 genomes, and most are single metagenomic contigs.
#
# Reversing any of these is one line, and the sequences are recorded in
# NOT_MODELLED.tsv rather than discarded, so the decision stays auditable.
NOT_MODELLED = {}
for _s in ("Cape Wrath virus", "Arkonam virus", "Yaquina Head virus",
           "Andasibe virus", "Bauline virus", "Gomoka virus",
           "Jacareacanga virus", "Kindia virus", "Mitchell River virus",
           "Monte Dourado virus", "Mykines virus", "Paroo River virus",
           "Seletar virus", "Sixgun city virus", "Tindholmur virus",
           "Llano Seco virus", "Mono Lake virus", "Huacho virus"):
    NOT_MODELLED[_s] = ("10 segments of 3.4-4.0 kb annotated VP1-VP7/NS1-NS4; "
                        "an orbivirus, not a togavirus")
for _s in ("Virgin River virus", "Arboledas virus", "Bobia virus",
           "Brus Laguna virus", "Kununurra virus", "Lake Clarendon virus",
           "Palestina virus", "Tsuruse virus"):
    NOT_MODELLED[_s] = ("3 segments (L/M/S) annotated N/M/L protein +/- NSs; "
                        "a bunyavirus, not a togavirus")
for _s in ("Rhizoctonia cerealis alphavirus-like virus",
           "Rhizoctonia solani alphavirus-like virus 1",
           "Rhizoctonia solani alphavirus-like 2",
           "Rhizoctonia solani alphavirus-like 3",
           "Rhizoctonia solani alphavirus-like 4",
           "Rhizoctonia solani alpha-like virus 5",
           "Rhizoctonia solani alpha-like virus 6",
           "Sclerotium rolfsii alphavirus-like virus 2",
           "Pleurotus ostreatus alpha-like virus",
           "Heterobasidion alpha-like virus 1",
           "Fusarium oxysporum f. sp. cubense alphavirus-like virus",
           "Plasmopara viticola lesion associated alpha-like virus 1",
           "Centipede toga-like virus", "Centipede toga-like virus 2",
           "Centipede toga-like virus 3", "Pseudoscorpian toga-like virus",
           "Mastotermes darwiniensis toga-like virus 1",
           "Grapevine toga-like virus", "Hattiesburg toga-like virus",
           "Shenzhen toga-like virus",
           "Flumine toga-like virus 1", "Flumine toga-like virus 2",
           "Flumine toga-like virus 3", "Flumine toga-like virus 4",
           "Flumine toga-like virus 5", "Flumine toga-like virus 6",
           "Togaviridae sp.", "Togaviridae sp. PSNE-Toga", "Alphavirus sp."):
    NOT_MODELLED[_s] = ("alphavirus-LIKE superfamily replicase only; no capsid, "
                        "E1, E2, E3 or 6K anywhere in the dump; 1-2 genomes")
# 14 records under taxon 11019 "Alphavirus" of 15-34 bp with no features at all.
NOT_MODELLED["(blank)"] = "15-34 bp records carrying no features"

# --------------------------------------------------------------------------
# FEATURES DELIBERATELY NOT BUILT
# --------------------------------------------------------------------------
# Real alphavirus products that this dump cannot train a profile for. Recorded
# here rather than declared in the JSON with nothing behind them, which is the
# failure mode the annotator cannot report.
NOT_BUILT = {
    "TF": ("transframe fusion protein: 4 sequences of 49-81 aa. Produced by a "
           "-1 ribosomal frameshift in 6K, which the annotator has no mechanism "
           "for (its `special` hooks are splice and transcript_edit only). "
           "Needs both a frameshift caller and ~50x more sequence."),
}

# --------------------------------------------------------------------------
# N-TERMINAL CURATION OVERRIDES -- feature_id -> residues to strip, with reason.
# Empty: no start-codon misassignment survived curation in this taxon. Kept so
# a future fix lands here rather than in a hand-edited fasta that the next
# rebuild erases.
# --------------------------------------------------------------------------
NTERM_TRIM = {}

# --------------------------------------------------------------------------
# TRIM A CDS COLLECTION TO ITS INITIATING MET
#
# A CDS feature's N-terminus IS an initiating Met -- that is what
# upstream_ext: 1 declares. Where a source CDS boundary is set early, the
# sequence carries leading residues that are not part of the protein, and a
# profile built from them begins upstream of the real start.
#
# That is not hypothetical here. Structural-polyprotein cluster 6 is 48 Western
# equine encephalitis records whose sequences begin
#
#     PFHEFAYESLKTRPAAPHKVPTIG | MFPYP...
#                                ^ the real start, Met in 48 of 48
#
# and `PFHEFAYESLKTRPAAPHKVPTIG` is nonstructural-polyprotein sequence -- the
# source CDS starts 25 codons early. The resulting profile aligned 25 residues
# upstream of the true start, into the junction region, where the match is noise
# and crosses a stop codon. The annotator then cropped at that stop and dropped
# the structural polyprotein ENTIRELY for those genomes: ~110 records emitted C
# and E3 but no polyprotein.
#
# So: for the CDS collections, trim each sequence forward to its first Met.
# Guarded three ways, because trimming forward to a Met is exactly the move the
# Met rule in references/curation.md warns against when applied blindly --
# it can delete real protein when the truncated form is the majority:
#   * only for keys in MET_TRIM_KEYS (CDS features, whose N-terminus is a start
#     codon; never a mature peptide, whose N-terminus is a cleavage site);
#   * the Met must lie within MET_TRIM_WINDOW residues of the start, so a
#     sequence whose first Met is deep inside is left alone;
#   * the trimmed length must still fall inside the feature's EXPECTED range.
# Anything failing a guard is left exactly as it was and reported.
# --------------------------------------------------------------------------
MET_TRIM_KEYS = {"STP", "NSP1234"}
MET_TRIM_WINDOW = 60

# --------------------------------------------------------------------------
# PER-SEQUENCE BINNING OVERRIDES -- feature_id -> (key or None to drop, reason).
# For source annotations that are simply wrong about which protein this is and
# that no regex can fix because the string is correct elsewhere.
# Populated from qc_cross_feature.py; see BUILD_NOTES.md.
# --------------------------------------------------------------------------
MISBINNED = {}

# ==========================================================================
# BINNING RULES  (feature_key, regex)  -- matched IN ORDER
# ==========================================================================
# Derived from the actual census in this dump (synonyms.tsv), not guessed.
#
# Order matters in two places and both are deliberate:
#
#  * the four nsP keys come BEFORE NSP1234, because "non-structural polyprotein
#    precursor P1234" and "RNA-directed RNA polymerase nsp4" both contain
#    "polyprotein"-adjacent words and the specific one must win.
#  * SP (structural polyprotein) comes AFTER the five structural mature
#    peptides for the same reason.
#
# What is deliberately NOT here, and why:
#
#  * `polyprotein` alone. In a rhabdovirus that word can only be L. In an
#    alphavirus it is either of the two polyproteins and the word does not say
#    which, so it goes to homology rescue instead of to a guess. (After the
#    species exclusions, 17 features carry it and most were alpha-like fungal
#    viruses.)
#  * `nonstructural protein 1`, `NS1`, `nonstructural polyprotein 4`,
#    `non-structural protein 4`, `NS5 polyprotein`. These are position labels,
#    the U-number trap wearing an alphavirus hat: the 15 CHIKV records saying
#    "nonstructural protein 1" are 136-137 aa and a real nsP1 is 535 aa, the 4
#    "nonstructural protein 4" are 60-65 aa against a 610 aa nsP4, and the 4
#    "NS5 polyprotein" are 87-88 aa in a virus that has no NS5. Every one is a
#    fragment mislabelled by ordinal. Binning them by number would train four
#    profiles on scraps.
#  * `envelope protein`, `envelope glycoprotein`, `surface glycoprotein`,
#    `envelope`, `membrane protein`. Ambiguous between E1 and E2 with no number
#    to settle it; homology rescue places them at 80% or leaves them out.
#  * `E2-E1 glycoprotein` (13, 88-218 aa) and `E3-E2 envelope protein`
#    (7, 110-247 aa). Named as precursors, sized as partials. Left unbinned;
#    see NOT_BUILT["PE2"].
# ==========================================================================
RULES = [
    # ---- non-structural mature peptides -------------------------------
    # nsP1: guanylyltransferase / methyltransferase. "methyltransferase"
    # is safe here because nsP1 is the only one in the genome.
    ("NSP1", r"\bns[Pp]?1\b|mRNA[- ]capping enzyme|methyltransferase"),
    # nsP2: protease + helicase. `protease` alone is safe -- the capsid is
    # also an autoprotease but is never annotated as one in this dump.
    ("NSP2", r"\bns[Pp]?2\b|protease\s+ns[Pp]?2"),
    ("NSP3", r"\bns[Pp]?3\b|non[- ]?structural protein P3|nonstructural protein P3"),
    # nsP4: the RdRp. The bare polymerase forms are admitted because EXPECTED
    # (500-700 aa) rejects anything that is really a P1234, and
    # reroute_length_outliers then re-homes it by homology -- the same
    # arrangement that rescues 2,100 aa "polymerase-associated" proteins in
    # Rhabdoviridae. "NA-directed" is a real typo in this dump.
    ("NSP4", r"\bns[Pp]?4\b|RNA[- ](directed|dependent)[- ]?RNA polymerase"
             r"|^N?A?-?directed RNA polymerase|^RNA polymerase$|^RdRp$|^replicase$"),
    # ---- non-structural polyprotein (the CDS) --------------------------
    ("NSP1234", r"non[- ]?structural polyprotein|nonstructural polyprotein"
                r"|non structural polyprotein|nonstructural protein P1(23)?4?$"
                r"|^Nonstructural protein P123$|^nonstructural protein P1$"
                r"|p2[37]0 nonstructural polyprotein|^Polyprotein 1$"),
    # ---- structural mature peptides ------------------------------------
    ("C",  r"capsid|^C protein$|^C$|^CP$|capside"),
    ("E3", r"\bE3\b"),
    ("E2", r"\bE2\b|\bPE2\b|\bp62\b"),
    ("6K", r"\b6K\b"),
    ("E1", r"\bE1\b|envelope (glyco)?protein 1$"),
    # ---- structural polyprotein (the CDS) ------------------------------
    ("STP", r"structural polyprotein|structural protein precus|structrual"
           r"|^structural protein$|^Polyprotein 2$"),
]
COMPILED = [(k, re.compile(p, re.I)) for k, p in RULES]

# --------------------------------------------------------------------------
# ANNOTATIONS THAT DECLARE THEMSELVES PARTIAL
#
# BV-BRC's alphavirus vocabulary marks incomplete products explicitly --
# "putative E1 envelope glycoprotein, fragment", "..., N-terminal",
# "..., C-terminal", "truncated polyprotein", "nonfunctional ... due to
# mutation". That is better evidence than any length window, and it catches the
# case a length window cannot: a declared C-terminal piece of nsP3 that happens
# to be 450 aa sits comfortably inside any range wide enough to hold Wenling
# fish alphavirus nsP3 at 421 aa.
#
# Before this rule existed, 89 Ross River records annotated "putative
# non-structural polyprotein precursor P1234, C-terminal" were being re-homed
# into the nsP3 collection at 98% identity -- they are polyprotein fragments
# that happen to span nsP3, and they would have taught the nsP3 profile two
# wrong termini.
#
# Declared partials are binned and written to <FEATURE>.partial.fasta, so they
# stay visible and countable, and never train a profile.
# --------------------------------------------------------------------------
PARTIAL_ANN = re.compile(
    r",\s*(fragment|N[- ]terminal|C[- ]terminal)\s*$"
    r"|^truncated\b|^Partial\b|nonfunctional\b", re.I)

# --------------------------------------------------------------------------
# Expected full-length ranges in amino acids. A sequence outside its range is
# NOT discarded -- it goes to <FEATURE>.outliers.fasta so it stays visible.
#
# Centres are the measured medians of the in-dump collections (CHIKV S27 as
# the yardstick: nsP1 535, nsP2 798, nsP3 530, nsP4 611, C 261, E3 64, E2 423,
# 6K 61, E1 439; P1234 2474, SP 1248). Bounds are widened to take in Salmonid
# alphavirus at the divergent end and Aura/Trocara at the long end, because a
# range fitted to CHIKV alone would park half the genus in .outliers.
# --------------------------------------------------------------------------
EXPECTED = {
    # bounds = (lowest value seen in any near-complete genome of a divergent
    #           species, rounded down) .. (highest, rounded up)
    "NSP1234": (2300, 2700),   # P1234 2321-2636. The 1813-1889 aa entries are
                               # P123 -- a real product of the un-read-through
                               # opal codon, not a partial -- but only 3 records
                               # name it, and a P123 profile could not be told
                               # from a P1234 one anyway. They land in .outliers.
    "NSP1":    (450,  600),    # 531-563 across the divergent species
    "NSP2":    (730,  880),    # 795 (seal) - 858 (salmonid)
    "NSP3":    (410,  650),    # 421 (Wenling fish) - 622 (Aura). The widest
                               # feature: nsP3's C-terminal hypervariable domain
                               # is what varies.
    "NSP4":    (560,  650),    # 605-611 everywhere; the tightest feature
    "STP":     (1150, 1400),   # 1241 (Eilat) - 1350 (Comber); salmonid 1319-1322
    "C":       (220,  320),    # 236 (salmonid) - 288
    "E3":      (45,   80),     # 53 (salmonid) - 71
    "E2":      (380,  470),    # 415 (salmonid) - 440 (Wenling fish)
    "6K":      (45,   80),     # 54 (Aura) - 68 (salmonid)
    "E1":      (400,  480),    # 431 (Caaingua) - 462 (salmonid)
    # --- multi-product intermediates, recovered from their flanking calls ---
    # Neither has a usable source collection (8 pE2 sequences in the whole
    # export, one of them full length; 2 for P123), so both come entirely from
    # Extracted/ via extract_combos.py.
    "PE2":     (380, 560),    # E3+E2; observed 382-506, median 487
    "P123":    (1650, 2100),  # nsP1+nsP2+nsP3; observed 1679-1992, median 1859
}
CORE_KEYS = ["NSP1234", "P123", "NSP1", "NSP2", "NSP3", "NSP4",
             "STP", "PE2", "C", "E3", "E2", "6K", "E1"]

# --------------------------------------------------------------------------
# Two homology thresholds, not one, and the difference matters in a polyprotein
# taxon.
#
# RESCUE_* adopts a protein no rule named into the feature it is homologous to.
# A partial submission should still be adopted and then land in <FEAT>.outliers
# on length, so query coverage stays at the documented 0.60.
#
# REHOME_* is the stricter test used when a rule DID name a protein and only its
# length disputes that. Here the claim being made is "this is not the feature it
# says, it is that one", and in Togaviridae that claim is easy to get wrong for
# a reason specific to polyproteins: every mature peptide is literally a
# substring of its precursor, so a 450 aa fragment of P1234 blasts to nsP3 at
# 98% identity over 100% of the fragment and nothing about identity or query
# coverage can tell it apart from a truncated nsP3. SUBJECT coverage can. A real
# nsP3 covers the nsP3 subject end to end; a slice of a polyprotein does not.
#
# Hence 0.80 on both sides for re-homing: the sequence must look like the whole
# of the other protein, not merely like part of it. Before this was added, 89 of
# the 100 re-homings were polyprotein fragments being promoted to mature
# peptides they would then have contaminated the alignments of.
# --------------------------------------------------------------------------
RESCUE_PIDENT = 80.0
RESCUE_COV = 0.60
RESCUE_SCOV = 0.50
REHOME_PIDENT = 80.0
REHOME_COV = 0.80
REHOME_SCOV = 0.80


def classify(ann):
    """Map an annotation string to a feature key.

    Matched twice: the string exactly as BV-BRC supplies it, then a normalised
    form with hyphens and underscores turned into spaces and whitespace
    collapsed. The second pass is strictly additive. It earns its place here on
    `non-structural`/`non structural` and on `nsP-1` style spellings.
    """
    if not ann:
        return None
    for key, prx in COMPILED:
        if prx.search(ann):
            return key
    norm = re.sub(r"[-_]", " ", ann)
    norm = re.sub(r"\s+", " ", norm).strip()
    if norm != ann:
        for key, prx in COMPILED:
            if prx.search(norm):
                return key
    return None


def _blast_best(queries, subjects_by_key, pident, cov, scov=0.0):
    """blastp every query against the named collections; return best hit per query.

    queries: list of sequences. subjects_by_key: {key: [seq, ...]}.
    `cov` is coverage of the QUERY, `scov` coverage of the SUBJECT.
    Returns {index: (key, pident, qcov, scov, bits)} for hits clearing all three.
    """
    if not queries or not shutil.which("blastp"):
        return {}
    tmp = tempfile.mkdtemp(prefix="toga.blast.")
    try:
        db = os.path.join(tmp, "named.faa")
        with open(db, "w") as fh:
            for k, seqs in sorted(subjects_by_key.items()):
                for i, s in enumerate(seqs):
                    fh.write(">%s@%d\n%s\n" % (k, i, s))
        subprocess.run(["makeblastdb", "-in", db, "-dbtype", "prot",
                        "-out", os.path.join(tmp, "db")], capture_output=True)
        q = os.path.join(tmp, "q.faa")
        with open(q, "w") as fh:
            for i, s in enumerate(queries):
                fh.write(">q%d\n%s\n" % (i, s))
        r = subprocess.run(
            ["blastp", "-query", q, "-db", os.path.join(tmp, "db"), "-outfmt",
             "6 qseqid sseqid pident length qlen slen bitscore", "-evalue", "1e-5",
             # 50, not the usual 5. In a polyprotein taxon the top hits for a
             # mature peptide are its own precursors -- 2,959 P1234 and 2,286
             # structural polyproteins sit in these collections -- and five
             # slots fill with them before a single E1 appears.
             "-max_target_seqs", "50", "-num_threads", "4"],
            capture_output=True, text=True)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    # Filter FIRST, then take the best of what qualifies -- not the other way
    # round. Taking the best hit and then thresholding it loses every rescue in
    # a polyprotein taxon: a 442 aa E1 aligns to a 1,248 aa structural
    # polyprotein and to a 442 aa E1 with the SAME bitscore, because it is the
    # same alignment, so which one is "best" is decided by tie-ordering. When
    # the precursor won, subject coverage came out at 442/1248 = 0.35, the hit
    # was dropped, and the E1 hit that would have passed had already been
    # discarded. Nine Everglades E1 proteins were lost that way, silently.
    best = {}
    for line in r.stdout.splitlines():
        qid, sid, pid, alen, qlen, slen, bits = line.split("\t")
        i = int(qid[1:])
        qcov, scov_ = int(alen) / max(int(qlen), 1), int(alen) / max(int(slen), 1)
        if float(pid) < pident or qcov < cov or scov_ < scov:
            continue
        cand = (sid.split("@")[0], float(pid), qcov, scov_, float(bits))
        if i not in best or cand[4] > best[i][4]:
            best[i] = cand
    return best


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    os.chdir(here)

    # ---- load ---------------------------------------------------------
    meta = {}
    for line in open("Toga.genome_meta.tsv"):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 7:
            meta[p[0]] = dict(name=p[2], length=int(p[3] or 0),
                              species=(p[6].strip() or "(blank)"))
    seq = {}
    for line in open("Toga.uniq.seq"):
        m, s = line.rstrip("\n").split("\t")
        seq[m] = s.strip().upper()
    ann_by_md5 = {}
    md5s = [l.strip() for l in open("Toga.uniq.md5") if l.strip()]
    for md5, line in zip(md5s, open("Toga.uniq.id_ann")):
        fid, a = (line.rstrip("\n").split("\t") + [""])[:2]
        ann_by_md5[md5] = (fid, a)

    # features carrying each md5, with their genome
    feats_by_md5 = collections.defaultdict(list)
    ftype_by_fid = {}
    for line in open("Toga.feat_full.tsv"):
        p = line.rstrip("\n").split("\t")
        if len(p) < 4 or not p[2].strip():
            continue
        g = re.match(r"^fig\|(\d+\.\d+)\.", p[0])
        if not g:
            continue
        feats_by_md5[p[2]].append((p[0], g.group(1), p[1]))
        ftype_by_fid[p[0]] = p[3]

    # ---- bin ----------------------------------------------------------
    bins = collections.defaultdict(dict)          # key -> {seq: header}
    partial = collections.defaultdict(dict)       # key -> {seq: header}, declared partial
    n_partial = 0
    syn = collections.Counter()                   # (ann, species) -> nfeat
    syn_bin = {}                                  # (ann, species) -> key
    unassigned = []                               # (md5, seq, fid, sp, ann, n)
    excluded = collections.Counter()              # species -> nfeat
    dropped_short, dropped_ambig = 0, 0
    met_trimmed = collections.Counter(); met_trim_skip = collections.Counter()
    met_trim_n = collections.defaultdict(list)

    for md5 in md5s:
        s = seq.get(md5, "")
        recs = feats_by_md5.get(md5, [])
        if not recs:
            continue
        fid, ann = ann_by_md5.get(md5, (recs[0][0], recs[0][2]))
        species = collections.Counter(meta.get(g, {}).get("species", "(blank)")
                                      for _f, g, _a in recs)
        # a sequence counts as modelled if ANY genome carrying it is one
        modelled = [sp for sp in species if sp not in NOT_MODELLED]
        sp_label = (sorted(modelled)[0] if modelled
                    else species.most_common(1)[0][0])
        nfeat = len(recs)
        syn[(ann, sp_label)] += nfeat
        if not modelled:
            excluded[sp_label] += nfeat
            syn_bin[(ann, sp_label)] = "NOT_MODELLED"
            continue

        for f, _g, _a in recs:
            if f in NTERM_TRIM:
                s = s[NTERM_TRIM[f]:]
        key = classify(ann)
        for f, _g, _a in recs:
            if f in MISBINNED:
                key = MISBINNED[f][0]
        syn_bin[(ann, sp_label)] = key or "UNASSIGNED"
        if key is None:
            if len(s.strip("X")) >= MIN_AA and not (set(s) & HARD_AMBIG) \
                    and s.strip("X").count("X") / max(len(s.strip("X")), 1) <= MAX_AMBIG_FRAC:
                unassigned.append((md5, s.strip("X"), fid, sp_label, ann, nfeat))
            continue
        s = s.strip("X")
        if len(s) < MIN_AA:
            dropped_short += 1
            continue
        if (set(s) & HARD_AMBIG) or s.count("X") / len(s) > MAX_AMBIG_FRAC:
            dropped_ambig += 1
            continue
        if key in MET_TRIM_KEYS and not s.startswith("M"):
            m = s.find("M", 0, MET_TRIM_WINDOW)
            lo_e, hi_e = EXPECTED.get(key, (0, 10 ** 6))
            if m > 0 and lo_e <= len(s) - m <= hi_e:
                met_trimmed[key] += 1
                met_trim_n[key].append(m)
                s = s[m:]
            elif m != 0:
                met_trim_skip[key] += 1
        gname = meta.get(recs[0][1], {}).get("name", "?")
        hdr = "%s [%s] %s | %s" % (fid, gname, sp_label, ann)
        if PARTIAL_ANN.search(ann):
            partial[key].setdefault(s, hdr)
            n_partial += 1
            continue
        bins[key].setdefault(s, hdr)

    # ---- homology rescue of the unnamed --------------------------------
    named = {k: list(v) for k, v in bins.items()}
    best = _blast_best([u[1] for u in unassigned], named, RESCUE_PIDENT,
                       RESCUE_COV, RESCUE_SCOV)
    adopted, lineage = [], []
    for i, (md5, s, fid, sp, ann, n) in enumerate(unassigned):
        if i in best:
            key, pid, cov, _sc, _b = best[i]
            tgt = partial if PARTIAL_ANN.search(ann) else bins
            tgt[key].setdefault(s, "%s [rescued] %s | %s" % (fid, sp, ann))
            adopted.append((key, fid, sp, len(s), pid, cov, ann))
            # record the rescue in synonyms.tsv. Without this the string-collapse
            # artifact reports every homology-rescued string as UNBINNED, which
            # understates the collapse and, worse, lists strings the module does
            # cover as strings it does not.
            syn_bin[(ann, sp)] = key
        else:
            lineage.append((fid, sp, len(s), ann))

    # ---- merge anything recovered from genomes -------------------------
    extracted = []
    exdir = os.path.join(here, "Extracted")
    if os.path.isdir(exdir):
        for fa in sorted(os.listdir(exdir)):
            if not fa.endswith(".faa") or "." not in fa[:-4]:
                continue
            mod, key = fa[:-4].split(".", 1)
            if mod != MODULE:
                continue
            n, hdr, buf = 0, None, []
            def flush():
                nonlocal n
                if hdr and buf:
                    sq = "".join(buf)
                    if sq not in bins[key]:
                        bins[key][sq] = hdr
                        n += 1
            for line in open(os.path.join(exdir, fa), errors="replace"):
                if line.startswith(";") or not line.strip():
                    continue          # provenance header lines, not sequence
                if line.startswith(">"):
                    flush(); hdr, buf = line[1:].strip(), []
                elif line.strip():
                    buf.append(line.strip())
            flush()
            if n:
                extracted.append((key, n))

    # ---- length gate, with homology re-homing for the strays -----------
    stray = []
    for k in list(bins):
        lo, hi = EXPECTED.get(k, (0, 10 ** 6))
        for s, h in bins[k].items():
            if not (lo <= len(s) <= hi):
                stray.append((k, s, h))
    inrange = {k: [s for s in v if EXPECTED.get(k, (0, 10**6))[0] <= len(s)
                   <= EXPECTED.get(k, (0, 10**6))[1]] for k, v in bins.items()}
    best = _blast_best([x[1] for x in stray], {k: v for k, v in inrange.items() if v},
                       REHOME_PIDENT, REHOME_COV, REHOME_SCOV)
    moved = []
    for i, (k, s, h) in enumerate(stray):
        if i not in best:
            continue
        tgt, pid, cov, _sc, _b = best[i]
        if tgt == k:
            continue
        lo, hi = EXPECTED.get(tgt, (0, 10 ** 6))
        if not (lo <= len(s) <= hi):
            continue
        del bins[k][s]
        bins[tgt].setdefault(s, h)
        moved.append((k, tgt, len(s), pid, h.split()[0]))

    # ---- write ---------------------------------------------------------
    outdir = os.path.join(here, "collections", MODULE)
    if os.path.isdir(os.path.join(here, "collections")):
        shutil.rmtree(os.path.join(here, "collections"))
    os.makedirs(outdir)
    report = []
    for k in sorted(bins, key=lambda x: (CORE_KEYS.index(x) if x in CORE_KEYS else 99, x)):
        lo, hi = EXPECTED.get(k, (0, 10 ** 6))
        keep = {s: h for s, h in bins[k].items() if lo <= len(s) <= hi}
        out = {s: h for s, h in bins[k].items() if not (lo <= len(s) <= hi)}
        for fn, d in ((k + ".fasta", keep), (k + ".outliers.fasta", out),
                      (k + ".partial.fasta", partial.get(k, {}))):
            if not d:
                continue
            with open(os.path.join(outdir, fn), "w") as fh:
                for s, h in sorted(d.items(), key=lambda kv: (-len(kv[0]), kv[1])):
                    fh.write(">%s\n%s\n" % (h, s))
        L = sorted(len(s) for s in keep) or [0]
        report.append((k, len(keep), len(out), len(partial.get(k, {})),
                       L[0], L[len(L) // 2], L[-1], lo, hi))

    with open(os.path.join(here, "collections", "report.tsv"), "w") as fh:
        fh.write("feature\tn_kept\tn_outliers\tn_declared_partial\tmin_aa\tmed_aa\tmax_aa\texp_lo\texp_hi\n")
        for r in report:
            fh.write("\t".join(str(x) for x in r) + "\n")

    with open(os.path.join(here, "synonyms.tsv"), "w") as fh:
        fh.write("count\tannotation\tspecies\tbins_to\n")
        for (a, sp), n in syn.most_common():
            fh.write("%d\t%s\t%s\t%s\n" % (n, a, sp, syn_bin.get((a, sp), "?")))

    with open(os.path.join(here, "collections", "UNASSIGNED_TRACKING.tsv"), "w") as fh:
        fh.write("feature_id\tspecies\taa\tannotation\tdisposition\n")
        for key, fid, sp, n, pid, cov, a in sorted(adopted, key=lambda x: x[0]):
            fh.write("%s\t%s\t%d\t%s\trescued to %s at %.1f%% id / %.0f%% cov\n"
                     % (fid, sp, n, a, key, pid, cov * 100))
        for fid, sp, n, a in sorted(lineage, key=lambda x: -x[2]):
            fh.write("%s\t%s\t%d\t%s\tleft out: no named feature at >=%.0f%% id\n"
                     % (fid, sp, n, a, RESCUE_PIDENT))

    with open(os.path.join(here, "collections", "NOT_MODELLED.tsv"), "w") as fh:
        fh.write("kind\tname\tn_features\treason\n")
        for sp, n in excluded.most_common():
            fh.write("species\t%s\t%d\t%s\n" % (sp, n, NOT_MODELLED[sp]))
        for k, why in sorted(NOT_BUILT.items()):
            fh.write("feature\t%s\t\t%s\n" % (k, why))

    # near-complete genomes missing each core protein
    complete = [g for g, d in meta.items()
                if d["length"] >= 10000 and d["species"] not in NOT_MODELLED]
    have = collections.defaultdict(set)
    for md5, recs in feats_by_md5.items():
        s = seq.get(md5, "")
        for f, g, a in recs:
            k = classify(a)
            if k:
                have[k].add(g)
    with open(os.path.join(here, "collections", "MISSING.tsv"), "w") as fh:
        fh.write("feature\tn_complete_genomes\tn_with\tn_missing\n")
        for k in CORE_KEYS:
            w = len([g for g in complete if g in have.get(k, ())])
            fh.write("%s\t%d\t%d\t%d\n" % (k, len(complete), w, len(complete) - w))

    # ---- console summary ----------------------------------------------
    print("module: %s" % MODULE)
    print("  unique sequences in dump   : %d" % len(md5s))
    print("  excluded (not modelled)    : %d features over %d species"
          % (sum(excluded.values()), len(excluded)))
    print("  dropped, <%d aa            : %d" % (MIN_AA, dropped_short))
    print("  dropped, >%.0f%% ambiguous   : %d" % (MAX_AMBIG_FRAC * 100, dropped_ambig))
    print("  declared partial (not used): %d" % n_partial)
    for k in sorted(met_trimmed):
        n = sorted(met_trim_n[k])
        print("  %-8s trimmed to the initiating Met: %d sequence(s), median %d "
              "residue(s) removed (%d left alone by a guard)"
              % (k, met_trimmed[k], n[len(n) // 2], met_trim_skip[k]))
    print("  rescued by homology        : %d" % len(adopted))
    print("  left out (no rule, no hit) : %d" % len(lineage))
    print("  re-homed by length+homology: %d" % len(moved))
    if extracted:
        print("  merged from Extracted/     : %s" % extracted)
    print()
    print("%-9s %7s %9s %8s %7s %7s %7s   %s" %
          ("feature", "kept", "outliers", "partial", "min", "med", "max", "expected"))
    for k, nk, no, npar, lo_, med, hi_, elo, ehi in report:
        print("%-9s %7d %9d %8d %7d %7d %7d   %d-%d"
              % (k, nk, no, npar, lo_, med, hi_, elo, ehi))
    for k, why in sorted(NOT_BUILT.items()):
        print("%-9s %7s %9s %8s   %s" % (k, "-", "-", "-", why.split(":")[0]))
    if moved:
        print("\nre-homed:")
        for a, b, n, pid, fid in moved[:20]:
            print("  %-8s -> %-8s %5d aa  %.1f%%  %s" % (a, b, n, pid, fid))
        if len(moved) > 20:
            print("  ... %d more" % (len(moved) - 20))


if __name__ == "__main__":
    main()
