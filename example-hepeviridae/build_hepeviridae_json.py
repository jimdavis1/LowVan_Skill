#!/usr/bin/env python3
"""Generate the Hepeviridae block of Viral_PSSM.json.

Lengths are measured from the collections, not asserted. Citations come from
CLAUDE_PMIDS / READ_PMIDS so the model-proposed flag cannot drift from reality:
delete an id from CLAUDE_PMIDS only when a human has read that paper.

Four features. The family's architecture was measured across every
near-complete genome (PARTITIONING.md): ORF1 + ORF2 + ORF3 in every genus, plus
a fourth small ORF in Rocahepevirus.

ORF1 IS DECLARED WHOLE. RefSeq annotates seven domains on it -- methyltransferase,
Y, papain-like protease, poly-proline hinge, X, helicase, RdRp -- and this module
deliberately does not. That annotation descends from Koonin 1992 (PMID 1518855),
which is titled "Computer-assisted assignment of functional domains" and calls
them "putative" throughout: a sequence-comparison prediction, never a
demonstration of cleavage. The experimental record since does not support
processing. Perttila 2013 (PMID 23255617) expressed full-length pORF1 and six
truncations in HeLa and Huh-7 under several vector systems with antisera to the
MT, HEL and POL domains, and concluded "the weight of evidence supports the
proposition that pORF1 is not subjected to specific proteolytic processing".
LeDesma 2023 (PMID 36852909) scanned the putative protease domain by mutagenesis
and found "ORF1 operates as a multifunctional protein, which is not subject to
proteolytic processing" -- the six essential cysteines coordinate divalent metal
ions and hold the fold together, and the domain rescues replication only in the
context of full-length ORF1.

The seven domain boundaries, as NC_001434 carries them, are recorded here so the
decision can be revisited without re-deriving them:

    methyltransferase      1- 181      X domain            712- 868
    Y domain             182- 397      helicase            869-1112
    papain-like protease 398- 556      RdRp               1113-1598
    poly-proline hinge   557- 622

ORF4 IS SPLIT, NOT SHARED. Paslahepevirus ORF4 and Rocahepevirus ORF4 carry the
same positional label and are not homologous: over every ORF4-labelled sequence
in the dump the best cross-genus BLASTp alignment is 5 residues at e-value 138.
Only the Rocahepevirus protein is declared. The Paslahepevirus one (the HEV-1
ER-stress factor, PMID 27035822) has two unique sequences in BV-BRC and is
retired to NOT_MODELLED.tsv.
"""
import json, os, collections, sys

MODULE = "Hepeviridae"
SEGMENT = "Single RNA Segment"

# ids a model proposed. They go in PMID *and* in PMID_claude_generated until a
# human reads the paper and deletes the entry here.
#  EVERY id here is model-proposed, so every one goes in PMID_claude_generated
#  and the PMID field stays empty on all four features.
#
#  An earlier version of this generator put four ids in PMID on the grounds that
#  I had read their abstracts. That was wrong and it is worth stating why:
#  reading a paper is not the same act as a curator confirming it establishes
#  what the feature claims. A DLIT is the curator's judgement, not the model's
#  reading. The flag is cleared by a human deleting the entry from this
#  registry, never by the model deciding it has done enough homework.
#  Same convention as the Matonaviridae module, where PMID is empty throughout.
CLAUDE_PMIDS = {
    "1518855":  "Koonin et al. 1992 PNAS, Computer-assisted assignment of functional domains in the nonstructural polyprotein of hepatitis E virus -- the source of the seven-domain annotation, and explicitly a prediction",
    "23255617": "Perttila, Spuul & Ahola 2013 J Gen Virol, Early secretory pathway localization and lack of processing for hepatitis E virus replication protein pORF1",
    "36852909": "LeDesma et al. 2023 eLife, Structural features stabilized by divalent cation coordination within hepatitis E virus ORF1 are critical for viral replication",
    "28096411": "Ding et al. 2017 PNAS, Hepatitis E virus ORF3 is a functional ion channel required for release of infectious particles",
    "34523963": "Sari et al. 2021 J Virol, The viral ORF3 protein is required for hepatitis E virus apical release and efficient growth in polarized hepatocytes and humanized mice",
    "19622744": "Guu et al. 2009 PNAS, Structure of the hepatitis E virus-like particle suggests mechanisms for virus assembly and receptor binding",
    "27035822": "Nair et al. 2016 PLoS Pathog, Endoplasmic reticulum stress induced synthesis of a novel viral factor mediates efficient replication of genotype-1 hepatitis E virus (cited on the retired Paslahepevirus ORF4)",
    "31399875": "Reuter et al. 2019 Arch Virol, Genomic and spatial variability of a European common vole hepevirus",
}
#  Every id above was checked to RESOLVE TO THE PAPER NAMED in PubMed. Resolving
#  an id is not reading the paper, and reading the paper is not curating it, so
#  all of them stay flagged.
READ_PMIDS = set()

FEATURES = {
    "ORF1": dict(
        anno="Nonstructural polyprotein", symbol="ORF1", ftype="CDS",
        bit=800, cov=0.65, up=1, down=1, copy=1,
        pmid=["1518855", "23255617", "36852909"]),
    "ORF2": dict(
        anno="Nucleocapsid protein", symbol="ORF2", ftype="CDS",
        bit=250, cov=0.65, up=1, down=1, copy=1,
        pmid=["19622744"]),
    #  bit_cutoff 90 is measured, not chosen by feel, per references/json-schema.md.
    #  On the 42-genome held-out panel:
    #
    #    noise ceiling                78.6   ORF3.4 returns a SECOND HSP, query
    #                                        residues 1-64, in ORF2's frame; with
    #                                        upstream_ext/downstream_ext it grew to
    #                                        ORF2's exact span and produced the one
    #                                        duplicate call in the whole evaluation
    #    typical cross-frame noise    26-35
    #    weakest genuine ORF3 call    99.0
    #
    #  So any cutoff in 79..99 works and 90 errs toward specificity.
    #
    #  What it costs, recorded rather than hidden: ORF3 is no longer called on
    #  three Chirohepevirus genomes that scored 62.6, 77.2 and 84.9. That is a
    #  COVERAGE hole, not a threshold problem -- the ORF3 collection holds 863
    #  Paslahepevirus sequences against 2 Chirohepevirus and 0 Piscihepevirus,
    #  so no threshold makes those genera callable. The fix is sequences.
    "ORF3": dict(
        anno="Membrane-associated viroporin required for virion release ORF3 protein",
        symbol="ORF3", ftype="CDS",
        bit=90, cov=0.65, up=1, down=1, copy=1,
        pmid=["28096411", "34523963"]),
    #  Plain uncharacterized string, no designation: 35 rows in the vocabulary
    #  already share it, and Reuter 2019 describes genomic variability, not a
    #  function for this protein. The ORF4 tag stays in the feature key and, in
    #  namespaced form, in the symbol -- following Sarbeco_ORF3a / Nobeco_NS7a /
    #  D_COV_NS7b, which is how the shipped vocabulary resolves a positional
    #  label that collides across taxa (here with Alphacoronavirus ORF4).
    "ORF4_ROCA": dict(
        anno="Uncharacterized lineage-specific protein", symbol="Roca_ORF4", ftype="CDS",
        bit=60, cov=0.65, up=1, down=1,
        pmid=["31399875"]),
}

def lengths(key):
    L = []; s = None
    p = os.path.join("collections", MODULE, key + ".fasta")
    for line in open(p):
        if line.startswith(">"):
            if s is not None: L.append(len(s))
            s = ""
        else: s += line.strip()
    if s is not None: L.append(len(s))
    return sorted(L)

def main():
    close = json.load(open("Rep-Contigs/close_genomes.json"))
    feats = {}
    print("%-12s %5s %6s %6s %8s %8s" % ("feature", "n", "min", "max", "min_len", "max_len"))
    for key, d in FEATURES.items():
        L = lengths(key)
        # widen the observed range a little, as references/json-schema.md says:
        # too tight and the quality tool flags healthy genomes
        lo = int(L[0] * 0.90)
        hi = int(L[-1] * 1.10)
        print("%-12s %5d %6d %6d %8d %8d" % (key, len(L), L[0], L[-1], lo, hi))
        f = collections.OrderedDict()
        f["anno"] = d["anno"]
        f["bit_cutoff"] = d["bit"]
        if "copy" in d: f["copy_num"] = d["copy"]
        f["coverage_cutoff"] = d["cov"]
        f["downstream_ext"] = d["down"]
        f["feature_type"] = d["ftype"]
        f["gene_symbol"] = d["symbol"]
        f["kmers"] = 1
        f["max_len"] = hi
        f["min_len"] = lo
        f["segment"] = SEGMENT
        f["upstream_ext"] = d["up"]
        #  DISJOINT, not nested. references/json-schema.md says a model-proposed
        #  id goes in both lists; check_dlits.py and all 36 shipped modules say
        #  they must not overlap. The code and the installed data agree with each
        #  other, so they win; the doc is wrong and is reported in SKILL_FEEDBACK.
        ids = [str(x) for x in d["pmid"]]
        proposed = [p for p in ids if p in CLAUDE_PMIDS]
        curated = [p for p in ids if p not in CLAUDE_PMIDS]
        if curated:  f["PMID"] = curated
        if proposed: f["PMID_claude_generated"] = proposed
        unknown = [p for p in ids if p not in CLAUDE_PMIDS and p not in READ_PMIDS]
        if unknown:
            sys.exit("PMID %s is in neither CLAUDE_PMIDS nor READ_PMIDS -- "
                     "say which it is before shipping" % unknown)
        feats[key] = f

    block = {MODULE: collections.OrderedDict([
        ("close_genomes", close),
        ("features", feats),
        ("segments", {SEGMENT: {"max_len": 8000, "min_len": 5500,
                                "replicon_geometry": "linear"}}),
    ])}
    out = json.dumps(block, indent=3, separators=(",", " : "),
                     sort_keys=True, ensure_ascii=False) + "\n"
    open("Hepeviridae_Viral_PSSM.json", "w").write(out)
    print("\nwrote Hepeviridae_Viral_PSSM.json  (%d features, %d rep contigs)"
          % (len(feats), len(close)))

main()
