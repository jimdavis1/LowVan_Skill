#!/usr/bin/env python3
"""Write the Orthoflavivirus module block.

The genome is a single ORF translated as one polyprotein and cleaved into ten
mature products, so the module follows the Togaviridae pattern: a polyprotein
CDS plus mat_peptides, rather than pretending each product is its own gene.

    C - prM/M - E - NS1 - NS2A - NS2B - NS3 - NS4A - 2K - NS4B - NS5

Extension flags follow the rule in references/json-schema.md: only the first
product of a precursor gets `upstream_ext: 1`, and every mature peptide gets
`downstream_ext: 0`, because a mature peptide never includes the stop codon
even when it runs to the protein's real end. The two products that begin at
the initiating Met are the capsid forms; everything downstream begins at a
cleavage site, which is neither a start codon nor a stop.

Note that Togaviridae diverges from that rule -- its E1 and NSP4 carry
`downstream_ext: 1`. This module follows the documented rule rather than the
older module.

NS1' IS DELIBERATELY ABSENT. It is real: 533 features, median 404 aa against
NS1's 352, which is exactly the 52-residue extension produced by the -1
ribosomal frameshift near the start of NS2A. Its carriers in this dump are
only the JEV serogroup -- West Nile 111, JEV 4, Usutu 3, Murray Valley 2,
Cacipacore 1 -- with zero dengue, Zika, yellow fever or TBE despite dengue
holding 4,121 ordinary NS1 sequences. Both published facts are reproduced
from the sequence data alone.

A PSSM cannot call it, because NS1' shares its first 352 residues with NS1
and a profile would simply call NS1. Per step 6b it needs
`special: transcript_edit` with a hand-curated nucleotide reference set,
which is a build cycle of its own. The collection is retained at
collections/Orthoflavivirus/NS1P.fasta so that work can start from it.
"""
import glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MODULE = "Orthoflavivirus"

#  Citations proposed by a model, never read by a curator. Disjoint from PMID.
CLAUDE_PMIDS = {
    "NS3": ["8107205", "33494395"],
    "NS5": ["2170676"],
    "PRM": ["18184701"],
    "E": ["7753193"],
    #  The three added here are the basis of the junction filter in
    #  build_collections.py, which is why they are recorded against the
    #  features whose boundaries they define:
    #    8445732   Lin, Amberg, Chambers & Rice, J Virol 1993;67(4):2327-2335.
    #              The NS4A/2K site is cleaved by NS2B-NS3 and sits exactly
    #              23 residues upstream of the 2K/NS4B signalase site. The
    #              cleaned 2K collection measures 23 at p10, median and p90,
    #              which reproduces this from sequence alone -- the filter
    #              knows only residue classes at the two ends, never a length.
    #    33494395  Wahaab et al., Pathogens 2021;10(2):102. NS2B-NS3 requires
    #              two basic residues (K-R, R-R, R-K, occasionally Q-R) at
    #              P2-P1 and a small residue (G, S, A) at P1'.
    #    2174669   Chambers, Hahn, Galler & Rice, Annu Rev Microbiol
    #              1990;44:649-688. The polyprotein cleavage map.
    "2K":    ["8445732"],
    "NS4A":  ["8445732"],
    "NS4B":  ["8445732"],
    "NS2B":  ["33494395"],
    "POLY":  ["2174669"],
}

#  min_len / max_len are the observed range of each CLEANED collection, not
#  a hand-drawn window. The originals were drawn from the contaminated
#  distributions and were far too wide: 2K was declared 15-193 wide enough to
#  admit an 18-residue chunk of the NS3 helicase, and NS4A 103-307 for a
#  protein that measures 119-161. Regenerate this table by re-running
#  build_collections.py and reading observed_bounds.json.
#
#  key: (annotation, symbol, type, min, max, upstream_ext, downstream_ext)
F = [
    #  Strings follow the vocabulary's mat_peptide convention: a cleaved
    #  product reads "Mature ...", a precursor keeps the precursor wording
    #  (as Togaviridae's "Envelope glycoprotein precursor pE2" does), and a
    #  signal peptide reads "Signal peptide of <symbol>". 2K takes the last
    #  form because that is literally what it is -- the signal peptide that
    #  targets NS4B to the ER membrane.
    #  Two strings are reused rather than invented: "Nucleocapsid protein",
    #  which Togaviridae already uses for a mat_peptide capsid, and
    #  "RNA-dependent RNA polymerase", which it already uses for nsP4.
    ("POLY",  "Genome polyprotein",                        "POLY",  "CDS",         2686, 4237, 1, 1),
    ("ANCHC", "Anchored capsid protein precursor ancC",    "ancC",  "mat_peptide",   80, 143, 1, 0),
    ("C",     "Nucleocapsid protein",                      "C",     "mat_peptide",   82, 124, 1, 0),
    ("PRM",   "Membrane glycoprotein precursor prM",       "prM",   "mat_peptide",  128, 185, 0, 0),
    #  pr is the one feature whose bounds must NOT come from the collection.
    #  BV-BRC annotates no pr peptide for the tick-borne clade at all -- the
    #  632 training sequences run 69-97 aa with 500 of them at exactly 91,
    #  all mosquito-borne. The profiles nonetheless call pr correctly in
    #  tick-borne genomes, at its true and shorter boundary: TBE 244 calls
    #  median 59, Kyasanur 164 at 59, Omsk 51 at 59, louping ill 30 at 59,
    #  all >=99% below the collection minimum. Tick-borne prM is shorter
    #  than mosquito-borne prM, so their pr genuinely is ~59 aa. A min_len
    #  of 69 flags a correct call on an entire clade, so the bound is set
    #  from what the feature IS across the taxon, not from what the
    #  collection happens to contain.
    ("PR",    "Mature peptide pr",                         "pr",    "mat_peptide",   55, 97, 0, 0),
    ("M",     "Mature membrane glycoprotein M",            "M",     "mat_peptide",   54, 76, 0, 0),
    ("E",     "Mature envelope glycoprotein E",            "E",     "mat_peptide",  352, 505, 0, 0),
    ("NS1",   "Mature non-structural protein NS1",         "NS1",   "mat_peptide",  263, 415, 0, 0),
    ("NS2A",  "Mature non-structural protein NS2A",        "NS2A",  "mat_peptide",  157, 239, 0, 0),
    ("NS2B",  "Mature protease cofactor NS2B protein",     "NS2B",  "mat_peptide",  118, 149, 0, 0),
    ("NS3",   "Mature protease and helicase NS3 protein",  "NS3",   "mat_peptide",  530, 623, 0, 0),
    ("NS4A",  "Mature non-structural protein NS4A",        "NS4A",  "mat_peptide",  119, 161, 0, 0),
    ("2K",    "Signal peptide of NS4B",                    "2K",    "mat_peptide",   19, 27, 0, 0),
    ("NS4B",  "Mature non-structural protein NS4B",        "NS4B",  "mat_peptide",  222, 268, 0, 0),
    ("NS5",   "RNA-dependent RNA polymerase",              "NS5",   "mat_peptide",  642, 1099, 0, 0),
]

#  bit_cutoff scales with length: a 23-residue 2K cannot clear the bar a
#  900-residue NS5 clears, and a high bar on a short product means the
#  feature is simply never called.
BIT = {"POLY": 400, "NS5": 400, "NS3": 300, "E": 200, "NS1": 150, "NS4B": 100,
       "NS2A": 80, "ANCHC": 60, "C": 60, "PRM": 60, "NS2B": 60, "NS4A": 60,
       "M": 40, "PR": 40, "2K": 25}


def close_genomes():
    p = os.path.join(HERE, "Rep-Contigs", "close_genomes.json")
    if os.path.exists(p):
        return json.load(open(p))
    out = {}
    for f in sorted(glob.glob(os.path.join(HERE, "Rep-Contigs", "*.dna"))):
        hdr = open(f).readline().lstrip(">").strip()
        acc = hdr.split()[0]
        out[os.path.basename(f)] = {"genome_ids": acc,
                                    "genome_name": " ".join(hdr.split()[1:]) or acc}
    return out


def main():
    feats = {}
    for key, anno, sym, ftype, mn, mx, up, down in F:
        blk = {
            "anno": anno,
            "gene_symbol": sym,
            "feature_type": ftype,
            "type": ftype,
            "bit_cutoff": BIT[key],
            "coverage_cutoff": 0.65,
            "min_len": mn,
            "max_len": mx,
            "upstream_ext": up,
            "downstream_ext": down,
            "kmers": 1,
            "copy_num": 1,
            "segment": "Single RNA Segment",
        }
        if key in CLAUDE_PMIDS:
            blk["PMID_claude_generated"] = CLAUDE_PMIDS[key]
        feats[key] = blk

    block = {MODULE: {
        "close_genomes": close_genomes(),
        "segments": {"Single RNA Segment":
                     {"min_len": 9000, "max_len": 13000, "replicon_geometry": "linear"}},
        "features": feats,
    }}
    out = os.path.join(HERE, "%s_Viral_PSSM.json" % MODULE)
    with open(out, "w") as fh:
        json.dump(block, fh, indent=3, sort_keys=True, separators=(",", " : "))
        fh.write("\n")
    print("wrote %s: %d features (%d CDS, %d mat_peptide), %d rep contigs"
          % (out, len(feats),
             sum(1 for f in feats.values() if f["feature_type"] == "CDS"),
             sum(1 for f in feats.values() if f["feature_type"] == "mat_peptide"),
             len(block[MODULE]["close_genomes"])))


if __name__ == "__main__":
    main()
