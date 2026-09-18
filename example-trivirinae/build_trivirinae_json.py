#!/usr/bin/env python3
"""Trivirinae module JSON.

The module is the half of Betaflexiviridae whose genomes move by a SINGLE
movement protein rather than a triple gene block: Trichovirus, Vitivirus,
Citrivirus, Tepovirus, Prunevirus, Chordovirus, Divavirus, Wamavirus. The TGB
genera are a separate module -- splitting them took routing from 75.3% lumped
to 91.6% here, which is the whole argument for partitioning on genome
organisation rather than phylogeny.

Five features, and only three of them are universal.

REP, MP and CP are called on 200/200 sampled genomes (CP 198/200). The
replicase is one uncleaved ORF carrying methyltransferase, helicase and RdRp
domains -- Betaflexiviridae encodes no protease, so "polyprotein" would be
wrong here for the same reason it was wrong for the tobamovirus replicase.

NABP and VITI_ORF2 are accessories, and their cutoffs are set from measured
cross-genus noise rather than from the score distribution alone:

  NABP        Vitivirus 77/79, Prunevirus 5/6, and nothing in Citrivirus,
              Divavirus or Wamavirus. The noise ceiling across genera that
              lack it is 39 bits (Tepovirus 39, Trichovirus 38, Prunevirus
              38); the weakest genuine Vitivirus signal is 75. bit_cutoff 60
              sits between them. One Trichovirus scores 261, far above noise
              -- that is a real NABP, and most likely a Vitivirus carrying a
              wrong genus label, which is a known hazard in this dump.

  VITI_ORF2   Vitivirus and nothing else: 56/79 Vitivirus, zero hits in any
              other genus at all. There is no cross-genus noise to exclude,
              so the cutoff is set low (50) for sensitivity.

VITI_ORF2 keeps the uncharacterized annotation string, and its gene_symbol
carries the lineage prefix -- Viti_ORF2, not ORF2 -- because a bare ORF2
already names a different protein in Hepeviridae. That is the same convention
S1 Table uses for Roca_ORF4 and Nobeco_NS7a. 235 of its 289 source
records say "hypothetical protein" and the rest give a mass ("19kDa-protein",
"22 kDa protein"); no function is established, so the positional name stays in
the feature key and gene_symbol and out of the annotation, per the U-number
rule.

NABP does take a functional string. That is a judgement worth flagging: the
only citations here are model-proposed, but the claim rests on the source
data's own consensus rather than on those papers -- 337 of 470 records call it
an RNA-binding or nucleic-acid-binding protein. If that is not the bar, the
fallback is "Uncharacterized lineage-specific protein" and a one-line change
below.
"""
import json, os, collections, sys

MODULE  = "Trivirinae"
SEGMENT = "Single RNA Segment"

#  Every id model-proposed; each was checked to resolve in PubMed to the paper
#  named, and one candidate was discarded for resolving to the Alphaflexiviridae
#  profile instead. Resolving is not reading, so all stay flagged and PMID is
#  absent on every feature.
CLAUDE_PMIDS = {
    "36399124": "Adams et al. 2022 J Gen Virol, Virus classification based on in-depth sequence analyses and development of demarcation criteria using the Betaflexiviridae as a case study",
    "15098118": "Adams et al. 2004 Arch Virol, The new plant virus family Flexiviridae and assessment of molecular criteria for species demarcation",
    "9125055":  "Minafra et al. 1997 Arch Virol, Grapevine virus A: nucleotide sequence, genome organization, and relationship in the Trichovirus genus",
    "16847135": "Zhou et al. 2006 J Gen Virol, Identification of an RNA-silencing suppressor in the genome of Grapevine virus A",
}
READ_PMIDS = set()

FEATURES = {
    "REP":       dict(anno="RNA-dependent RNA polymerase", symbol="L", ftype="CDS",
                      bit=400, cov=0.65, up=1, down=1, copy=1,
                      pmid=["36399124", "15098118"]),
    "MP":        dict(anno="Movement protein", symbol="Mov", ftype="CDS",
                      bit=50, cov=0.65, up=1, down=1, copy=1,
                      pmid=["36399124", "15098118"]),
    "CP":        dict(anno="Nucleocapsid protein", symbol="N", ftype="CDS",
                      bit=45, cov=0.65, up=1, down=1, copy=1,
                      pmid=["36399124", "15098118"]),
    "NABP":      dict(anno="Nucleic acid-binding protein", symbol="NABP", ftype="CDS",
                      bit=60, cov=0.65, up=1, down=1,
                      pmid=["9125055", "16847135"]),
    "VITI_ORF2": dict(anno="Uncharacterized lineage-specific protein", symbol="Viti_ORF2",
                      ftype="CDS", bit=50, cov=0.65, up=1, down=1,
                      pmid=["9125055"]),
}

#  NABP and VITI_ORF2 carry no copy_num on purpose: they are absent from whole
#  genera, and viral_genome_quality.pl treats a copy_num feature as essential
#  and flags every genome that legitimately lacks it.

def lengths(key):
    L, s = [], None
    for line in open(os.path.join("collections", MODULE, key + ".fasta")):
        if line.startswith(">"):
            if s is not None: L.append(len(s))
            s = ""
        else: s += line.strip()
    if s is not None: L.append(len(s))
    return sorted(L)

def main():
    close = json.load(open("Rep-Contigs/close_genomes.json"))
    feats = {}
    print("%-12s %6s %6s %6s %9s %9s %6s" % ("feature","n","min","max","min_len","max_len","bit"))
    for key, d in FEATURES.items():
        L = lengths(key); lo = int(L[0]*0.92); hi = int(L[-1]*1.08)
        print("%-12s %6d %6d %6d %9d %9d %6d" % (key,len(L),L[0],L[-1],lo,hi,d["bit"]))
        f = collections.OrderedDict()
        f["anno"] = d["anno"]; f["bit_cutoff"] = d["bit"]
        if "copy" in d: f["copy_num"] = d["copy"]
        f["coverage_cutoff"] = d["cov"]; f["downstream_ext"] = d["down"]
        f["feature_type"] = d["ftype"]; f["gene_symbol"] = d["symbol"]
        f["kmers"] = 1; f["max_len"] = hi; f["min_len"] = lo
        f["segment"] = SEGMENT; f["upstream_ext"] = d["up"]
        ids = [str(x) for x in d["pmid"]]
        unk = [p for p in ids if p not in CLAUDE_PMIDS and p not in READ_PMIDS]
        if unk: sys.exit("PMID %s in neither registry" % unk)
        prop = [p for p in ids if p in CLAUDE_PMIDS]
        cur  = [p for p in ids if p in READ_PMIDS]
        if cur:  f["PMID"] = cur
        if prop: f["PMID_claude_generated"] = prop
        feats[key] = f

    annos = [v["anno"] for k, v in feats.items() if v.get("copy_num")]
    if len(annos) != len(set(annos)):
        sys.exit("two essential features share an annotation string")

    block = {MODULE: collections.OrderedDict([
        ("close_genomes", close),
        ("features", feats),
        ("segments", {SEGMENT: {"max_len": 9500, "min_len": 5800,
                                "replicon_geometry": "linear"}}),
    ])}
    open("%s_Viral_PSSM.json" % MODULE, "w").write(
        json.dumps(block, indent=3, separators=(","," : "), sort_keys=True,
                   ensure_ascii=False) + "\n")
    print("\nwrote %s_Viral_PSSM.json (%d features, %d rep contigs)"
          % (MODULE, len(feats), len(close)))

main()
