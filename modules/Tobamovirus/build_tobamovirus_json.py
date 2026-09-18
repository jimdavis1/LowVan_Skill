#!/usr/bin/env python3
"""Generate the Tobamovirus block of Viral_PSSM.json.

Four features, and the pair that matters is REP126 / REP183.

THE READTHROUGH. Tobamovirus replicase is translated as a 126K protein that
terminates at a leaky amber, and a 183K product made when the ribosome reads
through it. Both are real proteins with different jobs, so both are declared --
exactly as Togaviridae declares P123 alongside NSP1234.

Measured in this dump, not assumed:

  full-length replicase-like sequences   1207
    1100-1199 aa   578      the terminated 126K form
    1200-1399 aa     5      the valley
    1500-1699 aa   607      the 183K readthrough form
  median short 1116, median long 1616    difference exactly 500 aa

  genomes carrying BOTH forms            481
  long form begins with the short form   481 of 481  (100%)
  modal extension                        499 codons

So the two are collinear and one is a prefix of the other. REP183 therefore
carries `internal_stop: 1`, and its alignments must span the amber or the
profile has never seen what follows it. They do: every one of the nine
first-pass REP183 alignments carries an X at 69% of its length, in the same
column, with ~1600 columns spanning it.

BV-BRC already encodes the readthrough position as X in 565 of 613 REP183
sequences, so `fasta-cluster-pssm-2.pl` needs `-x` on both replicase features
or it discards the entire sequence for containing an ambiguous residue -- which
would have left REP183 with 48 sequences instead of 613.

REP126's C-terminus IS the amber, so it gets `downstream_ext: 0`: there is no
codon further on for the annotator to find, and letting it scan would walk into
the readthrough region.
"""
import json, os, collections, sys

MODULE  = "Tobamovirus"
SEGMENT = "Single RNA Segment"

#  Every id model-proposed; PMID stays absent on every feature. Each was checked
#  to resolve to the paper named. Resolving is not reading, and reading is not
#  curating, so all stay flagged.
CLAUDE_PMIDS = {
    "16453503": "Pelham 1984 EMBO J, UAG readthrough during TMV RNA translation: isolation and sequence of two rabbit tRNA-Tyr species",
    "16453524": "1984 EMBO J, the molecular basis for the differential translation of TMV RNA",
    "6964389":  "Goelet et al. 1982 PNAS, Nucleotide sequence of tobacco mosaic virus RNA",
    "17794341": "Deom, Oliver & Beachy 1987 Science, The 30-kilodalton gene product of tobacco mosaic virus potentiates virus movement",
    "2769760":  "Namba, Pattanayek & Stubbs 1989 J Mol Biol, Refined structure of intact tobacco mosaic virus at 2.9 A by X-ray fiber diffraction",
    "28786782": "Adams et al. 2017 J Gen Virol, ICTV Virus Taxonomy Profile: Virgaviridae",
}
READ_PMIDS = set()

FEATURES = {
    "REP126": dict(anno="Methyltransferase and helicase replication protein", symbol="126K", ftype="CDS",
                   bit=600, cov=0.65, up=1, down=0, copy=1,
                   pmid=["6964389", "16453503", "28786782"]),
    "REP183": dict(anno="RNA-dependent RNA polymerase", symbol="183K", ftype="CDS",
                   bit=800, cov=0.65, up=1, down=1, copy=1, internal_stop=1,
                   pmid=["6964389", "16453503", "16453524", "28786782"]),
    "MP":     dict(anno="Movement protein", symbol="MP", ftype="CDS",
                   bit=120, cov=0.65, up=1, down=1, copy=1,
                   pmid=["17794341", "28786782"]),
    "CP":     dict(anno="Nucleocapsid protein", symbol="CP", ftype="CDS",
                   bit=80, cov=0.65, up=1, down=1, copy=1,
                   pmid=["2769760", "28786782"]),
}

def lengths(key):
    L=[]; s=None
    for line in open(os.path.join("collections", MODULE, key + ".fasta")):
        if line.startswith(">"):
            if s is not None: L.append(len(s))
            s=""
        else: s+=line.strip()
    if s is not None: L.append(len(s))
    return sorted(L)

def main():
    close=json.load(open("Rep-Contigs/close_genomes.json"))
    feats={}
    print("%-10s %5s %6s %6s %8s %8s"%("feature","n","min","max","min_len","max_len"))
    for key,d in FEATURES.items():
        L=lengths(key); lo=int(L[0]*0.92); hi=int(L[-1]*1.08)
        print("%-10s %5d %6d %6d %8d %8d"%(key,len(L),L[0],L[-1],lo,hi))
        f=collections.OrderedDict()
        f["anno"]=d["anno"]; f["bit_cutoff"]=d["bit"]
        if "copy" in d: f["copy_num"]=d["copy"]
        f["coverage_cutoff"]=d["cov"]; f["downstream_ext"]=d["down"]
        f["feature_type"]=d["ftype"]; f["gene_symbol"]=d["symbol"]
        if d.get("internal_stop"): f["internal_stop"]=1
        f["kmers"]=1; f["max_len"]=hi; f["min_len"]=lo
        f["segment"]=SEGMENT; f["upstream_ext"]=d["up"]
        ids=[str(x) for x in d["pmid"]]
        prop=[p for p in ids if p in CLAUDE_PMIDS]
        cur=[p for p in ids if p not in CLAUDE_PMIDS]
        if cur:  f["PMID"]=cur
        if prop: f["PMID_claude_generated"]=prop
        unk=[p for p in ids if p not in CLAUDE_PMIDS and p not in READ_PMIDS]
        if unk: sys.exit("PMID %s in neither registry"%unk)
        feats[key]=f
    block={MODULE: collections.OrderedDict([
        ("close_genomes",close),
        ("features",feats),
        ("segments",{SEGMENT:{"max_len":8000,"min_len":5500,"replicon_geometry":"linear"}}),
    ])}
    open("Tobamovirus_Viral_PSSM.json","w").write(
        json.dumps(block,indent=3,separators=(","," : "),sort_keys=True,ensure_ascii=False)+"\n")
    print("\nwrote Tobamovirus_Viral_PSSM.json (%d features, %d rep contigs)"%(len(feats),len(close)))
main()
