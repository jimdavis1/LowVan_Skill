#!/usr/bin/env python3
"""Merhavirus module JSON.

Split out of Alpharhabdovirinae so that CTRV's spliced L can be called.

Culex tritaeniorhynchus rhabdovirus carries a 76-nt GU-AG intron in its L
coding region (Kuwata et al. 2011, J Virol, PMID 21507977; antigenome
8648-8723). It replicates in the nucleus, which is what makes splicing
possible at all -- unique among Mononegavirales outside Bornaviridae. Without
the splice, the ordinary L profile emits L truncated at the in-frame stop at
8681: not a missing call but a wrong one.

`special: splice` is called by get_splice_variant_features.pl from
hand-curated nucleotide references in Splice-Variants/<Module>/<FEAT>.fasta.
That program is keyed on the module, so the feature could not work while it
lived in a module whose reference set covered all of Alpharhabdovirinae. It
needs a module of its own, which is what this is -- the same shape as
Orthopneumovirus_muris split from Orthopneumovirus.

Eight features, seven of them inherited from the parent unchanged. The
alignments are the parent's own clusters that contain Merhavirus sequences,
16 of them across G, G_MAT, G_SP, L, M, N and P, regenerated under this
module's name. The parent keeps its copies: Merhavirus sequences were part of
what those profiles were built from, and removing them would change the
parent's calls on everything else.

Rep contigs: 8 references cover 71 of 71 genomes with a contig >= 9 kb.
Merhavirus.2.dna is LC514403, CTRV, covering all 13 CTRV genomes.
Alpharhabdovirinae.8.dna (NC_025384.1, also CTRV) comes out of the parent's
set, or the genus routes to the module that cannot call its L.

CITATIONS. Every id here is inherited from Alpharhabdovirinae, where it sits
in `PMID`. Whether a curator read them cannot be established from the
repository -- the module predates it -- so they are carried as
`PMID_claude_generated` and listed in the build notes for promotion. That is
the conservative reading of the rule: PMID is for citations a human has
confirmed.
"""
import json, os, collections, sys

MODULE  = "Merhavirus"
SEGMENT = "Single RNA Segment"

CLAUDE_PMIDS = {
    "35723908": "Walker et al. 2022 J Gen Virol, ICTV Virus Taxonomy Profile: Rhabdoviridae (inherited from Alpharhabdovirinae)",
    "6897030":  "Anilionis, Wunner & Curtis 1981 Nature, Structure of the glycoprotein gene in rabies virus (inherited)",
    "21507977": "Kuwata et al. 2011 J Virol, RNA splicing in a new rhabdovirus from Culex mosquitoes",
}
READ_PMIDS = set()

FEATURES = collections.OrderedDict([
 ("N",        dict(anno="Nucleocapsid protein", symbol="N", ftype="CDS",
                   bit=100, cov=0.65, up=1, down=1, copy=1, kmers=1,
                   lo=312, hi=550, pmid=["35723908"])),
 ("P",        dict(anno="Phosphoprotein", symbol="P", ftype="CDS",
                   bit=80, cov=0.65, up=1, down=1, copy=1, kmers=1,
                   lo=174, hi=363, pmid=["35723908"])),
 ("M",        dict(anno="Matrix protein", symbol="M", ftype="CDS",
                   bit=60, cov=0.65, up=1, down=1, copy=1, kmers=1,
                   lo=146, hi=307, pmid=["35723908"])),
 ("G",        dict(anno="Envelope glycoprotein precursor", symbol="G", ftype="CDS",
                   bit=100, cov=0.65, up=1, down=1, copy=1, kmers=1,
                   lo=398, hi=706, pmid=["35723908"])),
 ("G_SP",     dict(anno="Signal peptide of G", symbol="SP", ftype="mat_peptide",
                   bit=30, cov=0.65, up=1, down=0, kmers=0,
                   lo=14, hi=33, pmid=["6897030"])),
 ("G_MAT",    dict(anno="Mature envelope glycoprotein", symbol="G_mature", ftype="mat_peptide",
                   bit=100, cov=0.65, up=0, down=1, kmers=0,
                   lo=423, hi=682, pmid=["6897030"])),
 ("L",        dict(anno="RNA-dependent RNA polymerase", symbol="L", ftype="CDS",
                   bit=100, cov=0.65, up=1, down=1, copy=1, kmers=1,
                   lo=1829, hi=2449, pmid=["35723908"])),
 #  No copy_num, deliberately. viral_genome_quality.pl keys %essential by the
 #  annotation string, and L already claims "RNA-dependent RNA polymerase";
 #  giving this one a copy_num too would make the two features overwrite each
 #  other's length windows. %anno_count still counts a spliced L against L's
 #  copy_num of 1, which is exactly the behaviour to test: the splice must
 #  REPLACE the truncated call, not add to it.
 ("L_SPLICED", dict(anno="RNA-dependent RNA polymerase", symbol="L", ftype="CDS",
                   lo=1910, hi=2336, special="splice", pmid=["21507977"])),
])

def main():
    close = json.load(open("Rep-Contigs/close_genomes.json"))
    feats = {}
    for key, d in FEATURES.items():
        f = collections.OrderedDict()
        f["anno"] = d["anno"]
        if "bit" in d: f["bit_cutoff"] = d["bit"]
        if "copy" in d: f["copy_num"] = d["copy"]
        if "cov" in d: f["coverage_cutoff"] = d["cov"]
        if "down" in d: f["downstream_ext"] = d["down"]
        f["feature_type"] = d["ftype"]; f["gene_symbol"] = d["symbol"]
        if "kmers" in d: f["kmers"] = d["kmers"]
        f["max_len"] = d["hi"]; f["min_len"] = d["lo"]
        f["segment"] = SEGMENT
        if d.get("special"): f["special"] = d["special"]
        if "up" in d: f["upstream_ext"] = d["up"]
        ids = [str(x) for x in d["pmid"]]
        unk = [p for p in ids if p not in CLAUDE_PMIDS and p not in READ_PMIDS]
        if unk: sys.exit("PMID %s in neither registry" % unk)
        prop = [p for p in ids if p in CLAUDE_PMIDS]
        cur  = [p for p in ids if p in READ_PMIDS]
        if cur:  f["PMID"] = cur
        if prop: f["PMID_claude_generated"] = prop
        feats[key] = f

    essential = [v["anno"] for v in feats.values() if v.get("copy_num")]
    if len(essential) != len(set(essential)):
        sys.exit("two essential features share an annotation string")

    block = {MODULE: collections.OrderedDict([
        ("close_genomes", close),
        ("features", feats),
        ("segments", {SEGMENT: {"max_len": 13000, "min_len": 9000,
                                "replicon_geometry": "linear"}}),
    ])}
    open("%s_Viral_PSSM.json" % MODULE, "w").write(
        json.dumps(block, indent=3, separators=(","," : "), sort_keys=True,
                   ensure_ascii=False) + "\n")
    print("wrote %s_Viral_PSSM.json (%d features, %d rep contigs)"
          % (MODULE, len(feats), len(close)))
    for k, v in feats.items():
        print("   %-12s %-10s %-34s %s" % (k, v["gene_symbol"], v["anno"][:34],
                                           v.get("special") or ""))
main()
