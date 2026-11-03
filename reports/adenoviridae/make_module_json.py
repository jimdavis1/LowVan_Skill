#!/usr/bin/env python3
"""Write <Module>_Viral_PSSM.json for the four Adenoviridae modules.

ANNOTATION STRINGS. The governing rule from the manuscript: the annotation
says what the protein DOES; the commonly used symbolic name goes in
`gene_symbol`, not in the string. So "Major capsid protein hexon" with symbol
"hexon", never "hexon protein" as the annotation.

Adenoviridae contributes no strings to the controlled vocabulary yet -- the
shipped TSV covers 35 taxa, all of them RNA viruses -- so these are coined,
following the house style exactly: "Mature ..." for a cleavage product,
"Uncharacterized lineage-specific protein" for a discovered accessory cluster
that the literature does not name.

WHERE A NUMBER IS KEPT ON PURPOSE. Aviadenovirus ORF9/10/11/14/17/20,
Siadenovirus ORF7/ORF8, Barthadenovirus RH0/RH3/RH5/LH1/E4.1 all carry a
position label rather than a function, because no paper assigns them one. A
functional string the literature does not support is worse than a number.

LENGTH WINDOWS are read from each collection's own distribution, padded, and
floored -- never from a reference genome, which is a database record and can
be wrong about a boundary (see dbcoords.json).
"""
import collections, json, os, sys

W = os.path.dirname(os.path.abspath(__file__)) or "."

#  ---- the 17 genes every adenovirus carries -------------------------------
CORE = {
 "HEXON":    ("Major capsid protein hexon", "hexon"),
 "PENTON":   ("Penton base protein", "penton"),
 "FIBER":    ("Receptor binding fiber protein", "fiber"),
 "PROTEASE": ("Adenain cysteine protease", "protease"),
 "POL":      ("DNA polymerase", "pol"),
 "PTP":      ("Terminal protein precursor", "pTP"),
 "DBP":      ("Single-stranded DNA-binding protein", "DBP"),
 "IVA2":     ("Genome packaging ATPase IVa2", "IVa2"),
 "K52":      ("Genome packaging protein 52/55K", "52K"),
 "K100":     ("Hexon assembly chaperone and host shutoff protein 100K", "100K"),
 "K33":      ("Alternative splicing factor 33K", "33K"),
 "K22":      ("Late transcription regulator 22K", "22K"),
 "PIIIA":    ("Capsid vertex protein IIIa precursor", "IIIa"),
 "PVI":      ("Capsid cement protein VI precursor", "pVI"),
 "PVII":     ("Core protein VII precursor", "pVII"),
 "PVIII":    ("Capsid cement protein VIII precursor", "pVIII"),
 "PX":       ("Core protein X precursor", "pX"),
 "UXP":      ("U exon protein", "UXP"),
}

MAST = {
 "PV":        ("Core protein V", "pV"),
 "PIX":       ("Capsid cement protein IX", "pIX"),
 "FIBER1":    ("Receptor binding long fiber protein", "fiber-1"),
 "FIBER2":    ("Receptor binding short fiber protein", "fiber-2"),
 "E1A":       ("Transcriptional activator and retinoblastoma protein binding E1A", "E1A"),
 "E1A_243R":  ("Transcriptional activator E1A 243R alternative form", "E1A"),
 "E1B19K":    ("Apoptosis inhibitor and Bcl-2 homolog E1B 19K", "E1B-19K"),
 "E1B55K":    ("p53 degradation and mRNA export factor E1B 55K", "E1B-55K"),
 "E3_12_5K":  ("E3 12.5K protein", "E3-12.5K"),
 "E3_GP19K":  ("Major histocompatibility complex class I retention glycoprotein E3 gp19K", "E3-gp19K"),
 "E3_6_7K":   ("Calcium homeostasis membrane protein E3 6.7K", "E3-6.7K"),
 "E3_CR1A":   ("E3 conserved region 1 alpha protein", "E3-CR1a"),
 "E3_CR1B":   ("Immunomodulatory glycoprotein E3 conserved region 1 beta", "E3-CR1b"),
 "E3_CR1G":   ("E3 conserved region 1 gamma protein", "E3-CR1g"),
 "E3_CR1D":   ("E3 conserved region 1 delta protein", "E3-CR1d"),
 "E3_RIDA":   ("Receptor internalization and degradation complex subunit alpha", "RID-alpha"),
 "E3_RIDB":   ("Receptor internalization and degradation complex subunit beta", "RID-beta"),
 "E3_14_7K":  ("Tumor necrosis factor antagonist E3 14.7K", "E3-14.7K"),
 "E3_ADP":    ("Adenovirus death protein", "ADP"),
 "E4ORF1":    ("PDZ-binding transforming protein E4orf1", "E4orf1"),
 "E4ORF2":    ("E4orf2 protein", "E4orf2"),
 "E4ORF3":    ("Promyelocytic leukemia body disrupting protein E4orf3", "E4orf3"),
 "E4ORF4":    ("Protein phosphatase 2A regulatory subunit E4orf4", "E4orf4"),
 "E4ORF5":    ("E4orf5 protein", "E4orf5"),
 "E4ORF6":    ("p53 degradation and mRNA export factor E4orf6", "E4orf6"),
 "E4ORF67":   ("E2F dimerization partner E4orf6/7", "E4orf6/7"),
 "I_LEADER":  ("i-leader protein", "i-leader"),
 "L1_13_6K":  ("L1 13.6K protein", "L1-13.6K"),
 "AGNO":      ("Agnoprotein", "agno"),
}

AVI = {
 "FIBER1":    ("Receptor binding fiber-1 protein", "fiber-1"),
 "FIBER2":    ("Receptor binding fiber-2 protein", "fiber-2"),
 "GAM1":      ("Anti-apoptotic and heat shock response inducer GAM-1", "GAM-1"),
 "DUTPASE":   ("Deoxyuridine triphosphate nucleotidohydrolase", "dUTPase"),
 "REPLIKE":   ("Parvovirus Rep-like protein with inactive helicase domain ORF2", "ORF2"),
 "LIPASE":    ("Lipase-like protein ORF19", "ORF19"),
 "ART":       ("Mono-ADP-ribosyltransferase-like protein ORF16", "ORF16"),
 "ORF22":     ("Retinoblastoma protein binding and E2F activating protein ORF22", "ORF22"),
 "ORF12":     ("Parvovirus Rep-like protein with inactive helicase domain ORF12", "ORF12"),
 "ORF13":     ("Parvovirus Rep-like protein with inactive helicase domain ORF13", "ORF13"),
 "ORF9":      ("ORF9 protein", "ORF9"),
 "ORF10":     ("ORF10 protein", "ORF10"),
 "ORF11":     ("ORF11 protein", "ORF11"),
 "ORF14":     ("ORF14 protein", "ORF14"),
 "ORF17":     ("ORF17 protein", "ORF17"),
 "ORF20":     ("ORF20 protein", "ORF20"),
}

SIA = {
 "SIALIDASE": ("Sialidase-like virion protein", "ORF1"),
 "HYD":       ("Hydrophobic protein", "hyd"),
 "E3":        ("E3 protein", "E3"),
 "ORF7":      ("ORF7 protein", "ORF7"),
 "ORF8":      ("ORF8 protein", "ORF8"),
}

BARTHA = {
 "P32K":      ("Structural protein p32K precursor", "p32K"),
 "P32K_MAT":  ("Mature structural protein p32K", "p32K"),
 "LH1":       ("LH1 protein", "LH1"),
 "LH2":       ("Structural protein LH2", "LH2"),
 "LH3":       ("Capsid cementing protein LH3", "LH3"),
 "RH0":       ("RH0 protein", "RH0"),
 "RH1":       ("F-box protein RH1", "RH1"),
 "RH2":       ("F-box protein RH2", "RH2"),
 "RH3":       ("RH3 protein", "RH3"),
 "RH4":       ("F-box protein RH4", "RH4"),
 "RH5":       ("RH5 protein", "RH5"),
 "RH6":       ("F-box protein RH6", "RH6"),
 "E4_1":      ("E4.1 protein", "E4.1"),
 "E4_2":      ("E4.2 protein", "E4.2"),
 "E4_3":      ("E4.3 protein", "E4.3"),
}

ANNO = {"Mastadenovirus":  dict(CORE, **MAST),
        "Aviadenovirus":   dict(CORE, **AVI),
        "Siadenovirus":    dict(CORE, **SIA),
        "Barthadenovirus": dict(CORE, **BARTHA)}

#  Whole-genome length window, in NUCLEOTIDES, for the segments block.
#  Measured from the dump's own genome lengths, padded.
SEGMENT = {"Mastadenovirus":  (24000, 40000),
           "Aviadenovirus":   (38000, 48000),
           "Siadenovirus":    (22000, 29000),
           "Barthadenovirus": (26000, 36000)}

#  Features that are mature peptides, not CDS. They sit INSIDE their parent
#  and must carry upstream_ext 0 -- with the default of 1 the annotator scans
#  upstream from the cleaved N-terminus, finds the precursor's start codon and
#  emits the mature peptide on the precursor's own coordinates. See
#  references/json-schema.md; it shipped that way on rabies G and on six
#  norovirus products at once.
MATPEP = {"P32K_MAT"}


def lengths(path):
    L, cur = [], 0
    for l in open(path):
        if l.startswith(">"):
            if cur:
                L.append(cur)
            cur = 0
        else:
            cur += len(l.strip())
    if cur:
        L.append(cur)
    return sorted(L)


def window(L):
    """Permissive window from the collection's own distribution.

    Trim the extreme 2% at each end so one mis-binned deposit cannot set the
    bound, then pad 15%. Permissive on purpose: a tight window silently
    DELETES any species whose protein is genuinely a different size, and that
    failure is invisible -- it shows up as a missing call, not an error.
    """
    if not L:
        return 20, 5000
    n = len(L)
    lo = L[max(0, int(0.02 * n))] if n >= 20 else L[0]
    hi = L[min(n - 1, int(0.98 * n))] if n >= 20 else L[-1]
    return max(20, int(lo * 0.85)), int(hi * 1.15) + 5


def main():
    for mod, anno in ANNO.items():
        cd = os.path.join(W, mod, "collections")
        feats, unc, missing = {}, 0, []
        for f in sorted(os.listdir(cd)):
            if not f.endswith(".fasta"):
                continue
            k = f[:-6]
            L = lengths(os.path.join(cd, f))
            if not L:
                continue
            if k in anno:
                a, sym = anno[k]
            elif k.startswith("UNC"):
                unc += 1
                a, sym = "Uncharacterized lineage-specific protein %s" % k[3:], k
            else:
                missing.append(k)
                continue
            lo, hi = window(L)
            d = {"anno": a, "gene_symbol": sym,
                 "feature_type": "mat_peptide" if k in MATPEP else "CDS",
                 "type": "mat_peptide" if k in MATPEP else "CDS",
                 "min_len": lo, "max_len": hi,
                 "coverage_cutoff": 0.60, "kmers": 1,
                 "upstream_ext": 0 if k in MATPEP else 1,
                 "downstream_ext": 1,
                 "segment": "Single DNA Segment",
                 "n_seqs": len(L), "median_aa": L[len(L) // 2]}
            feats[k] = d
        lo, hi = SEGMENT[mod]
        out = {mod: {"segments": {"Single DNA Segment":
                                  {"min_len": lo, "max_len": hi}},
                     "features": feats}}
        p = os.path.join(W, "%s_Viral_PSSM.json" % mod)
        json.dump(out, open(p, "w"), indent=3, sort_keys=True)
        print("%-18s %3d features (%d uncharacterized)%s"
              % (mod, len(feats), unc,
                 ("  NO ANNOTATION: " + ", ".join(missing)) if missing else ""))


main()
