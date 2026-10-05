#!/usr/bin/env python3
"""Write the module JSON for the four grouped Picornaviridae modules.

Windows are 0.8x/1.2x the collection median and bit cutoffs max(20, 0.45x the
median length) -- the same rule every other picornavirus feature was built
with. Annotation strings reuse the controlled vocabulary; the only judgement
here is which 2A string a group gets, since the slot holds a StopGo peptide in
some genera and a full protein in others, and that is decided by measuring the
collection rather than by the slot name.
"""
import collections, json, os

S = os.path.expanduser("~/lowvan_scratch/picorna")
GROUPS = json.load(open(S + "/rules/groups.json"))

ANNO = {
    "POLY": ("Genome polyprotein", "POLY"),
    "Ldr":  ("Leader protein", "Ldr"),
    #  ERBV L has C-terminal processing activity -- Hinton 2002 J Gen Virol
    #  83:3111 (PMID 12466488). Model-proposed citation.
    "Ldr_PRO":  ("Leader protease (Lpro)", "Ldr"),
    #  Sicinivirus carries a 462-residue leader against 64-174 for the rest of
    #  its group. No function is documented for it, so the vocabulary's
    #  uncharacterized form is used and the tag stays in the key and symbol.
    "Ldr_LONG": ("Uncharacterized lineage-specific protein", "Ldr-long"),
    "VP4": ("Capsid protein VP4 (1A)", "VP4"),
    "VP2": ("Capsid protein VP2 (1B)", "VP2"),
    "VP0": ("Capsid protein VP0 (1AB)", "VP0"),
    "VP3": ("Capsid protein VP3 (1C)", "VP3"),
    "VP1": ("Capsid protein VP1 (1D)", "VP1"),
    "2A":      (None, "2A"),          # decided from the median, below
    "2A_LONG": ("Protein 2A", "2A-long"),
    "2B":  ("Protein 2B, viroporin", "2B"),
    "2C":  ("Protein 2C, ATPase/helicase", "2C"),
    "3A":  ("Protein 3A", "3A"),
    "VPG": ("Genome-linked protein VPg (3B)", "VPg"),
    "3C":  ("Picornain 3C, protease (3C)", "3C"),
    "3D":  ("RNA-directed RNA polymerase 3D", "3D"),
}
PMID = {"Ldr_PRO": {"function": ["12466488"]}}

ORDER = ["POLY", "Ldr", "Ldr_PRO", "Ldr_LONG", "VP0", "VP4", "VP2", "VP3", "VP1",
         "2A", "2A_LONG", "2B", "2C", "3A", "VPG", "3C", "3D"]


def lens(p):
    L, cur = [], 0
    for line in open(p):
        if line.startswith(">"):
            if cur: L.append(cur)
            cur = 0
        else: cur += len(line.strip())
    if cur: L.append(cur)
    return sorted(L)


for G in GROUPS:
    C = f"{S}/collections/{G}"
    keys = [k for k in ORDER if os.path.exists(f"{C}/{k}.fasta") and lens(f"{C}/{k}.fasta")]
    feats = {}
    for k in keys:
        L = lens(f"{C}/{k}.fasta"); med = L[len(L) // 2]
        anno, sym = ANNO[k]
        if k == "2A":
            #  measured, not assumed: a StopGo peptide is tens of residues, a
            #  2A protein is low hundreds
            anno = "Protein 2A, StopGo peptide" if med < 60 else "Protein 2A"
        f = {
            "anno": anno, "gene_symbol": sym,
            "feature_type": "CDS" if k == "POLY" else "mat_peptide",
            "type": "CDS" if k == "POLY" else "mat_peptide",
            "bit_cutoff": max(20, int(round(0.45 * med))),
            "coverage_cutoff": 0.65,
            "min_len": max(1, int(round(0.8 * med))), "max_len": int(round(1.2 * med)),
            "upstream_ext": 1 if k == "POLY" else 0, "downstream_ext": 0,
            "kmers": 1 if k == "POLY" else 0,
            "segment": "Single RNA Segment",
        }
        if k == "POLY": f["copy_num"] = 1
        if k in PMID: f["PMID_claude_generated"] = PMID[k]
        feats[k] = f
    #  two features in one module must never share an annotation string: the
    #  quality tool keys %essential by the string, so the second read silently
    #  overwrites the first's length window (the Tobamovirus defect)
    seen = collections.Counter(f["anno"] for f in feats.values())
    dup = [a for a, n in seen.items() if n > 1]
    if dup: raise SystemExit("%s: duplicate annotation strings %s" % (G, dup))
    block = {G: {"features": feats,
                 "segments": {"Single RNA Segment": {"min_len": 5000, "max_len": 12000,
                                                     "replicon_geometry": "linear"}}}}
    def norm(x):
        if isinstance(x, dict): return {k: norm(v) for k, v in x.items()}
        if isinstance(x, list): return [norm(v) for v in x]
        if isinstance(x, float) and x == int(x): return int(x)
        return x
    p = f"{S}/{G}_Viral_PSSM.json"
    open(p, "w").write(json.dumps(norm(block), indent=3, separators=(",", " : "),
                                  sort_keys=True, ensure_ascii=False) + "\n")
    print("  %-16s %2d features: %s" % (G, len(feats), " ".join(keys)))
