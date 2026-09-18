#!/usr/bin/env python3
"""Triage Betaflexiviridae_MP into feature collections.

The movement-protein half of Betaflexiviridae: Trichovirus, Vitivirus,
Citrivirus, Prunevirus, Tepovirus, Chordovirus, Divavirus, Wamavirus. These
carry a single MP where the TGB half of the family carries a triple gene block,
which is why they are a separate module -- and, as it turned out, why they
route better apart than together (91.6% of genomes on 25 references, against
75.3% for the family lumped).

Capillovirus is deliberately NOT here. Its coat protein sits inside the ORF1
polyprotein (residues 1909-2100 of 2105, measured by BLASTp against Trichovirus
and Tepovirus CPs), so it has no separate CP CDS and cannot share this module's
CP feature.

Four features, separated by length as much as by string, because "replicase",
"RNA-dependent RNA polymerase" and "polyprotein" all name the same ~1800 aa
protein while "hypothetical protein" covers several different small ones:

    REP   ~1700-1990 aa   the replicase polyprotein
    MP     ~290-460 aa    movement protein
    CP     ~190-250 aa    coat protein
    NABP    ~94-128 aa    nucleic-acid binding protein, largely Vitivirus
"""
import re, os, sys, collections, argparse

MODULE = "Betaflexiviridae_MP"
DUMP   = "Betaflexiviridae_MP"

def rx(p): return re.compile(p, re.I)

#  VITI_ORF2 must be tested BEFORE the generic rules, because its records are
#  mostly "hypothetical protein" and the few that are named use a size label.
#  It is a real homology group, measured: 426 Vitivirus sequences of 140-230 aa,
#  67 of 70 sampled hit another member at a median 97% identity, and NONE of
#  them hits CP or NABP at all. No paper establishes a function for it, so it
#  ships as uncharacterized rather than with an invented one.
RULES = [
    ("VITI_ORF2", "Vitivirus", rx(r"^\d{2} ?kda.?protein$|19kda|^orf ?2$|^p19$|^p20$")),
    ("NABP", None, rx(r"nucleic.acid.bind|rna.bind|\bnabp\b|nucleic.acid.biding")),
    ("MP",   None, rx(r"movement|\bmp\b|cell.to.cell")),
    ("CP",   None, rx(r"\bcoat\b|capsid|\bcp\b")),
    ("REP",  None, rx(r"replicas|replicat|polymeras|\brdrp\b|polyprotein|helicase|methyltransferase")),
]
FORBIDDEN = [rx(r"^hypothetical protein"), rx(r"^unknown"), rx(r"^orf ?\d+$"),
             rx(r"^unnamed protein product$"), rx(r"^putative protein$"),
             rx(r"^p\d+$"), rx(r"^\d+ ?k(da)? protein$")]

EXPECTED = {
    "REP":  (1600, 2150),
    "MP":   ( 250,  480),
    "CP":   ( 150,  280),
    "NABP": (  80,  145),
    "VITI_ORF2": (140, 230),
}
MISBINNED = {}
#  Sequences the text rules cannot reach -- 395 of the 426 VITI_ORF2 records say
#  only "hypothetical protein" -- are recovered by the homology rescue, which is
#  the documented route for a protein whose annotation is not a name.
NOT_MODELLED = [
    ("CP (Capillovirus)", "Capillovirus", 0,
     "Capillovirus is a separate module: its coat protein is inside the ORF1 "
     "polyprotein (residues 1909-2100 of 2105) rather than a CDS of its own."),
]

def normalise(a):
    n = re.sub(r"[-_]", " ", a); return re.sub(r"\s+", " ", n).strip()

def classify(ann, genus):
    for cand in (ann, normalise(ann)):
        for key, g, prx in RULES:
            if g and g != genus: continue
            if prx.search(cand): return key
    return None


def genus_specific(ann, genus):
    """A genus-scoped rule beats the FORBIDDEN list.

    FORBIDDEN blocks size labels like "19 kDa protein" because in general a
    size is not a protein name. But for Vitivirus specifically it IS the name
    this protein is published under, and the homology test says the sequences
    carrying it are one group. So a rule that names both the genus and the
    string is evidence, not a guess, and is allowed through.
    """
    for cand in (ann, normalise(ann)):
        for key, g, prx in RULES:
            if not g or g != genus: continue
            if prx.search(cand): return key
    return None

def forbidden(a):
    return any(p.search(a) or p.search(normalise(a)) for p in FORBIDDEN)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--workdir",default="."); a=ap.parse_args()
    W=a.workdir
    g2g={}
    for l in open(os.path.join(W,DUMP+".id_name_gs")):
        p=l.rstrip("\n").split("\t"); g2g[p[0]]=p[3] or "?"
    seq=dict(l.rstrip("\n").split("\t") for l in open(os.path.join(W,DUMP+".uniq.seq")))
    md5=[l.strip() for l in open(os.path.join(W,DUMP+".uniq.md5"))]
    ann=[l.rstrip("\n").split("\t") for l in open(os.path.join(W,DUMP+".uniq.id_ann"))]
    out=collections.defaultdict(list); outl=collections.defaultdict(list)
    unass=collections.Counter(); syn=collections.Counter(); hit=collections.Counter()
    for m,(fid,aname) in zip(md5,ann):
        gid=fid.split("|",1)[1].rsplit(".",2)[0]; genus=g2g.get(gid,"?"); s=seq.get(m,"")
        if not s: continue
        key=genus_specific(aname,genus) or (None if forbidden(aname) else classify(aname,genus))
        if fid in MISBINNED: key=MISBINNED[fid][0]
        syn[(aname,genus,key or "UNASSIGNED")]+=1
        if key is None: unass[(aname,genus)]+=1; continue
        lo,hi=EXPECTED[key]; rec=">%s|%s|%s\n%s\n"%(fid,genus,m,s)
        if lo<=len(s)<=hi: out[key].append(rec); hit[key]+=1
        else: outl[key].append(rec)
    cdir=os.path.join(W,"collections",MODULE); os.makedirs(cdir,exist_ok=True)
    for key in EXPECTED:
        open(os.path.join(cdir,key+".fasta"),"w").writelines(out.get(key,[]))
        if outl.get(key): open(os.path.join(cdir,key+".outliers.fasta"),"w").writelines(outl[key])
    cd=os.path.join(W,"collections")
    with open(os.path.join(cd,"synonyms.tsv"),"w") as f:
        f.write("count\tannotation\tgenus\tbins_to\n")
        for (s_,g,k),c in sorted(syn.items(),key=lambda x:-x[1]): f.write("%d\t%s\t%s\t%s\n"%(c,s_,g,k))
    with open(os.path.join(cd,"UNASSIGNED_TRACKING.tsv"),"w") as f:
        f.write("count\tannotation\tgenus\n")
        for (s_,g),c in sorted(unass.items(),key=lambda x:-x[1]): f.write("%d\t%s\t%s\n"%(c,s_,g))
    with open(os.path.join(cd,"NOT_MODELLED.tsv"),"w") as f:
        f.write("feature\ttaxon\tn_sequences\treason\n")
        for r in NOT_MODELLED: f.write("%s\t%s\t%d\t%s\n"%r)
    print("%-8s %8s %10s   %s"%("feature","in-range","outliers","window"))
    for key in EXPECTED:
        print("%-8s %8d %10d   %d-%d aa"%(key,hit[key],len(outl.get(key,[])),*EXPECTED[key]))
    print("\nunassigned: %d sequences (%d strings)   total: %d"%(sum(unass.values()),len(unass),len(md5)))
main()
