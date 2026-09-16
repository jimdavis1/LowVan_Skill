#!/usr/bin/env python3
"""Recover the two multi-product intermediates from their flanking calls.

    PE2  (p62)  = E3.start   -> E2.end     covalent precursor, furin-cleaved late
    P123        = nsP1.start -> nsP3.end   the opal-terminated nonstructural product

Neither is in the BV-BRC export in useful numbers -- 8 sequences for pE2 (one of
them full length) and 2 for P123 -- but both flanking features are called
reliably, and the interval between them IS the protein. Same construction that
recovered full-length nsP3.

P123 contains the opal as an internal stop in readthrough-competent genomes; it
is written X for the same reasons as nsP3 (mafft turns '*' into a gap, psiblast
refuses an MSA containing one).
"""
import collections, glob, os

ANN="coverage_eval_f/ann"; EX="coverage_eval/ex"
#  key: (flanking features, min aa, max aa, must start at an initiating Met)
#  P123 begins at the polyprotein's own start codon, so a recovered P123 that
#  does not start with M is N-terminally truncated -- its left flank (nsP1) was
#  itself called short. Admitting those built a 23-member cluster whose
#  alignment began 190 residues in at `WVGFDTTPFMY...`, and that profile then
#  won on Chikungunya and placed P123 122 residues late. pE2 begins at a
#  cleavage site, not a start codon, so the same test must NOT be applied to it.
SPECS={"PE2": (("E3","E2"), 380, 560, False),
       "P123":(("nsP1","nsP3"), 1650, 2100, True)}
T={}; _b="TCAG"; _a="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"; i=0
for x in _b:
    for y in _b:
        for z in _b: T[x+y+z]=_a[i]; i+=1
tr=lambda s:"".join(T.get(s[j:j+3],"X") for j in range(0,len(s)-2,3))
rc=lambda s:s[::-1].translate(str.maketrans("ACGTN","TGCAN"))
meta={}
for l in open("Toga.genome_meta.tsv"):
    p=l.rstrip("\n").split("\t")
    if len(p)>=7: meta[p[0]]=dict(sp=(p[6].strip() or "?"),name=p[2])
os.makedirs("Extracted",exist_ok=True)
for key,((fa,fb),LO,HI,need_met) in SPECS.items():
    stats=collections.Counter(); out=[]; bysp=collections.Counter()
    for t in sorted(glob.glob(os.path.join(ANN,"*.feature.tbl"))):
        if not os.path.getsize(t): continue
        g=os.path.basename(t).replace(".feature.tbl","")
        f=collections.defaultdict(list)
        for line in open(t,errors="replace"):
            c=line.rstrip("\n").split("\t")
            if len(c)<10: continue
            try: lo,hi=int(c[7]),int(c[8])
            except ValueError: continue
            f[c[6]].append((lo,hi,c[9]))
        if len(f.get(fa,[]))!=1 or len(f.get(fb,[]))!=1:
            stats["flanks not called exactly once"]+=1; continue
        (a1,a2,sa),(b1,b2,sb)=f[fa][0],f[fb][0]
        if sa!=sb: stats["opposite strands"]+=1; continue
        p=os.path.join(EX,g+".fna")
        if not os.path.exists(p): stats["no sequence"]+=1; continue
        s="".join(x.strip() for x in open(p) if not x.startswith(">")).upper()
        if sa=="+": lo,hi=min(a1,a2),max(b1,b2); nt=s[lo-1:hi]
        else:       lo,hi=min(b1,b2),max(a1,a2); nt=rc(s[lo-1:hi])
        if hi<=lo: stats["inverted"]+=1; continue
        if len(nt)%3: stats["not a multiple of 3"]+=1; continue
        aa=tr(nt)
        if not (LO<=len(aa)<=HI): stats[f"length outside {LO}-{HI}"]+=1; continue
        ns=aa.count("*")
        if ns>1: stats["more than one stop -- broken"]+=1; continue
        if aa.count("X")/len(aa)>0.01: stats["too ambiguous"]+=1; continue
        if need_met and not aa.startswith("M"):
            stats["N-terminally truncated (no initiating Met)"]+=1; continue
        if ns==1: aa=aa.replace("*","X"); stats["recovered, carries the opal (-> X)"]+=1
        else: stats["recovered"]+=1
        bysp[meta.get(g,{}).get("sp","?")]+=1
        out.append((g,meta.get(g,{}).get("name","?"),meta.get(g,{}).get("sp","?"),aa))
    with open(f"Extracted/Togaviridae.{key}.faa","w") as fh:
        fh.write(f"; {key} recovered from the {fa}-{fb} interval. See extract_combos.py.\n")
        for g,nm,sp,aa in out: fh.write(f">extracted|{g} {nm} | {sp}\n{aa}\n")
    L=sorted(len(a[3]) for a in out)
    print(f"{key}: {len(out)} sequences, {len(bysp)} species, lengths {L[0]}-{L[-1]} (median {L[len(L)//2]})")
    for k,v in stats.most_common(4): print(f"     {k:44s} {v}")
