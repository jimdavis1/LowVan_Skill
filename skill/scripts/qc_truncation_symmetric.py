#!/usr/bin/env python3
"""Symmetric containment check over every surviving alignment of a feature.

The pipeline's own truncation QC asks "is cluster A's consensus a substring of
cluster B's" using only the row where A is the BLAST *query*. BLAST routinely
trims a few residues off a query-anchored alignment, so a consensus that is
fully contained can score just under `-qc_qcov 0.95` and survive.

Hepeviridae ORF3 cluster 8 is the case: 75 aa, 100% identical to residues 25-99
of a 122-aa ORF3, and it passed. The reverse row (cluster 1 as query) aligned
all 75 of its residues. This takes the better of the two directions.

Leftover profiles are included, which the pipeline's QC never sees at all --
they are built after it runs.
"""
import os, sys, subprocess, tempfile, argparse, collections

def read_fa(p):
    out=[]; h=None; s=[]
    for l in open(p):
        if l.startswith(">"):
            if h is not None: out.append((h,"".join(s)))
            h=l[1:].strip(); s=[]
        else: s.append(l.strip())
    if h is not None: out.append((h,"".join(s)))
    return out

def consensus(recs):
    if not recs: return ""
    w=len(recs[0][1]); out=[]
    for i in range(w):
        c=collections.Counter(r[1][i] for r in recs if i<len(r[1]))
        c.pop("-",None)
        if c: out.append(c.most_common(1)[0][0])
    return "".join(out)

ap=argparse.ArgumentParser()
ap.add_argument("--workdir",default=".")
ap.add_argument("--module",required=True)
ap.add_argument("--key",required=True)
ap.add_argument("--pident",type=float,default=75.0)
ap.add_argument("--qcov",type=float,default=0.95)
ap.add_argument("--lenratio",type=float,default=0.92)
ap.add_argument("--write",action="store_true",help="retire the contained alignments and their PSSMs")
a=ap.parse_args()

base=os.path.join(a.workdir,"Alignments",a.module,a.key)
# Every directory that can hold a live alignment. reclustered_alis/ holds the
# products of the N-terminal split and leftover profiles land in corrected_alis;
# scanning only corrected_alis silently exempts the N-term splits from the QC.
SRC_DIRS=["corrected_alis","reclustered_alis"]
alis={}; home={}
for d in SRC_DIRS:
    p=os.path.join(base,d)
    if not os.path.isdir(p): continue
    for fn in sorted(os.listdir(p)):
        if fn.endswith(".fa"):
            alis[fn[:-3]]=read_fa(os.path.join(p,fn)); home[fn[:-3]]=p
cdir=os.path.join(base,"corrected_alis")
alis={k:v for k,v in alis.items() if v}
cons={k:consensus(v) for k,v in alis.items()}
cons={k:v for k,v in cons.items() if v}
if len(cons)<2:
    print("  %-12s %d alignment(s) -- nothing to compare"%(a.key,len(cons))); sys.exit(0)

tmp=tempfile.mkdtemp()
open(tmp+"/c.faa","w").write("".join(">%s\n%s\n"%(k,v) for k,v in cons.items()))
subprocess.run(["makeblastdb","-in",tmp+"/c.faa","-dbtype","prot","-out",tmp+"/db"],check=True,stdout=subprocess.DEVNULL)
o=subprocess.run(["blastp","-query",tmp+"/c.faa","-db",tmp+"/db","-outfmt",
                  "6 qseqid sseqid pident length qlen slen","-evalue","1e-3",
                  "-max_target_seqs","50"],capture_output=True,text=True).stdout

# best alignment length seen for each unordered pair, from either direction
pair=collections.defaultdict(lambda:(0.0,0))
for l in o.splitlines():
    q,s,pid,ln,ql,sl=l.split("\t")
    if q==s: continue
    ln=int(ln); pid=float(pid)
    k=tuple(sorted((q,s)))
    if ln>pair[k][1]: pair[k]=(pid,ln)

drop={}
for (x,y),(pid,ln) in pair.items():
    lx,ly=len(cons[x]),len(cons[y])
    short,lng=(x,y) if lx<=ly else (y,x)
    ls,ll=min(lx,ly),max(lx,ly)
    if pid>=a.pident and ln/ls>=a.qcov and ls/ll<a.lenratio:
        prev=drop.get(short)
        if prev is None or ll>len(cons[prev[0]]):
            drop[short]=(lng,pid,ln/ls,ls/ll)

print("  %-12s %d alignments; %d contained"%(a.key,len(cons),len(drop)))
for k,(lng,pid,qc,lr) in sorted(drop.items()):
    print("      %-6s %4d aa  is %.1f%% identical over %.0f%% of itself to %-6s (%d aa, len ratio %.2f)  [%d seqs]"
          %(k,len(cons[k]),pid,100*qc,lng,len(cons[lng]),lr,len(alis[k])))
if a.write and drop:
    td=os.path.join(base,"truncated_alis"); os.makedirs(td,exist_ok=True)
    pd=os.path.join(base,"pssms")
    for k,(lng,pid,qc,lr) in drop.items():
        os.rename(os.path.join(home[k],k+".fa"), os.path.join(td,k+".fa"))
        for fn in list(os.listdir(pd)):
            if fn.endswith(".%s.pssm"%k):
                os.rename(os.path.join(pd,fn), os.path.join(td,fn))
    with open(os.path.join(base,"Truncation_Report_symmetric"),"w") as f:
        f.write("cluster\tlen\tcontained_in\tlen\tpercent_id\tcoverage_of_shorter\tlen_ratio\tn_seqs\n")
        for k,(lng,pid,qc,lr) in sorted(drop.items()):
            f.write("%s\t%d\t%s\t%d\t%.1f\t%.3f\t%.3f\t%d\n"%(k,len(cons[k]),lng,len(cons[lng]),pid,qc,lr,len(alis[k])))
    print("      retired %d alignment(s) to truncated_alis/"%len(drop))
