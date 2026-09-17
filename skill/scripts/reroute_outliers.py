#!/usr/bin/env python3
"""Re-route length-outliers that are really a different feature.

A text rule can only read the annotation, and some annotations are simply wrong
about which protein the sequence is. The length gate catches the consequence --
a 122 aa "nonstructural protein" is not a 1700 aa polyprotein, so it lands in
ORF1.outliers.fasta -- but it leaves the sequence parked under the wrong
feature, where qc_cross_feature.py then reports it as a MISLABEL.

This asks the cheap question the length gate could not: does the outlier match
some OTHER feature's in-range collection? Same threshold as the homology rescue
(>= 80% identity over >= 60% of the query), because SKILL.md requires the
adoption threshold to be identical everywhere or a sequence gets claimed twice.

Only moves an outlier when the better match is a different feature AND the
sequence falls inside that feature's length window. Everything else stays put.
"""
import os, sys, subprocess, tempfile, argparse, collections

def read_fa(p):
    out=[];h=None;s=[]
    for l in open(p):
        if l.startswith(">"):
            if h is not None: out.append((h,"".join(s)))
            h=l[1:].rstrip("\n");s=[]
        else: s.append(l.strip())
    if h is not None: out.append((h,"".join(s)))
    return out

ap=argparse.ArgumentParser()
ap.add_argument("--workdir",default=".")
ap.add_argument("--module",required=True)
ap.add_argument("--expected",required=True)
ap.add_argument("--min-id",type=float,default=80.0)
ap.add_argument("--min-qcov",type=float,default=0.60)
ap.add_argument("--write",action="store_true")
a=ap.parse_args()

EXP={}
for part in a.expected.split(","):
    k,r=part.split(":"); lo,hi=r.split("-"); EXP[k]=(int(lo),int(hi))
cdir=os.path.join(a.workdir,"collections",a.module)

inrange={k:read_fa(os.path.join(cdir,k+".fasta")) for k in EXP
         if os.path.exists(os.path.join(cdir,k+".fasta"))}
tmp=tempfile.mkdtemp()
with open(tmp+"/s.faa","w") as f:
    for k,v in inrange.items():
        for h,s in v: f.write(">%s@@%s\n%s\n"%(k,h.split()[0],s))
subprocess.run(["makeblastdb","-in",tmp+"/s.faa","-dbtype","prot","-out",tmp+"/db"],
               check=True,stdout=subprocess.DEVNULL)

moves=collections.defaultdict(list); summary=collections.Counter()
for k in EXP:
    op=os.path.join(cdir,k+".outliers.fasta")
    if not os.path.exists(op): continue
    recs=read_fa(op)
    with open(tmp+"/q.faa","w") as f:
        for i,(h,s) in enumerate(recs): f.write(">q%d\n%s\n"%(i,s))
    o=subprocess.run(["blastp","-query",tmp+"/q.faa","-db",tmp+"/db","-outfmt",
                      "6 qseqid sseqid pident length qlen","-evalue","1e-5",
                      "-max_target_seqs","5","-num_threads","8"],
                     capture_output=True,text=True).stdout
    best={}
    for l in o.splitlines():
        q,s,pid,ln,ql=l.split("\t"); pid=float(pid); ln=int(ln); ql=int(ql)
        if pid<a.min_id or ln/ql<a.min_qcov: continue
        feat=s.split("@@")[0]
        sc=pid*ln
        if q not in best or sc>best[q][1]: best[q]=(feat,sc,pid)
    keep=[]
    for i,(h,s) in enumerate(recs):
        b=best.get("q%d"%i)
        if b and b[0]!=k:
            lo,hi=EXP[b[0]]
            if lo<=len(s)<=hi:
                moves[b[0]].append((h,s)); summary[(k,b[0])]+=1; continue
        keep.append((h,s))
    moves.setdefault("__keep__"+k,keep)

print("outliers re-routed to a different feature:")
if not summary: print("   none")
for (frm,to),n in summary.most_common():
    print("   %-12s -> %-12s %4d"%(frm,to,n))
if a.write and summary:
    for k in EXP:
        kk="__keep__"+k
        if kk in moves:
            p=os.path.join(cdir,k+".outliers.fasta")
            with open(p,"w") as f:
                for h,s in moves[kk]: f.write(">%s\n%s\n"%(h,s))
    for k,v in moves.items():
        if k.startswith("__keep__") or not v: continue
        with open(os.path.join(cdir,k+".fasta"),"a") as f:
            for h,s in v: f.write(">%s\n%s\n"%(h,s))
        print("   appended %d to %s.fasta"%(len(v),k))
    with open(os.path.join(a.workdir,"collections","REROUTED.tsv"),"w") as f:
        f.write("from_feature\tto_feature\tcount\n")
        for (frm,to),n in summary.most_common(): f.write("%s\t%s\t%d\n"%(frm,to,n))
elif summary:
    print("\n(dry run -- rerun with --write)")
