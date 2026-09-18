#!/usr/bin/env python3
"""Coverage at a fixed rep-contig budget.

The wrong question is "how many references cover every genome" -- for a
divergent taxon that is one per species and the answer is always "too many".
The right question, and the one the project constrains, is: **with a budget of
N references chosen greedily largest-first, what fraction of the taxon's
genomes route?**

A long tail of singletons does not disqualify a module. It only means the tail
is uncovered, which is a documented gap. What disqualifies a module is failing
to reach most genomes within the budget.

    python3 repcontig_budget.py --family Virgaviridae --genera Tobamovirus \
            --min-len 6000 --budget 25
"""
import os,sys,subprocess,tempfile,argparse,collections,requests

ap=argparse.ArgumentParser()
ap.add_argument("--family",required=True)
ap.add_argument("--genera",default="")
ap.add_argument("--min-len",type=int,required=True)
ap.add_argument("--budget",type=int,default=25)
ap.add_argument("--mi",type=float,default=0.7)
ap.add_argument("--label",default=None)
ap.add_argument("--meta",default="Alsu.full.tsv")
a=ap.parse_args()
gen=set(x for x in a.genera.split(",") if x)

rows=[]
for l in open(a.meta):
    p=l.rstrip("\n").split("\t")
    if len(p)<10 or p[2]!=a.family: continue
    if gen and (p[3] or "") not in gen: continue
    try: ln=int(p[6])
    except: continue
    if ln>=a.min_len: rows.append(p[0])
if not rows: sys.exit("none")
label=a.label or (a.family if not gen else "+".join(sorted(gen)))

seqs={}
for i in range(0,len(rows),150):
    q="in(genome_id,(%s))&select(genome_id,sequence)&limit(3000)"%",".join(rows[i:i+150])
    r=requests.post("https://www.bv-brc.org/api/genome_sequence/",data=q,
      headers={"Content-Type":"application/rqlquery+x-www-form-urlencoded","Accept":"application/json"},timeout=300)
    for x in r.json(): seqs[x["genome_id"]]=seqs.get(x["genome_id"],"")+x["sequence"]

tmp=tempfile.mkdtemp(); fa=os.path.join(tmp,"g.fna")
with open(fa,"w") as f:
    for g,s in seqs.items(): f.write(">%s\n%s\n"%(g,s))
pre=os.path.join(tmp,"c")
subprocess.run(["mmseqs","easy-cluster",fa,pre,os.path.join(tmp,"t"),
                "--min-seq-id",str(a.mi),"-c","0.5","--cov-mode","1","--threads","6"],
               capture_output=True)
clu=collections.defaultdict(list)
for l in open(pre+"_cluster.tsv"):
    r_,m=l.rstrip("\n").split("\t"); clu[r_].append(m)
sizes=sorted((len(v) for v in clu.values()),reverse=True)
N=len(seqs); tot=len(sizes)
def cum(k): return sum(sizes[:k])
marks=[5,10,15,20,25,30,40,50]
print("%-26s %5d genomes, %4d clusters at %.0f%% identity"%(label,N,tot,100*a.mi))
print("%-26s coverage: %s"%("", "  ".join("%d:%.0f%%"%(k,100*cum(k)/N) for k in marks if k<=max(marks))))
c=100*cum(a.budget)/N
verdict = "SHIP" if c>=80 else ("MARGINAL" if c>=60 else "TABLE")
print("%-26s  -> %d references cover %.1f%% of genomes   %s"%("",a.budget,c,verdict))
print("%-26s     largest clusters: %s"%("", ", ".join(str(x) for x in sizes[:8])))
