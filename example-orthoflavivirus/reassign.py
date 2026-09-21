#!/usr/bin/env python3
"""Identify what each junction-dropped record actually is, by homology.

Two searches, not one. Every mat_peptide is a substring of the polyprotein,
so including POLY in the target set for a mat_peptide query is actively
harmful: POLY hits align at ~100% identity over the full query while covering
3-4% of the target, they outscore the true assignment on bits, and with
14,095 polyproteins they fill any --max-seqs window. A 139-residue anchored
capsid then looks like a partial polyprotein instead of a complete ancC.

So mat_peptide queries go against the 14 mat_peptide collections only, and
polyprotein-length queries go against POLY only.

Coverage is filtered BEFORE ranking by score, for the same reason.
"""
import collections, os, subprocess, tempfile, sys

BOX = "/Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Alsuviricetes/work/Flaviviridae"
MAT = ["ANCHC","C","PRM","PR","M","E","NS1","NS2A","NS2B","NS3","NS4A","2K","NS4B","NS5"]
ALL = MAT + ["POLY"]
SMALL=set("GSATC"); BASIC=set("RK"); MET=set("M")
J={"ANCHC":(MET,SMALL),"C":(MET,BASIC),"PRM":(None,SMALL),"PR":(None,BASIC),
   "M":(None,SMALL),"E":(None,SMALL),"NS1":(None,SMALL),"NS2A":(None,BASIC),
   "NS2B":(SMALL,BASIC),"NS3":(SMALL,BASIC),"NS4A":(SMALL,BASIC),"2K":(SMALL,SMALL),
   "NS4B":(None,BASIC),"NS5":(SMALL,None),"POLY":(MET,None)}
LONG = 2000   # a query this long is a polyprotein, not a mature peptide

def rd(p):
    d={};h=None;s=[]
    for line in open(p):
        if line.startswith('>'):
            if h: d[h]=''.join(s)
            h=line[1:].strip();s=[]
        else: s.append(line.strip())
    if h: d[h]=''.join(s)
    return d

clean={k: rd(os.path.join(BOX,'collections/Orthoflavivirus/%s.fasta'%k)) for k in ALL}
pre={}
for k in ALL:
    for h,v in rd(os.path.join(BOX,'collections_preX/Orthoflavivirus/%s.fasta'%k)).items():
        pre[(k,h.split()[0])]=v

drop=[l.rstrip('\n').split('\t') for l in open(os.path.join(BOX,'collections/JUNCTION_DROPPED.tsv'))][1:]
qs=[]
for k,fid,ln,a,b,bad in drop:
    v=pre.get((k,fid))
    if v is not None: qs.append((k,fid,v))
print("%d junction-dropped records to identify" % len(qs), flush=True)

tmp=tempfile.mkdtemp(prefix='ra.')
def search(qidx, feats, tag):
    if not qidx: return {}
    tgt=os.path.join(tmp,'t_%s.faa'%tag)
    with open(tgt,'w') as fh:
        for k in feats:
            for i,(h,v) in enumerate(clean[k].items()):
                fh.write(">%s|%d\n%s\n"%(k,i,v))
    q=os.path.join(tmp,'q_%s.faa'%tag)
    with open(q,'w') as fh:
        for i in qidx: fh.write(">q%d\n%s\n"%(i,qs[i][2]))
    res=os.path.join(tmp,'h_%s.m8'%tag)
    print("  [%s] %d queries vs %d targets ..." % (tag,len(qidx),sum(len(clean[k]) for k in feats)), flush=True)
    r=subprocess.run(["mmseqs","easy-search",q,tgt,res,os.path.join(tmp,'w_%s'%tag),
        "--format-output","query,target,fident,bits,qcov,tcov",
        "-s","5.7","--max-seqs","60","-e","1e-3","--threads","8"],
        capture_output=True,text=True)
    if not os.path.exists(res):
        print("  [%s] search produced no output\n%s"%(tag,r.stderr[-800:]), flush=True); return {}
    qual={}
    for line in open(res):
        f=line.rstrip().split('\t'); qi=int(f[0][1:]); feat=f[1].split('|')[0]
        ident,bits,qc,tc=float(f[2]),float(f[3]),float(f[4]),float(f[5])
        if qc<0.80 or tc<0.80 or ident<0.50: continue      # coverage first
        if qi not in qual or bits>qual[qi][1]:             # then best score
            qual[qi]=(feat,bits,ident,qc,tc)
    print("  [%s] %d of %d got a coverage-qualified hit" % (tag,len(qual),len(qidx)), flush=True)
    return qual

short=[i for i,(k,f,v) in enumerate(qs) if len(v)<LONG]
long_ =[i for i,(k,f,v) in enumerate(qs) if len(v)>=LONG]
qual={}
qual.update(search(short, MAT, "mat"))
qual.update(search(long_, ["POLY"], "poly"))

moves=collections.Counter(); why=collections.Counter(); ex=collections.defaultdict(list); out=[]
for i,(k,fid,v) in enumerate(qs):
    if i not in qual: why['no coverage-qualified hit']+=1; continue
    feat,bits,ident,qc,tc=qual[i]
    if feat==k: why['best hit is its own feature (right protein, broken end)']+=1; continue
    fs,ls=J[feat]
    if (fs and v[0] not in fs) or (ls and v[-1] not in ls):
        why['fails the target feature junction rule']+=1; continue
    moves[(k,feat)]+=1; out.append((k,feat,fid,len(v),ident,qc,tc))
    if len(ex[(k,feat)])<3: ex[(k,feat)].append((fid,len(v),ident,qc,tc))

print("\n--- reassignments ---", flush=True)
for (a,b),c in moves.most_common():
    print("  %-6s -> %-6s %5d" % (a,b,c))
    for fid,l,ii,qc,tc in ex[(a,b)]:
        print("       %-30s len %4d ident %.2f qcov %.2f tcov %.2f" % (fid,l,ii,qc,tc))
print("  TOTAL MOVED: %d of %d" % (sum(moves.values()), len(qs)))
print("\n--- stay dropped ---")
for k,c in why.most_common(): print("  %-52s %5d" % (k,c))
with open(os.path.join(BOX,'collections/REASSIGNED.tsv'),'w') as fh:
    fh.write("from_feature\tto_feature\tfeature_id\tlength\tidentity\tqcov\ttcov\n")
    for r in sorted(out): fh.write("%s\t%s\t%s\t%d\t%.3f\t%.2f\t%.2f\n"%r)
print("\nDONE -- wrote collections/REASSIGNED.tsv (%d rows)" % len(out))
