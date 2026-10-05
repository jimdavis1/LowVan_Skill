#!/usr/bin/env python3
"""Collections + measured length windows for one picornavirus genus.

Windows come from the per-species medians inside the genus, padded 20%, not
from the observed extremes -- the tails are fragments at one end and uncleaved
precursors at the other.
"""
import collections,re,sys,os,json,importlib
G=re.compile(r"^fig\|(\d+\.\d+)\.(\w+)\.\d+$")
POLY=re.compile(r"^(polyprotein|genome polyprotein|ORF1 protein|ORF1|putative genome polyprotein|"
                r"precursor polyprotein|polyprotein precursor|\w+ polyprotein|POLYPROTEIN|"
                r"structural polyprotein|polyprotein \(?P\d+\)?|polyprotein P1)$",re.I)
PREC_BASE=(r"^(P1|P2|P3|P1-P2|VP4-2C|P1 protein|P1 polyprotein|capsid polyprotein|"
           r"capsid proteins? (P1 )?precursor( P\d+)?|%s" r"capsid|capsid proteins?)$")
#  VP0 is the uncleaved VP4+VP2 PRECURSOR in the genera that cleave it, and the
#  MATURE product in those that do not. Treating it as a precursor globally
#  gutted Kobuvirus (18 of 200 VP0 records kept) and Parechovirus (269 of 305).
def prec_re(vp0_is_mature):
    return re.compile(PREC_BASE % ("" if vp0_is_mature else "VP0|"), re.I)

def main(genus, modname, feats, vp0_mature=False):
    PREC = prec_re(vp0_mature)
    sys.path.insert(0,"rules"); R=importlib.import_module(modname)
    gids={l.split("\t")[0] for l in open("meta/%s.tsv"%genus)}
    sp={}
    for l in open("meta/%s.tsv"%genus):
        f=l.rstrip("\n").split("\t"); sp[f[0]]=f[4] if len(f)>4 else "?"
    seq=dict(l.rstrip("\n").split("\t")[:2] for l in open("Picorna.uniq.seq") if "\t" in l)
    md5=[l.strip() for l in open("Picorna.uniq.md5")]
    ann=[l.rstrip("\n").split("\t") for l in open("Picorna.uniq.id_ann")]
    rows=collections.defaultdict(list); pool=[]; prec=[]
    for i,m in enumerate(md5):
        f=ann[i]; mm=G.match(f[0] if f else "")
        if not mm or mm.group(1) not in gids or not seq.get(m): continue
        a=(f[1] if len(f)>1 else "").strip(); s=seq[m]; rec=(f[0],a,s,sp.get(mm.group(1),"?"))
        if PREC.match(a): prec.append(rec); continue
        k="POLY" if POLY.match(a) else R.triage(a)[0]
        (rows[k] if k else pool).append(rec)
    # windows from per-species medians
    win={}
    for k in feats:
        bysp=collections.defaultdict(list)
        for _,_,s,species in rows.get(k,[]): bysp[species].append(len(s))
        meds=[sorted(v)[len(v)//2] for v in bysp.values() if len(v)>=3]
        if not meds:
            allL=sorted(len(r[2]) for r in rows.get(k,[]))
            if not allL: continue
            meds=[allL[len(allL)//2]]
        win[k]=(max(1,int(min(meds)*0.80)), int(max(meds)*1.20))
    os.makedirs("collections/%s"%genus,exist_ok=True)
    print("%-6s %7s %7s  %-12s %s"%("feat","kept","outlier","window","per-species medians"))
    for k in feats:
        if k not in win: print("%-6s   NOTHING"%k); continue
        lo,hi=win[k]; keep=[r for r in rows[k] if lo<=len(r[2])<=hi]
        with open("collections/%s/%s.fasta"%(genus,k),"w") as fh:
            for fid,a,s,_ in keep: fh.write(">%s %s\n%s\n"%(fid,a,s))
        bysp=collections.defaultdict(list)
        for _,_,s,species in rows[k]: bysp[species].append(len(s))
        meds=sorted(sorted(v)[len(v)//2] for v in bysp.values() if len(v)>=3)
        print("%-6s %7d %7d  %-12s %s"%(k,len(keep),len(rows[k])-len(keep),
              "%d-%d"%(lo,hi), "%d-%d"%(meds[0],meds[-1]) if meds else "-"))
    for nm,v in (("PRECURSORS",prec),("POOL",pool)):
        with open("collections/%s/%s.fasta"%(genus,nm),"w") as fh:
            for fid,a,s,_ in v: fh.write(">%s %s\n%s\n"%(fid,a,s))
    json.dump(win,open("rules/%s_windows.json"%genus.lower(),"w"),indent=1)
    print("\nprecursors held out %d | pool %d"%(len(prec),len(pool)))
    print("pool:",collections.Counter(r[1] for r in pool).most_common(6))

if __name__=="__main__":
    genus,mod=sys.argv[1],sys.argv[2]; feats=sys.argv[3].split(",")
    #  4th arg: "VP0MATURE" for the genera whose capsid stays as VP0
    main(genus,mod,feats,len(sys.argv)>4 and sys.argv[4]=="VP0MATURE")