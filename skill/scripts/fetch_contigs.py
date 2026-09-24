#!/usr/bin/env python3
"""Download per-genome contig FASTAs for a taxon, batched.

Replaces `New-annotate-viral-taxon.pl -download-only`, which issues one request
per genome and (per SKILL_FEEDBACK) fails silently. Writes one <genome_id>.fna
per genome plus a metadata TSV in the shape merge_segments.py expects.
"""
import os, sys, time, argparse, collections
import requests
from concurrent.futures import ThreadPoolExecutor

URL="https://www.bv-brc.org/api/genome_sequence/"
HDR={"Content-Type":"application/rqlquery+x-www-form-urlencoded","Accept":"application/json"}

def fetch(gids, retries=5):
    q="in(genome_id,(%s))&select(genome_id,accession,sequence,length,description)&limit(25000)"%",".join(gids)
    for a in range(retries):
        try:
            r=requests.post(URL,data=q,headers=HDR,timeout=300)
            if r.status_code==200: return r.json()
            sys.stderr.write("HTTP %d\n"%r.status_code)
        except Exception as e: sys.stderr.write("ERR %s\n"%e)
        time.sleep(2*(a+1))
    raise RuntimeError("batch failed")

ap=argparse.ArgumentParser()
ap.add_argument("--meta", required=True, help="Alsu.full.tsv-shaped file already filtered to the taxon")
ap.add_argument("--out", required=True)
ap.add_argument("--min", type=int, default=1000)
ap.add_argument("--max", type=int, default=50000)
ap.add_argument("--batch", type=int, default=200)
ap.add_argument("--jobs", type=int, default=10)
a=ap.parse_args()

rows=[]
for l in open(a.meta):
    p=l.rstrip("\n").split("\t")
    try: ln=int(p[6])
    except: continue
    if a.min<=ln<=a.max: rows.append(p)
gids=[p[0] for p in rows]
meta={p[0]:p for p in rows}
os.makedirs(a.out,exist_ok=True)
batches=[gids[i:i+a.batch] for i in range(0,len(gids),a.batch)]
sys.stderr.write("%d genomes in %d batches\n"%(len(gids),len(batches)))

n=0; got=collections.defaultdict(list)
with ThreadPoolExecutor(max_workers=a.jobs) as ex:
    for res in ex.map(fetch,batches):
        for x in res: got[x["genome_id"]].append(x)
        n+=1
        if n%20==0: sys.stderr.write("  %d/%d\n"%(n,len(batches)))
#  Deduplicate contig ids within a genome. BV-BRC occasionally returns the
#  same sequence record twice -- 3 of 2,006 Pestivirus genomes and 3 of 670
#  Pegivirus -- and rast-create-genome then dies with
#      Error -32603 invoking add_contigs: Attempt to add duplicate contig id
#  losing the whole genome rather than the duplicate. Keyed on accession
#  because that is what becomes the contig id below; a genuinely different
#  segment has a different accession and is kept.
w=0; deduped=0
for gid,recs in got.items():
    seen=set(); keep=[]
    for x in recs:
        acc=x.get("accession",gid)
        if acc in seen: deduped+=1; continue
        seen.add(acc); keep.append(x)
    with open(os.path.join(a.out,gid+".fna"),"w") as f:
        for x in keep:
            f.write(">%s   %s\n%s\n"%(x.get("accession",gid),x.get("description","") or meta[gid][1],x["sequence"]))
    w+=1
if deduped: sys.stderr.write("dropped %d duplicate contig record(s)\n"%deduped)
missing=[g for g in gids if g not in got]
with open(a.out+".metadata","w") as f:
    f.write("genome_id\tgenome_name\tfamily\tgenus\tspecies\tlength\taccession\n")
    for gid in gids:
        if gid not in got: continue
        p=meta[gid]
        f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(gid,p[1],p[2],p[3],p[4],p[6],
                 ",".join(x.get("accession","") for x in got[gid])))
sys.stderr.write("wrote %d genome files, %d requested had no sequence\n"%(w,len(missing)))
