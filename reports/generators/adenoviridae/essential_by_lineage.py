#!/usr/bin/env python3
"""copy_num candidates weighted by LINEAGE, not by genome count.
For each module: every feature's call rate per 70% nt cluster (gclust/g70_cluster.tsv),
over full-panel genomes routed to the module (panel/<M>.anno/*.ann.gto). A feature is a
lineage-wide candidate if the MEDIAN... no: if it is called in >= 85% of genomes in >= 85%
of clusters with >= 2 genomes (singleton clusters counted too, at their one genome).
Also flags features whose annotation string is shared with another feature (never essential).
usage: essential_by_lineage.py <Module>"""
import json,glob,os,sys,collections
m=sys.argv[1]; W=os.path.expanduser('~/lowvan_scratch/adeno_v3')
cl={}
for l in open(f'{W}/gclust/g70_cluster.tsv'):
    rep,mem=l.split()[:2]; cl[mem.split('.con')[0].split('|')[0]]=rep.split('.con')[0].split('|')[0]
F=json.load(open(f'{W}/build/{m}/{m}_Viral_PSSM.json'))[m]['features']
anno_n=collections.Counter(v['anno'] for v in F.values())
calls=collections.defaultdict(lambda: collections.defaultdict(int)); gen=collections.defaultdict(int)
miss=0; G=[]
for f in glob.glob(f'{W}/panel/{m}.anno/*.ann.gto'):
    g=json.load(open(f))
    if g.get('viral_family')!=m: continue
    G.append((f,g,sum(len(c['dna']) for c in g['contigs'])))
med=sorted(x[2] for x in G)[len(G)//2]; nfull=0
for f,g,L in G:
    if L<0.9*med: continue          # near-complete genomes only: partial deposits are not lineages
    nfull+=1
    gid=os.path.basename(f)[:-8]; c=cl.get(gid)
    if c is None: miss+=1; c='NA:'+gid
    gen[c]+=1
    for k in {x['family_assignments'][0][1].split('.')[1] for x in g['features'] if x.get('family_assignments')}:
        calls[c][k]+=1
print(f'{m}: median {med} nt, {nfull}/{len(G)} near-complete; {sum(gen.values())} genomes in {len(gen)} clusters ({miss} without a cluster id)')
print(f'  {"feature":12s} copy_num  genomes%  clusters>=50%   shared-anno')
for k in sorted(F, key=lambda k: -sum(calls[c][k] for c in gen)):
    v=F[k]; tot=sum(gen.values()); gr=sum(calls[c][k] for c in gen)/tot
    cr=sum(1 for c in gen if calls[c][k]/gen[c]>=0.5)/len(gen)
    sh='SHARED' if anno_n[v['anno']]>1 else ''
    verdict='keep' if v.get('copy_num') and cr>=0.85 and not sh else ('DROP' if v.get('copy_num') else '')
    if v.get('copy_num') or cr>=0.85:
        print(f'  {k:12s} {str(v.get("copy_num")):8s} {100*gr:6.1f}   {100*cr:6.1f}        {sh:7s} {verdict}')
