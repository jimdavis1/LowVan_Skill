#!/usr/bin/env python3
"""synonyms.tsv for a v3 adeno module (count/annotation/genus/bins_to), one row per
(deposited string, genus, destination), counting deposited FEATURES.
Destination = the shipped feature whose final collection holds the sequence (md5);
else the homology verdict's feature from <M>/ASSIGN.tsv if that key is declared;
PARTIAL_BY_NAME or not in the pool (submitter-marked partial) -> DISCARDED;
other pool sequences -> UNASSIGNED.  usage: make_synonyms.py <Module>"""
import os,sys,csv,glob,hashlib,collections,json
m=sys.argv[1]; W=os.path.expanduser('~/lowvan_scratch/adeno_v3'); A=os.path.expanduser('~/lowvan_scratch/adeno')
sys.path.insert(0,'/Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Alsuviricetes/LowVan_Skill/skill/scripts')
import synonyms_table
mod={r['genome']:r for r in csv.DictReader(open(f'{W}/fclass/genome_module.tsv'),delimiter='\t')}
J=json.load(open(f'{W}/build/{m}/{m}_Viral_PSSM.json'))[m]['features']
def seqs(p):
    h=None;s=[]
    for l in open(p):
        if l.startswith('>'):
            if h is not None: yield h,''.join(s)
            h=l[1:].strip();s=[]
        else: s.append(l.strip())
    if h is not None: yield h,''.join(s)
final={}
for f in glob.glob(f'{W}/build/{m}/collections/{m}/*.fasta'):
    k=os.path.basename(f)[:-6]
    if k not in J: continue
    for h,s in seqs(f): final.setdefault(hashlib.md5(s.upper().replace('*','').encode()).hexdigest(),k)
assign={}
for r in csv.DictReader(open(f'{W}/{m}/ASSIGN.tsv'),delimiter='\t'): assign[r['md5']]=(r['verdict'],r['feature'])
prod=dict(l.rstrip('\n').split('\t',1) for l in open(f'{A}/A.id_product') if '\t' in l)
C=collections.Counter(); hit=collections.Counter()
for l in open(f'{A}/A.id_md5'):
    p=l.split()
    if len(p)<2: continue          # feature with no protein md5 (RNA etc.)
    fid,md5=p[:2]
    gid=fid.split('|')[1].rsplit('.',2)[0]
    r=mod.get(gid)
    if not r or r['module']!=m: continue
    genus=r['seq_genus']
    if md5 in final: d=final[md5]; hit['final collection']+=1
    elif md5 in assign:
        v,k=assign[md5]
        if v=='PARTIAL_BY_NAME': d='DISCARDED'
        elif v=='ASSIGNED' and k in J: d=k; hit['assigned, later filtered']+=1
        elif v=='ASSIGNED': d='NOT_MODELLED'
        else: d='UNASSIGNED'
    else: d='DISCARDED'
    hit[d if d in ('DISCARDED','UNASSIGNED','NOT_MODELLED') else 'bound']+=1
    C[(prod.get(fid,''),genus,d)]+=1
out=f'{W}/build/{m}/synonyms.tsv'
st=synonyms_table.write(out,[(n,a,g,d) for (a,g,d),n in C.items()])
print(m,dict(hit),st)
