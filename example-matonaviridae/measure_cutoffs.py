#!/usr/bin/env python3
"""bit_cutoff window per feature, measured at the protein level.

The rule (references/json-schema.md): the cutoff sits ABOVE the noise ceiling
-- the best score a feature's profiles reach on sequence that is not that
feature -- and BELOW the weakest true signal.

"Not that feature" needs care in a polyprotein taxon.  A mature peptide lies
inside its precursor, so P150's profile hitting a P200 sequence is ON TARGET,
not noise.  Only cross-ORF and sibling-peptide hits are noise:

    P200 contains P150, P90        STP contains C, E2, E1

so (P200,P150) is related and (P200,C) is not.
"""
import glob,os,re,subprocess,tempfile,collections,json
FEATS=['P200','P150','P90','P110','C','E2','E1']
CONTAINS={'P200':{'P150','P90'},'P110':{'C','E2','E1'}}
RELATED=set()
for w,ps in CONTAINS.items():
    for p in ps:
        RELATED.add((w,p)); RELATED.add((p,w))
gn={}
for l in open('Matona.id_name_gs'):
    p=l.rstrip('\n').split('\t'); gn[p[0]]=p[1]
def sp(fid):
    m=re.match(r'(?:>)?(fig\|(\d+\.\d+))',fid)
    n=gn.get(m.group(2),'?') if m else '?'
    if n.startswith('Rubella') or n.startswith('Rubivirus rubellae'): return 'rubella'
    if n.startswith('Rubivirus strelense') or n.startswith('Rustrela'): return 'Rustrela'
    if n.startswith('Ruhugu'): return 'Ruhugu'
    return 'other'
def readfa(p):
    d={};h=None
    for line in open(p):
        if line.startswith('>'): h=line[1:].split()[0]; d[h]=''
        elif h: d[h]+=line.strip().replace('-','')
    return d
# ---- every sequence, tagged with its feature
seqs_of={}
for f in FEATS:
    b=f'Alignments/Matonaviridae/{f}'
    d={}
    for a in glob.glob(f'{b}/corrected_alis/*.fa')+glob.glob(f'{b}/reclustered_alis/*.fa'):
        d.update(readfa(a))
    lo=f'{b}/Leftover_Seqs.aa'
    if os.path.exists(lo) and os.path.getsize(lo): d.update(readfa(lo))
    seqs_of[f]=d
tmp=tempfile.mkdtemp()
allfa=f'{tmp}/all.faa'
owner={}
with open(allfa,'w') as fh:
    for f,d in seqs_of.items():
        for k,v in d.items():
            tag=f"{f}@@{k}"
            owner[tag]=(f,k)
            fh.write(f">{tag}\n{v}\n")
# ---- score every profile against everything
best=collections.defaultdict(dict)      # feat -> {tag: bits}  (best over its profiles)
for f in FEATS:
    for pssm in sorted(glob.glob(f'Alignments/Matonaviridae/{f}/pssms/*.pssm')):
        out=subprocess.run(['psiblast','-in_pssm',pssm,'-subject',allfa,
                            '-outfmt','6 sseqid bitscore','-comp_based_stats','0',
                            '-evalue','100','-max_target_seqs','100000'],
                           capture_output=True,text=True).stdout
        for line in out.splitlines():
            s,b=line.split('\t'); b=float(b)
            if b>best[f].get(s,0): best[f][s]=b
res={}
print(f"{'feat':<6}{'nseq':>5}{'weakest true':>14}{'noise ceiling':>15}{'window':>10}   {'uncovered':>9}")
print("-"*76)
for f in FEATS:
    own=[t for t in owner if owner[t][0]==f]
    sig=[(best[f].get(t,0.0),t) for t in own]
    noise=[(best[f].get(t,0.0),t) for t in owner
           if owner[t][0]!=f and (f,owner[t][0]) not in RELATED]
    sig.sort(); noise.sort()
    uncov=[t for b,t in sig if b==0]
    weak=min((b for b,_ in sig if b>0), default=0)
    ceil_=max((b for b,_ in noise), default=0)
    res[f]=dict(weakest=weak,noise=ceil_,n=len(own),
                uncovered=[owner[t][1] for t in uncov],
                noise_top=[(owner[t][0],owner[t][1],b) for b,t in noise[-3:]])
    print(f"{f:<6}{len(own):>5}{weak:>14.0f}{ceil_:>15.0f}{weak-ceil_:>10.0f}   {len(uncov):>9}")
print("-"*76)
for f in FEATS:
    r=res[f]
    if r['uncovered']:
        print(f"\n{f}: {len(r['uncovered'])} sequence(s) hit NO profile of their own feature:")
        for u in r['uncovered'][:8]:
            print(f"    {u}  [{sp(u)}]  {len(seqs_of[f][u])} aa")
    if r['noise_top']:
        print(f"{f}: top cross-feature (noise) hits: " +
              ', '.join(f"{a}/{c:.0f}" for a,_b,c in r['noise_top']))
json.dump(res,open('cutoff_measurements.json','w'),indent=1)
