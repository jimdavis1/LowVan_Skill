#!/usr/bin/env python3
"""Every number the five coverage audits use, from measurement files only -> pages/<M>.measure.json"""
import json,glob,os,csv,re,collections,sys
W=os.path.expanduser('~/lowvan_scratch/adeno_v3')
mods=['Siadenovirus','Barthadenovirus','Aviadenovirus','Mastadenovirus_MultiFiber','Mastadenovirus']
rt=list(csv.DictReader(open(f'{W}/audit/full/routing.tsv'),delimiter='\t'))
reroute={os.path.basename(f)[:-8]:json.load(open(f)).get('viral_family') for f in glob.glob(f'{W}/panel/reroute.anno/*.ann.gto')}
gm={r['genome']:r for r in csv.DictReader(open(f'{W}/fclass/genome_module.tsv'),delimiter='\t')}
cl={}
for l in open(f'{W}/gclust/g70_cluster.tsv'):
    a,b=l.split()[:2]; cl[b.split('.con')[0]]=a.split('.con')[0]
for m in mods:
    M={}
    J=json.load(open(f'{W}/build/{m}/{m}_Viral_PSSM.json'))[m]; F=J['features']; FJ=json.load(open(f'{W}/build/{m}/features.json'))
    M['features']=len(F); M['special']=sorted(k for k,v in F.items() if v.get('special'))
    M['mat']=sorted(k for k,v in F.items() if v.get('feature_type')=='mat_peptide')
    M['essential']=sorted(k for k,v in F.items() if v.get('copy_num'))
    M['uncharacterized']=sorted(k for k,v in F.items() if v['anno'].startswith('Uncharacterized'))
    M['pssms']=sum(len(glob.glob(os.path.expanduser(f'~/lowvan_scratch/VA_v3/Viral-PSSMs/{m}.pssms/{k}/*'))) for k in F)
    M['reps']=len(J['close_genomes'])
    M['cutoffs']=[(k,v.get('bit_cutoff'),v.get('cutoff_source','')) for k,v in sorted(FJ.items()) if v.get('cutoff_source')]
    asg=[r for r in gm.values() if r['module']==m]; M['genomes']=len(asg)
    M['clusters70']=len({cl.get(r['genome'],r['genome']) for r in asg})
    M['seq_genus']=collections.Counter(r['seq_genus'] for r in asg).most_common()
    # routing (full panel, with reroute applied)
    rows=[r for r in rt if r['panel']==m]
    routed=lambda r: reroute.get(r['genome'],r['routed'])
    M['routing']={'panel':len(rows),'own_before':sum(r['routed']==m for r in rows),'own_after':sum(routed(r)==m for r in rows),
                  'other_after':collections.Counter(routed(r) for r in rows if routed(r)!=m)}
    # per-feature call rates over near-complete routed genomes (full panel)
    G=[]
    for f in glob.glob(f'{W}/panel/{m}.anno/*.ann.gto'):
        g=json.load(open(f))
        gid=os.path.basename(f)[:-8]
        if gid in reroute: g=json.load(open(f'{W}/panel/reroute.anno/{gid}.ann.gto'))
        if g.get('viral_family')!=m: continue
        G.append((sum(len(c['dna']) for c in g['contigs']),{x['family_assignments'][0][1].split('.')[1] for x in g['features'] if x.get('family_assignments')}))
    med=sorted(L for L,_ in G)[len(G)//2]; NC=[s for L,s in G if L>=0.9*med]
    M['nearcomplete']=len(NC); M['routed']=len(G); M['median_len']=med
    M['call_rate']={k:round(100*sum(k in s for s in NC)/len(NC),1) for k in F}
    # quality (subset, rules)
    q=collections.Counter(); fl=collections.Counter()
    for f in glob.glob(f'{W}/panel/{m}.rules.anno/*.qual.gto'):
        g=json.load(open(f)); bad=bool(g.get('genome_quality_flags')) or any(c.get('contig_quality_flags') for c in g.get('contigs',[]))   # run_gto_eval's 'good' = no genome AND no contig flags
        q['poor' if bad else 'good']+=1
        for x in g.get('genome_quality_flags',[])+[y for c in g.get('contigs',[]) for y in (c.get('contig_quality_flags') or [])]: fl[x.split(':')[0].split(',')[0]]+=1
    M['quality']=dict(q); M['quality_flags']=fl.most_common()
    # registry
    R=json.load(open(f'{W}/build/{m}/registry.json'))
    M['registry']=R if len(json.dumps(R))<200000 else None
    # extras / stacks / cds stop
    M['calls_audit']=[l.rstrip() for l in open(f'{W}/audit/full/calls_audit_{m}.out') if 'extra calls' in l and 'single-copy' not in l]
    M['cds_stop_dext']=[l.rstrip() for l in open(f'{W}/audit/full/cds_stop.dext.txt') if l.startswith(m+' ')]
    M['matpep_ends_dext']=[l.rstrip() for l in open(f'{W}/audit/full/matpep_ends.dext.txt') if l.startswith(m+' ')]
    # vocabulary
    agg=json.load(open(f'{W}/build/{m}/agg.json')); unb=json.load(open(f'{W}/build/{m}/unbinned.json'))
    M['collapse']={'targets':len(agg),'strings':sum(a['nstrings'] for a in agg),'features':sum(a['total'] for a in agg),
                   'unbinned_strings':len(unb),'unbinned_features':sum(u.get('count',u.get('n',0)) if isinstance(u,dict) else 0 for u in unb)}
    rare=[l.strip() for l in open(f'{W}/pages/logs/{m}.rare.log') if 'fold' in l or 'pooled' in l]; M['rarefaction']=rare
    syn=collections.Counter()
    for l in list(open(f'{W}/build/{m}/synonyms.tsv'))[1:]:
        c,a,g,d=l.rstrip('\n').split('\t'); syn[d if d in ('UNASSIGNED','DISCARDED','NOT_MODELLED') else 'bound']+=int(c)
    M['synonyms']=dict(syn)
    M['top_unassigned']=[]
    for l in list(open(f'{W}/build/{m}/synonyms.tsv'))[1:]:
        c,a,g,d=l.rstrip('\n').split('\t')
        if d=='UNASSIGNED' and len(M['top_unassigned'])<12: M['top_unassigned'].append((a,int(c)))
    json.dump(M,open(f'{W}/pages/{m}.measure.json','w'),indent=1,default=list)
    print(m,'genomes',M['genomes'],'clusters',M['clusters70'],'routed',M['routing'],'nearcomplete',M['nearcomplete'],'quality',M['quality'],'collapse',M['collapse'],'syn',M['synonyms'])
