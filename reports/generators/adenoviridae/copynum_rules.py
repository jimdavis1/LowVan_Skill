#!/usr/bin/env python3
"""Filter apply_copy_num's essential list by the curator's rules (16 Nov 2026).
EXCLUDE: uncharacterized (or shared annotation string), special, lineage-specific
(called in >=50% of near-complete genomes in <85% of 70% clusters), and list every
'struggle' (explicit/retuned cutoff, extra calls, inseparable, same-frame stacks with a
non-precursor feature) for review. usage: copynum_rules.py <Module>"""
import json,re,os,sys,subprocess,collections,csv
m=sys.argv[1]; W=os.path.expanduser('~/lowvan_scratch/adeno_v3')
J=json.load(open(f'{W}/build/{m}/{m}_Viral_PSSM.json'))[m]['features']
FJ=json.load(open(f'{W}/build/{m}/features.json.pre_copynum'))
pairs={tuple(x) for x in json.load(open(next(os.path.join(r,f) for r,_,fs in os.walk(f'{W}/build/{m}/collections') for f in fs if f=='protected_pairs.json')))}
annon=collections.Counter(v['anno'] for v in J.values())
# lineage presence
lin={}
for l in subprocess.run(['python3',f'{W}/pages/essential_by_lineage.py',m],capture_output=True,text=True).stdout.splitlines():
    p=l.split()
    if len(p)>=4 and p[0] in J:
        try: lin[p[0]]=float(p[3])
        except: pass
# extras / inseparable (full panel)
ext={}
for l in open(f'{W}/audit/full/calls_audit_{m}.out'):
    mm=re.match(r'(\S+)\s+genomes\s+\d+\s+extra calls\s+(\d+).*?(INSEPARABLE|SEPARABLE)?',l)
    if mm and int(mm.group(2))>0: ext[mm.group(1)]=f'{mm.group(2)} extra' + (' INSEPARABLE' if 'INSEPARABLE' in l else '')
# stacks (dext) excluding product-on-own-precursor
stk=collections.defaultdict(list)
for l in open(f'{W}/audit/full/stacks_{m}.dext.out'):
    p=l.split()
    if len(p)>5 and p[0] in J and 'on:' in l:
        for on in l.split('on:')[1].split(','):
            o=on.strip().split('x')[0]
            if (p[0],o) in pairs or (o,p[0]) in pairs: continue
            stk[p[0]].append(on.strip())
rows=[]
for k,v in sorted(J.items()):
    if not v.get('copy_num'): continue
    excl=[]; rev=[]
    if v['anno'].startswith('Uncharacterized'): excl.append('uncharacterized')
    elif annon[v['anno']]>1: excl.append('shared annotation string')
    if v.get('special'): excl.append('special')
    if k in lin and lin[k]<85: excl.append(f'lineage-specific ({lin[k]:.0f}% of lineages)')
    f=FJ.get(k,{})
    if 'bit_cutoff' in f: rev.append(f"cutoff set {f['bit_cutoff']}" + (f" ({f['cutoff_source'][:40]})" if f.get('cutoff_source') else ''))
    if 'coverage_cutoff' in f and f['coverage_cutoff']!=0.65: rev.append(f"coverage {f['coverage_cutoff']}")
    if k in ext: rev.append(ext[k])
    if stk.get(k): rev.append('stacks '+';'.join(stk[k]))
    print(f"{k:11s} {'EXCLUDE: '+', '.join(excl) if excl else ('REVIEW: '+' | '.join(rev) if rev else 'keep')}" + (f"   [also: {' | '.join(rev)}]" if excl and rev else ''))
