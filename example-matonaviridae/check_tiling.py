#!/usr/bin/env python3
"""Prove the binning without reference to any database annotation.

If p150 and p90 really are the two products of p200, then in any genome
carrying all three the precursor must EQUAL the concatenation of its products.
Same for p110 = C + E2 + E1.  This is a sequence identity test, so it is
independent of whatever the submitter called anything.
"""
import collections, re, os
W=os.path.dirname(os.path.abspath(__file__))
md5=[l.rstrip('\n') for l in open(f'{W}/Matona.uniq.md5')]
seq=[l.rstrip('\n').split('\t')[1] if '\t' in l else '' for l in open(f'{W}/Matona.uniq.seq')]
m2s=dict(zip(md5,seq))
# which feature did build_collections assign each sequence to?
feat_of={}
for f in ['P200','P150','P90','P110','C','E2','E1']:
    for suf in ['','.partial','.outliers']:
        p=f'{W}/collections/Matonaviridae/{f}{suf}.fasta'
        if not os.path.exists(p): continue
        for line in open(p):
            if line.startswith('>'): cur=line
            else: pass
        # re-read properly
        h=None
        for line in open(p):
            if line.startswith('>'): h=line[1:].split()[0]; buf=[]
            elif h: feat_of.setdefault(h,f)
# genome -> {feat: seq}  using the id_md5 map so we get per-genome context
per=collections.defaultdict(dict)
for l in open(f'{W}/Matona.id_md5'):
    fid,m=l.rstrip('\n').split('\t')
    g=re.match(r'fig\|(\d+\.\d+)',fid)
    if not g or m not in m2s: continue
    f=feat_of.get(fid)
    if f: per[g.group(1)][f]=m2s[m]
def test(name, whole, parts):
    ok=bad=0; examples=[]
    for gid,d in per.items():
        if whole in d and all(p in d for p in parts):
            cat=''.join(d[p] for p in parts)
            if d[whole]==cat: ok+=1
            else:
                bad+=1
                if len(examples)<3:
                    examples.append((gid,len(d[whole]),len(cat)))
    print(f"  {name:<26} precursor == concat(products):  {ok} exact, {bad} mismatched")
    for e in examples: print(f"      mismatch {e[0]}  precursor {e[1]} aa vs concat {e[2]} aa")
    return ok,bad
print("=== sequence-level proof of the mature-peptide binning ===")
a=test('p200 = p150 + p90','P200',['P150','P90'])
b=test('p110 = C + E2 + E1','P110',['C','E2','E1'])
print("\n=== N-terminal agreement (does each product start where it should?) ===")
for whole,parts in [('P200',['P150','P90']),('P110',['C','E2','E1'])]:
    for gid,d in list(per.items()):
        if whole in d and all(p in d for p in parts):
            w=d[whole]; off=0; line=[]
            for p in parts:
                line.append(f"{p} starts at {off+1} ({d[p][:6]}...)")
                off+=len(d[p])
            print(f"  {gid} {whole} {len(w)} aa: " + '; '.join(line))
            break
