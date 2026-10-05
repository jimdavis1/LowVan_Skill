#  Build the four grouped Picornaviridae modules.
source /Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Picornaviridae/env.sh
set -e
S=$HOME/lowvan_scratch/picorna; SK=$LOWVAN_KIT/skill/scripts; PY=/Users/jdavis/anaconda3/bin/python3
cd $S
echo "=== $(date +%H:%M) fetch contigs"
for G in PICO_LDR_VP4 PICO_NOLDR_VP4 PICO_LDR_VP0 PICO_NOLDR_VP0; do
  [ -d Contigs.$G ] || $PY fetch_contigs.py --meta meta/$G.tsv --out Contigs.$G --min 5000 --jobs 4 2>&1 | tail -1 | sed "s/^/  $G: /"
done
echo "=== $(date +%H:%M) declare features and windows"
$PY $S/rules/make_group_json.py
echo "=== $(date +%H:%M) build profiles"
for G in PICO_LDR_VP4 PICO_NOLDR_VP4 PICO_LDR_VP0 PICO_NOLDR_VP0; do
  B=$S/build/$G; mkdir -p $B/Alignments/$G
  rsync -a collections/$G $B/collections/
  cp $S/${G}_Viral_PSSM.json $B/
  cd $B
  $PY - "$G" <<'PYG' > features.json
import json,sys,os
G=sys.argv[1]; F=json.load(open("%s_Viral_PSSM.json"%G))[G]["features"]
spec={}
for k,f in sorted(F.items()):
    if k=="POLY": continue
    p="collections/%s/%s.fasta"%(G,k)
    L=[];cur=0
    for l in open(p):
        if l.startswith(">"):
            if cur: L.append(cur)
            cur=0
        else: cur+=len(l.strip())
    if cur: L.append(cur)
    n=len(L); med=sorted(L)[n//2] if n else 0
    if n<10 or med<20: spec[k]={"anno":f["anno"],"params":{"m":1,"mi":0.6,"mc":0.6}}
    elif n<20:         spec[k]={"anno":f["anno"],"params":{"m":2}}
    else:              spec[k]={"anno":f["anno"],"params":{}}
print(json.dumps(spec,indent=1))
PYG
  echo "######## $G"
  python3 $SK/build_features.py --workdir . --module $G --features features.json 2>&1 | tail -3
  for K in $($PY -c "import json;print(' '.join(json.load(open('features.json'))))"); do
    python3 $SK/build_leftover_pssms.py --workdir . --module $G --key $K --write 2>&1 | grep -E "built [1-9]" | sed "s/^/   $K: /" || true
  done
  cd $S
done
echo "=== $(date +%H:%M) rep contigs"
for G in PICO_LDR_VP4 PICO_NOLDR_VP4 PICO_LDR_VP0 PICO_NOLDR_VP0; do
  B=$S/build/$G
  python3 $SK/build_rep_contigs.py --contigs Contigs.$G --meta meta/$G.tsv --module $G \
     --budget 25 --min-len 5000 --out $B/Rep-Contigs --threads 4 2>&1 | tail -3 | sed "s/^/  $G: /"
  $PY - "$G" <<'PYC'
import json,sys,os
G=sys.argv[1]; S=os.path.expanduser("~/lowvan_scratch/picorna"); B=f"{S}/build/{G}"
cg=json.load(open(f"{B}/Rep-Contigs/close_genomes.json"))
p=f"{B}/{G}_Viral_PSSM.json"; d=json.load(open(p)); d[G]["close_genomes"]=cg
def norm(x):
    if isinstance(x,dict): return {k:norm(v) for k,v in x.items()}
    if isinstance(x,list): return [norm(v) for v in x]
    if isinstance(x,float) and x==int(x): return int(x)
    return x
open(p,"w").write(json.dumps(norm(d),indent=3,separators=(","," : "),sort_keys=True,ensure_ascii=False)+"\n")
print("   close_genomes: %d"%len(cg))
PYC
done
echo "=== $(date +%H:%M) install"
for G in PICO_LDR_VP4 PICO_NOLDR_VP4 PICO_LDR_VP0 PICO_NOLDR_VP0; do
  python3 $SK/install_module.py --workdir $S/build/$G --repo $HOME/lowvan_scratch/VA_picorna --only $G 2>&1 | grep -E "declared|PSSMs ->|rep contigs|problem"
  cp $S/build/$G/${G}_Viral_PSSM.json $S/
done
echo BUILD_GROUPS_COMPLETE
