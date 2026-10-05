source /Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Picornaviridae/env.sh
S=$HOME/lowvan_scratch/picorna; SK=$LOWVAN_KIT/skill/scripts; PY=/Users/jdavis/anaconda3/bin/python3
for G in PICO_LDR_VP4 PICO_NOLDR_VP4 PICO_LDR_VP0 PICO_NOLDR_VP0; do
  B=$S/build/$G; cd $B
  MISS=$($PY -c "
import json,glob,os
G='$G'; B=os.getcwd()
F=json.load(open(G+'_Viral_PSSM.json'))[G]['features']
print(' '.join(k for k in F if 'special' not in F[k] and not glob.glob('Alignments/%s/%s/pssms/*.pssm'%(G,k))))")
  [ -z "$MISS" ] && { echo "  $G: nothing missing"; continue; }
  echo "######## $G rebuilding: $MISS"
  $PY - "$G" $MISS <<'PYG' > fix.json
import json,sys
G=sys.argv[1]; keys=sys.argv[2:]
F=json.load(open("%s_Viral_PSSM.json"%G))[G]["features"]
spec={}
for k in keys:
    L=[];cur=0
    for l in open("collections/%s/%s.fasta"%(G,k)):
        if l.startswith(">"):
            if cur: L.append(cur)
            cur=0
        else: cur+=len(l.strip())
    if cur: L.append(cur)
    n=len(L); med=sorted(L)[n//2] if n else 0
    #  POLY and other long/diverse collections need a lower identity floor to
    #  cluster at all; the documented -mi 0.8 leaves a polyprotein set in
    #  singletons across a multi-genus module.
    spec[k]={"anno":F[k]["anno"],"params":{"m":2,"mi":0.6,"mc":0.6}}
print(json.dumps(spec,indent=1))
PYG
  python3 $SK/build_features.py --workdir . --module $G --features fix.json 2>&1 | tail -4
  for K in $MISS; do python3 $SK/build_leftover_pssms.py --workdir . --module $G --key $K --write 2>&1 | grep -E "built [1-9]" | sed "s/^/   $K: /" || true; done
  cd $S
done
echo "=== verify and reinstall"
for G in PICO_LDR_VP4 PICO_NOLDR_VP4 PICO_LDR_VP0 PICO_NOLDR_VP0; do
  $PY -c "
import json,glob,os,sys
G='$G'; B=os.path.expanduser('~/lowvan_scratch/picorna/build/'+G)
F=json.load(open(B+'/'+G+'_Viral_PSSM.json'))[G]['features']
bad=[k for k in F if 'special' not in F[k] and not glob.glob(B+'/Alignments/%s/%s/pssms/*.pssm'%(G,k))]
print('  %-16s FEATURES WITH NO PSSM: %s'%(G, bad if bad else 'none'))"
  python3 $SK/install_module.py --workdir $S/build/$G --repo $HOME/lowvan_scratch/VA_picorna --only $G 2>&1 | grep -E "declared|problem"
  cp $S/build/$G/${G}_Viral_PSSM.json $S/
done
echo FIX_GROUPS_COMPLETE
