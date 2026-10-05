#  Rebuild a genus module from collections enlarged by project_matpeptides.py.
#  Usage: rebuild_genus.sh <Genus>     (expects fin/<G>/Extracted and patch/<G>.json if any)
source /Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Picornaviridae/env.sh
G=$1; S=$HOME/lowvan_scratch/picorna; SK=$LOWVAN_KIT/skill/scripts
P=/Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Picornaviridae
B=$S/build/$G; rm -rf $B; mkdir -p $B/Alignments/$G; cd $B
echo "######## $(date '+%H:%M') rebuild $G"
rsync -a $S/fin/$G/collections $S/fin/$G/Extracted $B/
rsync -a $P/work/$G/Rep-Contigs $B/
cp $LOWVAN_KIT/modules/$G/${G}_Viral_PSSM.json $B/
rsync -a $S/fin/$G/Alignments/$G/POLY $B/Alignments/$G/      # POLY collection unchanged: keep its build
python3 $S/merge_and_json.py $B $G $( [ -f $S/patch/$G.json ] && echo $S/patch/$G.json )
python3 - $G <<'PY' > features.json
import json,sys,os
G=sys.argv[1]; F=json.load(open("%s_Viral_PSSM.json"%G))[G]["features"]
spec={}
for k,f in sorted(F.items()):
    if k=="POLY": continue
    seqs=[]; cur=0
    for l in open("collections/%s/%s.fasta"%(G,k)):
        if l.startswith(">"):
            if cur: seqs.append(cur)
            cur=0
        else: cur+=len(l.strip())
    if cur: seqs.append(cur)
    n=len(seqs); med=sorted(seqs)[n//2] if n else 0
    #  Thin or very short collections cannot cluster at the documented -mi 0.8.
    #  Senecavirus 2A is 8 sequences of 9 residues (the StopGo NPG peptide): at
    #  -m 2 -mi 0.8 it built nothing and shipped as a declared feature with no
    #  profile. The Oct 2 build had used -m 1 -mi 0.6; restore that floor.
    if n < 10 or med < 20: spec[k]={"anno":f["anno"],"params":{"m":1,"mi":0.6,"mc":0.6}}
    elif n < 20:           spec[k]={"anno":f["anno"],"params":{"m":2}}
    else:                  spec[k]={"anno":f["anno"],"params":{}}
print(json.dumps(spec,indent=1))
PY
cat features.json | python3 -c "import json,sys;print('  features:',{k:v['params'] for k,v in json.load(sys.stdin).items()})"
python3 $SK/build_features.py --workdir . --module $G --features features.json 2>&1 | grep -vE '^\s*$' | tail -16
echo "=== leftover pass"
for K in $(python3 -c "import json;print(' '.join(json.load(open('features.json'))))"); do
  python3 $SK/build_leftover_pssms.py --workdir . --module $G --key $K --write 2>&1 | grep -E "leftover sequence|built|nothing" | sed "s/^/  $K: /"
done
echo "=== install check"
python3 $SK/install_module.py --workdir . --repo $HOME/lowvan_scratch/VA_picorna --only $G --check 2>&1 | tail -6
#  a declared feature with no profile is never called and is silent at runtime
python3 - "$G" <<'PYCHK'
import json,os,sys,glob
G=sys.argv[1]; F=json.load(open("%s_Viral_PSSM.json"%G))[G]["features"]
bad=[k for k in F if "special" not in F[k]
     and not glob.glob("Alignments/%s/%s/pssms/*.pssm"%(G,k))]
print("  FEATURES DECLARED WITH NO PSSM: %s"%(bad if bad else "none"))
if bad: sys.exit("REFUSING: build or undeclare these before shipping")
PYCHK
python3 $SK/json_canon.py --check ${G}_Viral_PSSM.json | tail -1
echo "REBUILD_${G}_COMPLETE"
