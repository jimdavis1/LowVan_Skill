#  Grouped modules, second attempt at the member-floor problem.
#
#  A first attempt lowered -mi to 0.6 on the main pass. NEVER DO THAT: a lower
#  identity MERGES genera into one cluster instead of giving each its own,
#  and it left VP0 and VPG with no alignments at all. Profile count went
#  28 -> 27, i.e. nothing gained.
#
#  Correct mechanism: keep the documented main pass (-m 5 -mi 0.8), which
#  clusters the dominant genus properly, then run the LEFTOVER pass at a
#  member floor of 1. The leftover pool is by definition the sequences that
#  reached no cluster -- in these modules that is exactly the small genera --
#  and the skill's own reasoning is that those are the divergent members a
#  profile set most needs to cover. The identity floor is NOT lowered.
source /Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Picornaviridae/env.sh
S=$HOME/lowvan_scratch/picorna; SK=$LOWVAN_KIT/skill/scripts; PY=/Users/jdavis/anaconda3/bin/python3
for G in PICO_NOLDR_VP0 PICO_NOLDR_VP4 PICO_LDR_VP0 PICO_LDR_VP4; do
  B=$S/build/$G; cd $B
  echo "######## $(date +%H:%M) $G"
  rm -rf Alignments/$G; mkdir -p Alignments/$G
  $PY - "$G" <<'PYG' > features.json
import json,sys
G=sys.argv[1]; F=json.load(open("%s_Viral_PSSM.json"%G))[G]["features"]
spec={}
for k,v in sorted(F.items()):
    n=sum(1 for l in open("collections/%s/%s.fasta"%(G,k)) if l.startswith(">"))
    spec[k]={"anno":v["anno"],"params":({} if n>=20 else {"m":2} if n>=6 else {"m":1})}
print(json.dumps(spec,indent=1))
PYG
  python3 $SK/build_features.py --workdir . --module $G --features features.json 2>&1 | tail -2
  echo "   --- leftover pass at a member floor of 1"
  for K in $($PY -c "import json;print(' '.join(json.load(open('features.json'))))"); do
    python3 $SK/build_leftover_pssms.py --workdir . --module $G --key $K --min-seqs 1 --min-len-frac 0.4 --write 2>&1 \
      | grep -E "leftover sequence|built|nothing" | tr '\n' ' ' | sed "s/^/   $K: /"; echo
  done
  python3 $SK/install_module.py --workdir $B --repo $HOME/lowvan_scratch/VA_picorna --only $G 2>&1 | grep -E "declared|FATAL|REFUS"
  cd $S
done
echo REBUILD2_COMPLETE
