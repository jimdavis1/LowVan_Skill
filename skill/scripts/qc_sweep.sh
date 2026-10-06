#  The QC checks from the skill's "what done looks like", across all 14
#  Picornaviridae modules. Read-only: reports, changes nothing.
source /Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Picornaviridae/env.sh
S=$HOME/lowvan_scratch/picorna; SK=$LOWVAN_KIT/skill/scripts; K=$LOWVAN_KIT
VA=$HOME/lowvan_scratch/VA_picorna
M="Avihepatovirus Teschovirus Cardiovirus Sapelovirus Senecavirus Kobuvirus Parechovirus Hepatovirus Aphthovirus Enterovirus PICO_LDR_VP4 PICO_NOLDR_VP4 PICO_LDR_VP0 PICO_NOLDR_VP0"
while pgrep -f "eval_groups.sh|nogenus.sh" >/dev/null 2>&1; do sleep 60; done
cd $S
echo "################ 1. install --check (every module)"
for G in $M; do
  B=$S/build/$G; [ -d $B ] || B=$S/fin/$G
  python3 $SK/install_module.py --workdir $B --repo $VA --only $G --check 2>&1 | grep -E "declared|FATAL|problem" | sed "s/^/  /"
done
echo; echo "################ 2. json_canon --check"
for G in $M; do python3 $SK/json_canon.py --check $K/modules/$G/${G}_Viral_PSSM.json 2>&1 | tail -1; done
echo; echo "################ 3. check_annotations (vocabulary)"
for G in $M; do echo "--- $G"; python3 $SK/check_annotations.py --json $K/modules/$G/${G}_Viral_PSSM.json 2>&1 | tail -4; done
echo; echo "################ 4. check_dlits (citations)"
for G in $M; do echo "--- $G"; python3 $SK/check_dlits.py --json $K/modules/$G/${G}_Viral_PSSM.json 2>&1 | tail -4; done
echo; echo "################ 5. qc_cross_feature (mislabelled clusters)"
for G in $M; do
  B=$S/build/$G; [ -d $B ] || B=$S/fin/$G
  echo "--- $G"; (cd $B && python3 $SK/qc_cross_feature.py --workdir . 2>&1 | tail -5)
done
echo; echo "################ 6. rebuild_pssms audit (stale profiles)"
for G in $M; do
  B=$S/build/$G; [ -d $B ] || B=$S/fin/$G
  echo "--- $G"; (cd $B && python3 $SK/rebuild_pssms.py 2>&1 | tail -3)
done
echo QC_SWEEP_COMPLETE
