#!/bin/bash
#  Held-out panel evaluation for the four Adenoviridae modules.
#
#  ROUTING SCOPE: VA_adeno carries 150 rep-contig files across 38 modules, of
#  which 72 are Adenoviridae. That is the right breadth -- a genome that
#  misroutes needs its neighbours present to be visible, and a module-only
#  repo absorbs the misroute silently (curator, 2 Nov, on Parvoviridae).
#  150 reference files is ~300 process launches per genome, well under the
#  ~2,600 that exhausted IPC vouchers and panicked the kernel at 1,271 refs.
#
#  Specials are NOT skipped: a run that omits them understates the module.
source /Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Retroviridae/env.sh
P=$HOME/lowvan_scratch/adeno; S=$LOWVAN_KIT/skill/scripts
REPO=$HOME/lowvan_scratch/VA_adeno
export LOWVAN_DATA_DIR=$REPO
JOBS=${JOBS:-4}
for m in Barthadenovirus Siadenovirus Aviadenovirus Mastadenovirus; do
  want=$(ls $P/panel/$m.gto 2>/dev/null | wc -l | tr -d ' ')
  echo "=== $m  $want genome(s)  $(date +%H:%M:%S)"
  python3 $S/run_gto_eval.py --gto-dir $P/panel/$m.gto --repo $REPO \
          --out $P/panel/$m.anno --report $P/panel/$m.report.tsv \
          --tbl-dir $P/panel/$m.tbl --jobs $JOBS \
          > $P/logs/$m.panel.log 2>&1
  echo "    exit=$? tbl=$(ls $P/panel/$m.tbl 2>/dev/null | wc -l | tr -d ' ') of $want"
done
echo "PANEL EVAL DONE $(date +%H:%M:%S)"
