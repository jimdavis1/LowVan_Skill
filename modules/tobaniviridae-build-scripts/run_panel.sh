#!/bin/bash
#  Held-out panel for one module: GTOs -> annotate -> quality score.
set -eu
M=$1
W=$HOME/lowvan_scratch/tobani
cd "$W"
mkdir -p evalrep tbl gto gto_out
if [ ! -d "gto/$M" ]; then
  ./r python3 kit/skill/scripts/make_gto.py --fasta-dir "panel/$M" \
      --metadata "panel/$M.tsv" --out "gto/$M" --jobs 8 | tail -2
fi
./r python3 kit/skill/scripts/run_gto_eval.py --gto-dir "gto/$M" \
    --repo "$W/VA_tob" --out "gto_out/$M" --jobs 4 \
    --report "evalrep/$M.tsv" --tbl-dir "tbl/$M"
