#!/bin/bash
#  Render the coverage audit. gen_coverage_audit.py resolves a figure block's
#  `from` relative to the CURRENT directory, so it must run in pages/.
set -eu
W=$HOME/lowvan_scratch/tobani
export LOWVAN_KIT=$W/kit
export PATH="/Applications/BV-BRC.app/runtime/bin:/Users/jdavis/anaconda3/envs/lowvan/bin:$PATH"
cd "$W/pages"
python3 "$W/kit/skill/scripts/gen_coverage_audit.py" \
   --facts tobaniviridae_audit.json --out tobaniviridae-coverage-audit.html
python3 "$W/kit/skill/scripts/check_artifact_html.py" \
   tobaniviridae-coverage-audit.html tobaniviridae-pssm-registry.html \
   tobaniviridae-string-collapse.html tobaniviridae-vocabulary-saturation.html
