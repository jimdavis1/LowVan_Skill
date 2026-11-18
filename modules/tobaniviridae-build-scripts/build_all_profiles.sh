#!/bin/bash
#  Profiles for one module: the clustered pass, then the leftover pass at the
#  same -mi, then one profile per still-uncovered sequence.
set -eu
M=$1
W=$HOME/lowvan_scratch/tobani
cd "$W"
./r python3 kit/skill/scripts/build_features.py --workdir "$M" --module "$M" \
   --features "$M/features.json" | tail -20
KEYS=$(python3 -c "
import json,sys
f=json.load(open('$M/features.json'))
print(' '.join(k for k,v in f.items() if not v.get('special')))")
for k in $KEYS; do
  ./r python3 kit/skill/scripts/build_leftover_pssms.py --workdir "$M" \
      --module "$M" --key "$k" --write 2>/dev/null | tail -1 || true
  ./r python3 scripts/build_singleton_pssms.py "$M" "$k" --write | tail -1
done
echo -n "$M profiles: "; find "$M/Alignments" -name '*.pssm' | wc -l
