#!/bin/bash
#  Annotate Berne virus alone and show every M profile's score, to size the
#  false 1,070 aa Membrane glycoprotein call against the genuine ones.
set -eu
W=$HOME/lowvan_scratch/tobani
export LOWVAN_DATA_DIR=$W/VA_tob
export PERL5LIB=$W/kit/vendor/bv-brc/lib:/Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/LowVan_Skill_Development/_infrastructure/Viral_Annotation:${PERL5LIB:-}
export PATH=/Applications/BV-BRC.app/runtime/bin:$PATH
mkdir -p "$W/probe" && cd "$W/probe"
cp "$W/Contigs/11156.6.fna" .
perl "$W/VA_tob/annotate_by_viral_pssm.pl" -i 11156.6.fna -p berne \
     -c "$W/VA_tob/Viral-Rep-Contigs" -max 500000 -min 200 2> berne.err > /dev/null || true
grep -A12 'Torovirinae	M	Membrane' berne.err | head -16
echo '--- called features ---'
awk -F'\t' '{print $7, $8, $9, $10}' berne.feature.tbl
