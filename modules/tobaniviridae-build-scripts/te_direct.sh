#!/bin/bash
#  Run the transcript-edit caller on one genome's annotated GTO and show its
#  own diagnosis on stderr.
set -eu
M=$1; G=$2
W=$HOME/lowvan_scratch/tobani
cd "$W"
export LOWVAN_DATA_DIR=$W/VA_tob
export PERL5LIB=$W/kit/vendor/bv-brc/lib:/Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/LowVan_Skill_Development/_infrastructure/Viral_Annotation:${PERL5LIB:-}
export PATH=/Applications/BV-BRC.app/runtime/bin:$PATH
perl "$W/VA_tob/get_transcript_edited_features.pl" \
   -i "gto_out/$M/$G.ann.gto" -o /dev/null 2>&1 | tail -25
