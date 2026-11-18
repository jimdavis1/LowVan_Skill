#!/bin/bash
#  The 26 exemplars the GTO pipeline could not annotate failed at
#  rast-create-genome with "403 Forbidden" from the BV-BRC ID server -- they
#  were never offered to a module, so calling them UNROUTED (as the coverage
#  report does) reads as a routing failure when it is an infrastructure one.
#  annotate_by_viral_pssm.pl needs no ID server, so routing can be measured
#  for them directly. It cannot produce special features or quality scores.
set -eu
W=$HOME/lowvan_scratch/tobani
export LOWVAN_DATA_DIR=$W/VA_tob
export PERL5LIB=$W/kit/vendor/bv-brc/lib:/Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/LowVan_Skill_Development/_infrastructure/Viral_Annotation:${PERL5LIB:-}
export PATH=/Applications/BV-BRC.app/runtime/bin:$PATH
O=$W/cov/flat
mkdir -p "$O"
cd "$O"
for f in "$W"/cov/work/ex_retry3/*; do
  g=$(basename "$f"); g=${g%.fna}; g=${g%.fasta}
  [ -f "$g.feature.tbl" ] && continue
  perl "$W/VA_tob/annotate_by_viral_pssm.pl" -i "$f" -p "$g" \
       -c "$W/VA_tob/Viral-Rep-Contigs" -max 500000 -min 200 \
       > /dev/null 2> "$g.err" || true
done
n=0; r=0
for f in "$W"/cov/work/ex_retry3/*; do
  g=$(basename "$f"); g=${g%.fna}; g=${g%.fasta}
  n=$((n+1))
  if [ -s "$g.feature.tbl" ]; then r=$((r+1)); fi
done
echo "flat annotator: $r of $n ID-server-failed exemplars routed and annotated"
