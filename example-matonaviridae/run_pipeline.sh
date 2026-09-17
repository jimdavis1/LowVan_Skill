#!/bin/bash
# Build clusters / alignments / PSSMs for one or more Matonaviridae features.
#
# Layout produced (what install_module.py and rebuild_pssms.py expect):
#   Alignments/Matonaviridae/<FEAT>/{clusters,alis,corrected_alis,pssms}/
#
# Parameters are the documented defaults throughout. Any departure is recorded
# in Alignments/Matonaviridae/<FEAT>/BUILD_PARAMS.
set -u
KIT=${LOWVAN_KIT:?set LOWVAN_KIT to the LowVan_Skill clone}
WORK=${LOWVAN_WORK:-$PWD}
DEFAULTS="-m 5 -f 0.33 -n 0.75 -c 0 -fd 0.20 -e 3 -efo 0.15 -mi 0.8 -mc 0.8 \
 -p_nterm 65 -n_nterm 3 -nterm_eval_length 10 -nterm-mmseq-id 0.7"

feat=$1; anno=$2; shift 2
extra="$*"
out=$WORK/Alignments/Matonaviridae/$feat
rm -rf "$out"; mkdir -p "$out"
cd "$out" || exit 1
echo "=== $feat : $anno"
# -tmp names the working subdirectory the pipeline creates; without it the name
# is random, and the outputs land somewhere different on every run.
cat "$WORK/collections/Matonaviridae/$feat.fasta" \
  | perl "$KIT/build/fasta-cluster-pssm-2.pl" -a "$anno" -p "Matonaviridae.$feat" \
        -tmp build $DEFAULTS $extra > run.stdout 2> run.stderr
rc=$?
# lift the run's outputs to <FEAT>/ , which is the layout install_module.py and
# rebuild_pssms.py read
if [ -d build ]; then
  find build -maxdepth 1 -mindepth 1 ! -name 'build.*' -exec mv {} . \; 2>/dev/null
  rm -rf build
fi
{ echo "feature      $feat"; echo "annotation   $anno";
  echo "defaults     $DEFAULTS";
  echo "departures   ${extra:-none}";
  echo "built        $(date -u +%Y-%m-%dT%H:%M:%SZ)"; } > BUILD_PARAMS
nclu=$(ls clusters/*.fasta 2>/dev/null | wc -l | tr -d ' ')
nali=$(ls corrected_alis/*.fa 2>/dev/null | wc -l | tr -d ' ')
npssm=$(ls pssms/*.pssm 2>/dev/null | wc -l | tr -d ' ')
nleft=$(grep -c '^>' Leftover_Seqs.aa 2>/dev/null || echo 0)
echo "    rc=$rc clusters=$nclu corrected_alis=$nali pssms=$npssm leftovers=$nleft"
