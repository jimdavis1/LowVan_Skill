#!/bin/bash
set -u
# LOWVAN_KIT: the clone of this repository. LOWVAN_WORK: your build directory.
KIT=${LOWVAN_KIT:?set LOWVAN_KIT to the LowVan_Skill clone}
WORK=${LOWVAN_WORK:-$PWD}
D="-m 5 -f 0.33 -n 0.75 -c 0 -fd 0.20 -e 3 -efo 0.15 -mi 0.8 -mc 0.8 -p_nterm 65 -n_nterm 3 -nterm_eval_length 10 -nterm-mmseq-id 0.7"
k=$1; anno=$2
out=$WORK/Alignments/Togaviridae/$k; rm -rf "$out"; mkdir -p "$out"; cd "$out" || exit 1
cat "$WORK/collections/Togaviridae/$k.build.fasta" \
 | perl "$KIT/build/fasta-cluster-pssm-2.pl" -a "$anno" -p "Togaviridae.$k" -tmp build -x $D > run.stdout 2> run.stderr
[ -d build ] && { find build -maxdepth 1 -mindepth 1 ! -name 'build.*' -exec mv {} . \; 2>/dev/null; rm -rf build; }
{ echo "feature      $k"; echo "annotation   $anno"; echo "defaults     $D";
  echo "departures   -x  (keep non-standard residues). Input is $k.build.fasta,";
  echo "             which admits X ONLY in the recovered sequences. Required";
  echo "             because P123 carries the readthrough opal as X and";
  echo "             FCP_Main_Utils.pm:244 drops /(B|J|X|Z)/i without -x.";
  echo "source       Extracted/Togaviridae.$k.faa (extract_combos.py) -- the";
  echo "             BV-BRC export has no usable collection for this feature.";
  echo "built        $(date -u +%Y-%m-%dT%H:%M:%SZ)"; } > BUILD_PARAMS
echo "  $k: clusters=$(ls clusters/*.fasta 2>/dev/null|wc -l|tr -d ' ') alis=$(ls corrected_alis/*.fa 2>/dev/null|wc -l|tr -d ' ') pssms=$(ls pssms/*.pssm 2>/dev/null|wc -l|tr -d ' ') leftovers=$(grep -c '^>' Leftover_Seqs.aa 2>/dev/null||echo 0)"
