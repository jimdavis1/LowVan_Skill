#  Ship a rebuilt genus: runtime (Box), repo module, Box work tree. Usage: ship_genus.sh <Genus>
source /Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Picornaviridae/env.sh
G=$1; S=$HOME/lowvan_scratch/picorna; SK=$LOWVAN_KIT/skill/scripts
P=/Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Picornaviridae
B=$S/build/$G; M=$LOWVAN_KIT/modules/$G; W=$P/work/$G; STAMP=v1_20261002
set -e
echo "=== runtime"
python3 $SK/install_module.py --workdir $B --repo $LOWVAN_DATA_DIR --only $G 2>&1 | grep -E "PSSMs ->|alignments ->|rep contigs|Viral_PSSM.json:"
echo "=== repo module"
cp $B/${G}_Viral_PSSM.json $M/
[ -f $S/notes/$G.md ] && cp $S/notes/$G.md $M/BUILD_NOTES.md
for F in $(ls $B/Alignments/$G); do
  D=$M/PSSM-Alignments/$F; mkdir -p $D; A=$B/Alignments/$G/$F
  for d in corrected_alis reclustered_alis alis clusters truncation_qc; do [ -d $A/$d ] && rsync -a --delete $A/$d/ $D/$d/; done
  for f in Curation_Report Truncation_Report FLAGGED_NTERM BUILD_PARAMS Leftover_Seqs.aa nterm_evaluation.log; do [ -s $A/$f ] && cp $A/$f $D/; done
done
#  feature dirs the rebuild no longer has are removed from the repo (git keeps them)
for D in $(ls $M/PSSM-Alignments); do [ -d $B/Alignments/$G/$D ] || rm -rf $M/PSSM-Alignments/$D; done
echo "  $(ls $M/PSSM-Alignments | wc -l) feature dirs"
echo "=== Box work tree (previous build kept as *.$STAMP)"
[ -d $W/Alignments.$STAMP ] || mv $W/Alignments $W/Alignments.$STAMP
[ -d $W/collections.$STAMP ] || mv $W/collections $W/collections.$STAMP
#  Extracted/ exists only for genera that were projected
rsync -a --exclude 'tmp/' --exclude '.L-renamed-to-*' $B/Alignments $B/collections $W/
[ -d $B/Extracted ] && rsync -a $B/Extracted $W/
cp $B/${G}_Viral_PSSM.json $W/
[ -f $B/features.json ] && cp $B/features.json $W/
[ -f $S/notes/$G.md ] && cp $S/notes/$G.md $W/BUILD_NOTES.md
TE=$HOME/lowvan_scratch/VA_picorna/Transcript-Editing/$G
if [ -d "$TE" ]; then
  mkdir -p $LOWVAN_DATA_DIR/Transcript-Editing/$G $M/Transcript-Editing $W/Transcript-Editing
  rsync -a $TE/ $LOWVAN_DATA_DIR/Transcript-Editing/$G/
  rsync -a $TE/ $M/Transcript-Editing/
  rsync -a $TE/ $W/Transcript-Editing/
  echo "  transcript-edit refs: $(grep -c '>' $TE/*.fasta) sequences -> runtime, repo, work tree"
fi
echo "SHIP_${G}_COMPLETE"
