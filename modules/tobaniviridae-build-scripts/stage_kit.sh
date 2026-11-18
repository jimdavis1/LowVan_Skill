#!/bin/bash
#  Copy the four modules into the kit layout: JSON, rep contigs, the alignments
#  behind every shipped profile, the transcript-edit references, the audits and
#  the per-module build inputs.
set -eu
W=$HOME/lowvan_scratch/tobani
K=${1:?usage: stage_kit.sh <kit worktree>}
for M in Torovirinae Piscanivirinae Remotovirinae Serpentovirinae; do
  D="$K/modules/$M"
  mkdir -p "$D/Rep-Contigs" "$D/PSSM-Alignments"
  cp "$W/$M/${M}_Viral_PSSM.json" "$D/"
  cp "$W/$M/features.json" "$D/features_all.json"
  cp "$W/$M/length_refs.json" "$D/"
  cp "$W/$M/build_collections.py" "$D/"
  cp "$W/$M/LENGTH_AUDIT.md" "$W/$M/BIN_AUDIT.md" "$D/"
  [ -f "$W/$M/length_accept.tsv" ] && cp "$W/$M/length_accept.tsv" "$D/" || true
  cp "$W/$M"/Rep-Contigs/*.dna "$D/Rep-Contigs/"
  cp "$W/$M"/Rep-Contigs/close_genomes.json "$D/Rep-Contigs/"
  mkdir -p "$K/Transcript-Editing/$M"
  cp "$W/$M/Transcript-Editing/$M/ORF1AB.fasta" "$K/Transcript-Editing/$M/"
  cp "$W/$M/Transcript-Editing/$M/ORF1AB.build.tsv" "$D/"
  #  the alignments behind the shipped profiles, which are the provenance
  for F in "$W/$M"/Alignments/"$M"/*/; do
    f=$(basename "$F")
    [ -d "$F/pssms" ] || continue
    mkdir -p "$D/PSSM-Alignments/$f"
    for sub in corrected_alis reclustered_alis; do
      [ -d "$F/$sub" ] && cp -R "$F/$sub" "$D/PSSM-Alignments/$f/" || true
    done
    [ -f "$F/BUILD_PARAMS" ] && cp "$F/BUILD_PARAMS" "$D/PSSM-Alignments/$f/" || true
  done
  #  the tracking files the engine wrote
  mkdir -p "$D/collections"
  for t in synonyms.tsv UNASSIGNED_TRACKING.tsv NOT_MODELLED.tsv EXPLICIT_WINDOWS.tsv \
           RUNTS_DROPPED.tsv OVERLONG_DROPPED.tsv AMBIGUOUS_DROPPED.tsv; do
    [ -f "$W/$M/collections/$t" ] && cp "$W/$M/collections/$t" "$D/collections/" || true
  done
  n=$(find "$D/PSSM-Alignments" -name '*.fa' | wc -l | tr -d ' ')
  p=$(find "$W/$M/Alignments" -name '*.pssm' | wc -l | tr -d ' ')
  echo "$M: $n alignment(s) staged for $p installed profile(s), $(ls "$D/Rep-Contigs"/*.dna | wc -l | tr -d ' ') rep contig(s)"
done
