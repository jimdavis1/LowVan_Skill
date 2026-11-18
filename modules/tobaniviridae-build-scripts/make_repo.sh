#!/bin/bash
#  A runtime repo holding the Tobaniviridae family and nothing else: the four
#  modules of this family are exactly what SKILL.md step 11 means by "the
#  family's rep contigs, nothing wider or narrower".
set -eu
W=$HOME/lowvan_scratch/tobani
R=$W/VA_tob
rm -rf "$R"
mkdir -p "$R/Viral-PSSMs" "$R/Viral-Rep-Contigs" "$R/PSSM-Alignments" \
         "$R/Transcript-Editing" "$R/Splice-Variants" "$R/Other_Scripts"
cp "$W"/kit/annotate/*.pl "$W"/kit/annotate/*.pm "$R"/
cp "$W"/kit/build/*.pl "$W"/kit/build/*.pm "$R"/ 2>/dev/null || true
cp -R "$W"/kit/evaluate/. "$R/Other_Scripts/" 2>/dev/null || true
echo '{}' > "$R/Viral_PSSM.json"
ls "$R"
