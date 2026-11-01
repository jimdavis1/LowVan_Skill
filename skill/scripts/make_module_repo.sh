#!/bin/bash
#  Build a repo whose Viral-Rep-Contigs holds ONE module's references.
#
#  Routing runs one blastn per reference file. VA_parvo_route cut that from
#  1,271 to 176 and that is what made the family-wide passes survivable; this
#  cuts it again, to the 25 or so a single module ships, for a re-run of that
#  module's OWN panel.
#
#  It is exact for that panel, and only because it was measured. The caller
#  keeps the max-bit reference, so restricting to a subset that contains the
#  full-set winner returns the same winner. logs/route_screen.tsv holds the
#  winner for every sampled genome against all 1,271 references: every one of
#  Protoparvovirus's 300 wins on a Protoparvovirus.*.dna file, 24 distinct
#  ones. Check that before using this for another module -- a panel that wins
#  on a neighbour's contig will route differently here, and silently.
#
#  A panel does not always win on its own module's contigs:
#  DensovirinaeAmbi's 184 split 152 / 25 DensovirinaeMono / 2
#  Hamaparvovirinae / 2 Protoparvovirus / 2 Gammacoronavirus / 1
#  Trivirinae. Pass every module the screen names, or the ones you
#  leave out will be re-routed somewhere else without a word.
#
#  Usage: bash make_module_repo.sh <Module>[,<Module>...] [dest]
set -euo pipefail
M=$1
MODS=$(echo "$M" | tr ',' ' ')
SRC=$HOME/lowvan_scratch/VA_parvo_route
DST=${2:-$HOME/lowvan_scratch/VA_only_$(echo "$M" | tr ',' '_')}
rm -rf "$DST"; mkdir -p "$DST/Viral-Rep-Contigs"
for f in "$SRC"/* "$SRC"/.[!.]*; do
  b=$(basename "$f")
  [ "$b" = "Viral-Rep-Contigs" ] && continue
  [ -e "$DST/$b" ] && continue
  ln -s "$(cd "$(dirname "$f")" && readlink "$b" 2>/dev/null || echo "$f")" "$DST/$b" 2>/dev/null || ln -s "$f" "$DST/$b"
done
#  A bare module name takes all of that module's contigs; a name containing
#  a dot is taken as one reference FILE, which is how a panel that wins on
#  two files of a neighbouring module avoids dragging in that module's other
#  twenty-three.
n=0
for m in $MODS; do
  case "$m" in
    *.dna)
      f="$SRC/Viral-Rep-Contigs/$m"
      [ -e "$f" ] && { ln -s "$(readlink "$f" 2>/dev/null || echo "$f")" "$DST/Viral-Rep-Contigs/$m"; n=$((n+1)); }
      continue;;
  esac
  for f in "$SRC"/Viral-Rep-Contigs/"$m".*.dna; do
    [ -e "$f" ] || continue
    ln -s "$(readlink "$f" 2>/dev/null || echo "$f")" "$DST/Viral-Rep-Contigs/$(basename "$f")"
    n=$((n+1))
  done
done
echo "$DST  --  $n reference contig(s) for $M (source had $(ls "$SRC"/Viral-Rep-Contigs | wc -l | tr -d ' '))"
