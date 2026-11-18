#!/bin/bash
#  One genome map per module, from a panel genome where every declared feature
#  is called exactly once. The map is what shows a reader what the module does
#  to a genome: where the features sit, which are essential, and which are
#  cleavage products of which precursor.
set -eu
W=$HOME/lowvan_scratch/tobani
cd "$W"
mkdir -p pages
make_one () {
  M=$1; G=$2; CAP=$3
  ./r python3 kit/skill/scripts/gen_genome_map.py \
      --gto "gto_out/$M/$G.ann.gto" \
      --module-json "$M/${M}_Viral_PSSM.json" --module "$M" \
      --out "pages/tobaniviridae-genome-map-$M.html" --caption "$CAP"
}
make_one Torovirinae 74501.274 'Bovine torovirus 12151, the type genus of the family.'
make_one Piscanivirinae 889873.3 'Fathead minnow nidovirus, one of the two bafiniviruses whose genome was described.'
make_one Remotovirinae 1986112.9 'Bovine nidovirus 1 NSWL2; this genus has one genome paper and the smallest genome in the family.'
make_one Serpentovirinae 2016400.8 'Morelia viridis nidovirus BH171/14-7, a green tree python serpentovirus.'
ls -la pages/tobaniviridae-genome-map-*.html
