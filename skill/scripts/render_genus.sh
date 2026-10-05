#  Render the four pages for a genus. Usage: render_genus.sh <Genus> "<date>"
#  Needs reports/generators/<slug>_notes.json and <slug>_prose.json (hand-written).
source /Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Picornaviridae/env.sh
G=$1; DATE=$2; S=$HOME/lowvan_scratch/picorna; SK=$LOWVAN_KIT/skill/scripts
R=$LOWVAN_KIT/reports; s=$(echo $G | tr 'A-Z' 'a-z'); F=$S/fin/$G
cd $F
python3 $SK/gen_registry.py --registry registry.json --taxon $G --notes $R/generators/${s}_notes.json \
   --date "$DATE" --out $R/$s-pssm-registry.html | tail -1
python3 $SK/gen_collapse.py --taxon $G --date "$DATE" --out $R/$s-string-collapse.html | tail -1
python3 $SK/gen_rarefaction.py --rarefaction rarefaction.json --taxon $G --date "$DATE" --out $R/$s-vocabulary-saturation.html | tail -1
#  gen_facts failing silently once left a stale audit JSON in place and the
#  page was rendered from thirteen-feature numbers after two more were added
if ! python3 $S/gen_facts.py $G $R/generators/${s}_prose.json $R/generators/${s}_audit.json > $F/facts_check.json; then
  echo "  gen_facts FAILED for $G -- not rendering the audit from a stale facts file"; exit 1
fi
cd $R && python3 $SK/gen_coverage_audit.py --facts $R/generators/${s}_audit.json --out $R/$s-coverage-audit.html | tail -1
ls -la $R/$s-*.html
