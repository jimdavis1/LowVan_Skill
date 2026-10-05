#  Per-genus measurements for the four artifacts. Usage: finish_genus.sh <Genus>
source /Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Picornaviridae/env.sh
G=$1; S=$HOME/lowvan_scratch/picorna; SK=$LOWVAN_KIT/skill/scripts
P=/Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Picornaviridae
F=$S/fin/$G; mkdir -p $F; cd $F
echo "######## $(date '+%H:%M') $G finish"
#  a rebuilt module (build/<G>) supersedes the shipped one
if [ -d $S/build/$G ]; then SRC=$S/build/$G; J=$SRC/${G}_Viral_PSSM.json; else SRC=$P/work/$G; J=$LOWVAN_KIT/modules/$G/${G}_Viral_PSSM.json; fi
rm -rf $F/Alignments $F/collections   # a stale tree here once added 3 dead PSSMs to the registry
rsync -a --exclude 'tmp/' --exclude '.L-renamed-to-*' $SRC/Alignments $SRC/collections $F/
cp $J $F/
#  L.fasta is the pre-rename copy of Ldr.fasta; it is not a declared feature
[ -f $F/collections/$G/Ldr.fasta ] && rm -f $F/collections/$G/L.fasta
#  genus-only dump, so the collapse denominator is this genus and not the family
python3 $S/subset_dump.py --src $S --prefix Picorna --out $F --name G --genera $G 2>&1 | tail -2
echo "=== registry"
python3 $SK/collect_registry.py --workdir . --json ${G}_Viral_PSSM.json --out registry.json --threads 4 2>&1 | tail -4
echo "=== synonyms + string collapse"
python3 $SK/synmap_from_collections.py --workdir . --module $G --dump G 2>&1 | tail -4
cp collections/synonyms.tsv . 2>/dev/null
python3 $SK/collect_synmap.py . 2>&1 | tail -4
echo "=== source products for the coverage exemplars"
#  rarefaction needs enough genomes: use the coverage exemplars, or for a genus with
#  fewer than 100 of them, annotate a random 200 of its genomes (all, if fewer)
ANN=$S/coverage/wk.$G/ann
if [ $(ls $ANN/*.feature.tbl | wc -l) -lt 100 ]; then
  mkdir -p rare_contigs; python3 -c "import glob,random,shutil;f=sorted(glob.glob('$S/Contigs.$G/*.fna'));random.seed(1);[shutil.copy(x,'rare_contigs/') for x in random.sample(f,min(200,len(f)))]"
  python3 $SK/evaluate_coverage.py --contigs rare_contigs --repo $HOME/lowvan_scratch/VA_picorna \
     --json ${G}_Viral_PSSM.json --work rare --pident 1.0 --cov 1.0 --threads 6 --out rare.json 2>&1 | grep -E "exemplar|routed"
  ANN=$F/rare/ann
fi
ls $ANN/*.feature.tbl | xargs -n1 basename | sed 's/\.feature\.tbl$//' > exemplars.txt
/Users/jdavis/anaconda3/bin/python3 $S/fetch_products.py exemplars.txt bvbrc_products.tsv
echo "=== rarefaction"
python3 $SK/annotation_rarefaction.py --new $ANN --old bvbrc_products.tsv \
   --out rarefaction.json --replicates 100 2>&1 | tail -16
mkdir -p $P/work/$G/artifacts && cp registry.json agg.json unbinned.json typos.json synonyms.tsv \
   rarefaction.json bvbrc_products.tsv $P/work/$G/artifacts/ 2>/dev/null
echo "FINISH_${G}_COMPLETE"
