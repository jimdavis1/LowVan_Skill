#  Build Cardiovirus LSTAR (overlapping ORF, ordinary PSSM) and 2BSTAR
#  (-1 PRF transframe product, special: transcript_edit, no PSSM).
#
#  L*   Kong & Roos 1991 PMID 2033677; Chen 1995 PMID 7585219
#  2B*  Loughran 2011 PMID 22025686; Finch 2015 PMID 26063423
source /Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Picornaviridae/env.sh
set -e
S=$HOME/lowvan_scratch/picorna; SK=$LOWVAN_KIT/skill/scripts
G=Cardiovirus; B=$S/build/$G; SP=$S/special/$G
PY=/Users/jdavis/anaconda3/bin/python3

echo "=== extract from all $(ls $SP/all/ann/*.feature.tbl | wc -l) annotated genomes"
$PY $S/special_cardio.py $SP/all/ann $S/Contigs.$G $SP/out

echo "=== stage the module build tree"
rm -rf $B; mkdir -p $B/collections/$G $B/Alignments/$G
rsync -a --exclude 'tmp/' $S/fin/$G/collections/$G/ $B/collections/$G/
rm -f $B/collections/$G/L.fasta   # pre-rename duplicate of Ldr.fasta, not a declared feature
rsync -a --exclude 'tmp/' $S/fin/$G/Alignments/$G/  $B/Alignments/$G/
rsync -a /Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Picornaviridae/work/$G/Rep-Contigs $B/
cp $LOWVAN_KIT/modules/$G/${G}_Viral_PSSM.json $B/

#  L* is an ordinary collection: one sequence per distinct ORF
cp $SP/out/Lstar.faa $B/collections/$G/LSTAR.fasta
echo "  LSTAR collection: $(grep -c '>' $B/collections/$G/LSTAR.fasta) sequences"

echo "=== build the LSTAR profiles"
cd $B
cat > feat_lstar.json <<'EOF'
{"LSTAR": {"anno": "L* protein", "params": {"m": 5, "mi": 0.8, "mc": 0.8}}}
EOF
python3 $SK/build_features.py --workdir . --module $G --features feat_lstar.json 2>&1 | tail -4
python3 $SK/build_leftover_pssms.py --workdir . --module $G --key LSTAR --write 2>&1 \
  | grep -E "leftover sequence|built|nothing" | sed 's/^/  /'

echo "=== 2BSTAR transcript-editing reference set"
TE=$HOME/lowvan_scratch/VA_picorna/Transcript-Editing/$G
mkdir -p $TE
cp $SP/out/2Bstar_constructs.fna $TE/2BSTAR.fasta
echo "  $(grep -c '>' $TE/2BSTAR.fasta) nucleotide references (one per genome)"

echo "=== declare both features"
$PY - <<'PY'
import json, os, collections
G = "Cardiovirus"
B = os.path.expanduser("~/lowvan_scratch/picorna/build/%s" % G)
SP = os.path.expanduser("~/lowvan_scratch/picorna/special/%s/out" % G)
p = os.path.join(B, "%s_Viral_PSSM.json" % G)
d = json.load(open(p)); F = d[G]["features"]

def lens(fa):
    L, cur = [], 0
    for line in open(fa):
        if line.startswith(">"):
            if cur: L.append(cur)
            cur = 0
        else: cur += len(line.strip())
    if cur: L.append(cur)
    return sorted(L)

#  L*: its own AUG and its own stop, overlapping POLY in another frame.
#  CDS convention in this project is upstream_ext 1 / downstream_ext 0.
L = lens(os.path.join(SP, "Lstar.faa")); med = L[len(L)//2]
F["LSTAR"] = {
    "anno": "L* protein", "gene_symbol": "L*",
    "feature_type": "CDS", "type": "CDS",
    "bit_cutoff": max(20, int(round(0.45*med))), "coverage_cutoff": 0.65,
    "min_len": int(round(0.8*med)), "max_len": int(round(1.2*med)),
    "upstream_ext": 1, "downstream_ext": 0,
    "kmers": 0,                      # overlaps POLY in another frame
    "segment": "Single RNA Segment",
    "PMID_claude_generated": {"coordinates": ["2033677"], "function": ["7585219"]},
}
#  2B*: -1 PRF, not collinear -> no PSSM. N-terminus is the 2A/2B cleavage, so
#  mat_peptide, exactly as Togaviridae TF. Lengths include the trailing stop.
P = lens(os.path.join(SP, "2Bstar.faa"))
F["2BSTAR"] = {
    "anno": "TransFrame protein", "gene_symbol": "2B*",
    "feature_type": "mat_peptide", "type": "mat_peptide",
    "min_len": P[0] - 2, "max_len": P[-1] + 2,
    "kmers": 0, "segment": "Single RNA Segment",
    "special": "transcript_edit",
    "PMID_claude_generated": {"coordinates": ["22025686"], "function": ["26063423"]},
}
def norm(x):
    if isinstance(x, dict): return {k: norm(v) for k, v in x.items()}
    if isinstance(x, list): return [norm(v) for v in x]
    if isinstance(x, float) and x == int(x): return int(x)
    return x
open(p, "w").write(json.dumps(norm(d), indent=3, separators=(",", " : "),
                              sort_keys=True, ensure_ascii=False) + "\n")
print("  LSTAR  n=%d median %d  window %d-%d  bit %d"
      % (len(L), med, F["LSTAR"]["min_len"], F["LSTAR"]["max_len"], F["LSTAR"]["bit_cutoff"]))
print("  2BSTAR n=%d lengths %d-%d (incl stop)  window %d-%d"
      % (len(P), P[0], P[-1], F["2BSTAR"]["min_len"], F["2BSTAR"]["max_len"]))
PY

echo "=== validate and install"
python3 $SK/json_canon.py --check ${G}_Viral_PSSM.json | tail -1
python3 $SK/install_module.py --workdir $B --repo $HOME/lowvan_scratch/VA_picorna --only $G 2>&1 | tail -8
cp $B/${G}_Viral_PSSM.json $S/
echo BUILD_CARDIO_SPECIAL_COMPLETE
