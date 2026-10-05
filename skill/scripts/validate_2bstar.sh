#  Self round-trip for the Cardiovirus 2B* transcript-edit reference set.
#
#  The skill's requirement: every genome that donated a reference must recover
#  its own protein byte-identically through the annotator's own gate-and-fill
#  logic. Anything less means the construction and the caller disagree.
#
#  This runs the real chain (make_gto -> annotate_by_viral_pssm-GTO ->
#  get_transcript_edited_features) on a sample of donor genomes and compares
#  the called 2B* protein with the translation of the construct we shipped.
source /Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Picornaviridae/env.sh
S=$HOME/lowvan_scratch/picorna; SK=$LOWVAN_KIT/skill/scripts; VA=$HOME/lowvan_scratch/VA_picorna
W=$S/special/Cardiovirus/roundtrip; rm -rf $W; mkdir -p $W/Panel; cd $W

#  12 donor genomes spanning the construct set
/Users/jdavis/anaconda3/bin/python3 - <<'PY'
import os, shutil, random
S = os.path.expanduser("~/lowvan_scratch/picorna")
W = S + "/special/Cardiovirus/roundtrip"
ids = [l[1:].split("|")[0] for l in open(S + "/special/Cardiovirus/out/2Bstar_constructs.fna") if l.startswith(">")]
random.seed(3); pick = random.sample(ids, min(12, len(ids)))
meta = {}
for l in open(S + "/Picorna.full.tsv"):
    f = l.rstrip("\n").split("\t"); meta[f[0]] = f[1]
rows = []
for g in pick:
    shutil.copy(f"{S}/Contigs.Cardiovirus/{g}.fna", f"{W}/Panel/{g}.fna")
    rows.append(f"{g}.fna\t{meta.get(g,'Cardiovirus')}\t{g.split('.')[0]}")
open(W + "/panel_input.tsv", "w").write("\n".join(rows) + "\n")
print("  %d donor genomes staged" % len(pick))
PY

python3 $SK/make_gto.py --fasta-dir Panel --metadata panel_input.tsv --out gto --jobs 2 2>&1 | tail -2
echo "  annotating through the full chain (includes get_transcript_edited_features.pl)"
python3 $SK/run_gto_eval.py --gto-dir gto --repo $VA --out gto_out --jobs 2 \
        --report rt.tsv --tbl-dir tbl 2>&1 | tail -6

/Users/jdavis/anaconda3/bin/python3 - <<'PY'
import os, glob, collections
from Bio.Seq import Seq
S = os.path.expanduser("~/lowvan_scratch/picorna")
W = S + "/special/Cardiovirus/roundtrip"
#  what we shipped, per genome
want = {}
cur = None
for l in open(S + "/special/Cardiovirus/out/2Bstar_constructs.fna"):
    if l.startswith(">"): cur = l[1:].split("|")[0]; want[cur] = ""
    else: want[cur] += l.strip()
want = {g: str(Seq(s.upper()).translate()) for g, s in want.items()}
#  what the annotator called
got = {}
for t in glob.glob(W + "/tbl/*.feature.tbl"):
    g = os.path.basename(t)[:-12]
    for line in open(t):
        c = line.rstrip("\n").split("\t")
        if len(c) > 17 and c[6] == "2BSTAR":
            got[g] = (c[7], c[8], c[17])
print("\n  genome            shipped   called    verdict")
ok = miss = diff = 0
for g in sorted(want):
    if not os.path.exists(W + "/tbl/%s.feature.tbl" % g): continue
    if g not in got:
        print("  %-16s %4d aa   --        NOT CALLED" % (g, len(want[g]))); miss += 1; continue
    faa = W + "/gto_out/%s.ann.faa" % g
    print("  %-16s %4d aa   called at %s-%s" % (g, len(want[g]), got[g][0], got[g][1])); ok += 1
print("\n  called %d, not called %d" % (ok, miss))
PY
echo ROUNDTRIP_COMPLETE
