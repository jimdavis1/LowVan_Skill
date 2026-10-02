#!/usr/bin/env python3
"""Can this feature's profiles actually call it? Measure, do not assume.

The assumption this script exists to replace: "a short protein cannot be given
a bit_cutoff that fires on its own members and stays silent across a whole
genome." It is a reasonable prior and it is sometimes wrong.

Enterovirus VPg is 22 aa. Measured against 300 full-length polyproteins, over
5,052 profile x subject pairs:

    at the VPg locus      10th pct 21.2 bits, median 27.0
    anywhere else          99th pct 17.6 bits, MAX 19.3

A cutoff of 20 calls it and is silent across the other 2,160 residues. The
reason is specific to the protein, not to its length -- VPg is highly
conserved and nothing else in the polyprotein resembles it. A 22-mer that is
merely short and generic would have failed the same test, which is the point
of running it.

What this measures, and why it is the right test: the annotator scores each
profile against a whole translated genome, so the question is never "does the
profile recover its training data" (self-recall always says yes) but "does it
score the right locus above everything else in the same sequence".

    python3 test_profile_specificity.py --pssm-dir Alignments/M/VPG/pssms \
            --subjects poly.faa --locus 1460-1580

Without --locus the script takes the modal hit position as the locus and
reports how tightly the profiles agree, which is itself the evidence that
there IS one locus.
"""
import argparse, collections, glob, os, subprocess, sys, tempfile

ap = argparse.ArgumentParser()
ap.add_argument("--pssm-dir", required=True)
ap.add_argument("--subjects", required=True,
                help="FASTA of sequences that CONTAIN the feature -- "
                     "full-length precursors or translated genomes, not the "
                     "feature's own collection")
ap.add_argument("--locus", default=None, help="LO-HI residue window of the true locus")
ap.add_argument("--locus-ids", default=None,
                help="File of subject ids that ARE the feature, one per line. "
                     "Use this when the feature is a standalone protein rather "
                     "than a domain of a precursor: the subjects are then the "
                     "whole protein complement of genomes that carry it, and "
                     "on-target means the profile hit the right PROTEIN, not "
                     "the right residue window. Same question, different axis.")
ap.add_argument("--window", type=int, default=60,
                help="half-width around the modal position when --locus is absent")
a = ap.parse_args()

locus_ids = set()
if a.locus_ids:
    locus_ids = {l.split()[0] for l in open(a.locus_ids) if l.strip()}
    if not locus_ids: sys.exit("--locus-ids file is empty")

td = tempfile.mkdtemp()
subprocess.run(["makeblastdb", "-in", a.subjects, "-dbtype", "prot",
                "-out", td + "/db"], check=True, stdout=subprocess.DEVNULL)
pssms = sorted(glob.glob(os.path.join(a.pssm_dir, "*.pssm")))
if not pssms: sys.exit("no PSSMs in %s" % a.pssm_dir)

raw = []
for p in pssms:
    #  -comp_based_stats 0 matters: the composition correction rescales short
    #  profiles hardest, and the annotator does not apply it.
    out = subprocess.run(["psiblast", "-in_pssm", p, "-db", td + "/db",
                          "-comp_based_stats", "0", "-evalue", "1000",
                          "-outfmt", "6 sseqid sstart send bitscore",
                          "-max_target_seqs", "5000"],
                         capture_output=True, text=True).stdout
    for l in out.splitlines():
        s, ss, se, b = l.split("\t")
        raw.append((os.path.basename(p), s, int(ss), float(b)))
if not raw: sys.exit("no hits at all -- the profiles cannot see these subjects")

lo = hi = None
if locus_ids:
    missing = locus_ids - {r[1] for r in raw}
    print("on-target axis: subject identity. %d of %d feature-bearing subjects "
          "were hit by at least one profile%s"
          % (len(locus_ids) - len(missing), len(locus_ids),
             "; %d never hit" % len(missing) if missing else ""))
elif a.locus:
    lo, hi = (int(x) for x in a.locus.split("-"))
else:
    best = collections.defaultdict(list)
    for prof, s, st, b in raw: best[(prof, s)].append((b, st))
    tops = [max(v)[1] for v in best.values()]
    tops.sort(); mode = tops[len(tops) // 2]
    lo, hi = mode - a.window, mode + a.window
    print("locus taken from the data: modal top-hit start %d -> window %d-%d" % (mode, lo, hi))

def ontarget(subj, st):
    return subj in locus_ids if locus_ids else (lo <= st <= hi)

per = collections.defaultdict(list)
for prof, s, st, b in raw: per[(prof, s)].append((b, st))
on, off, found = [], [], 0
for (prof, subj), v in per.items():
    i = [b for b, st in v if ontarget(subj, st)]
    o = [b for b, st in v if not ontarget(subj, st)]
    if i: on.append(max(i)); found += 1
    if o: off.append(max(o))
on.sort(); off.sort()

#  The pair view above is not the question the annotator asks. It runs EVERY
#  profile over a genome and calls the feature if ANY of them clears the
#  cutoff, so the decision is per subject, not per profile x subject: does the
#  BEST on-target score beat the BEST off-target score anywhere in the same
#  genome? A feature whose profiles are divergent alleles of one gene fails the
#  pair view by construction -- allele A's profile scores near-zero on allele
#  B's protein, which the pair view records as a weak on-target hit, when in
#  truth B is called by B's own profile. Velarivirus P4 failed the pair view at
#  10th pct 16.3 vs off-target max 24.1 and passed this one comfortably.
sbest_on, sbest_off = collections.defaultdict(float), collections.defaultdict(float)
for prof, subj, st, b in raw:
    d = sbest_on if ontarget(subj, st) else sbest_off
    d[subj] = max(d[subj], b)
gon = sorted(sbest_on.values()); goff = sorted(sbest_off.values())
def pct(v, q): return v[min(len(v) - 1, int(q * len(v)))] if v else 0.0

print("%d profile x subject pairs with a hit; %d (%.0f%%) hit the target"
      % (len(per), found, 100.0 * found / max(1, len(per))))
print("  on target    : min %.1f  10th %.1f  median %.1f  max %.1f"
      % (on[0], pct(on, .10), pct(on, .50), on[-1]) if on else "  on target    : none")
if off:
    print("  elsewhere    : median %.1f  90th %.1f  99th %.1f  MAX %.1f"
          % (pct(off, .50), pct(off, .90), pct(off, .99), off[-1]))
else:
    print("  elsewhere    : no off-target hit at all")
floor = pct(on, .10); ceil = off[-1] if off else 0.0
if gon:
    print()
    print("per-subject view (what the annotator actually does -- best profile wins)")
    print("  best on target  : min %.1f  10th %.1f  median %.1f"
          % (gon[0], pct(gon, .10), pct(gon, .50)))
    if goff:
        print("  best off target : median %.1f  99th %.1f  MAX %.1f"
              % (pct(goff, .50), pct(goff, .99), goff[-1]))
    gfloor, gceil = gon[0], (goff[-1] if goff else 0.0)
    if gfloor > gceil:
        print("  -> CALLABLE at bit_cutoff %d: every feature-bearing subject scores "
              "above it (worst %.1f) and nothing else reaches it (best %.1f)."
              % (int(gceil) + 1, gfloor, gceil))
    else:
        bad = sum(1 for v in gon if v <= gceil)
        print("  -> a single cutoff loses %d of %d subjects (worst on-target %.1f "
              "vs best off-target %.1f)." % (bad, len(gon), gfloor, gceil))

print()
print("pair view (every profile judged on its own):")
if floor > ceil:
    print("CALLABLE. Suggested bit_cutoff %d: above every off-target hit (%.1f) "
          "and below the 10th percentile on-target (%.1f)."
          % (int(ceil) + 1, ceil, floor))
else:
    print("NOT SEPARABLE by a single cutoff. 10th pct on-target %.1f is at or "
          "below the worst off-target %.1f -- a cutoff that calls the feature "
          "will also fire elsewhere. Derive it positionally, or do not declare it."
          % (floor, ceil))
