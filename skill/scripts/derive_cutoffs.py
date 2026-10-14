#!/usr/bin/env python3
"""Set each feature's bit_cutoff from what its own profiles actually score.

`gen_module_json.py` defaults the cutoff from length alone -- 60 bits at 250
residues, 30 at 80, 25 below -- which is a floor for "callable at all", not a
threshold that discriminates. Scaling a shipped module's cutoffs by length is
no better: Sapovirus got cutoffs scaled from the homologous Norovirus features
and SIX OF NINE sat above their own weakest training member, so the feature
could not call sequences it was built from. NS1 was called on 136 of 218
genomes and NS2 on 177, for no reason but the threshold.

A cutoff has to sit below the weakest member the profile must still find, and
above whatever noise the taxon generates. This measures the first directly:
every profile is run against its own collection and the distribution of
self-scores is read off.

**The weakest self-hit is not the weakest member.** A profile scores its own
training data many times over, and some of those hits are short spurious
matches against a repeat or a low-complexity run elsewhere in a longer
sequence. Sapovirus NS4's weakest raw self-hit was 20.4 bits; its weakest hit
that actually covered the protein was 61.6, three times higher. Taking the raw
minimum would have set the cutoff at 12 and invited exactly the failure this
project already hit once, where a 49-residue self-match at 69.6 bits sat
beside a real one at 2,627. So a hit counts only if it spans `--min-cov` of
the subject.

    python3 derive_cutoffs.py --workdir . --module Sapovirus            # report
    python3 derive_cutoffs.py --workdir . --module Sapovirus \
            --features features_all.json --write

`--margin` is how far below the weakest genuine member to sit. **0.9, and the
default used to be 0.6, which was wrong in a way worth recording** because the
failure looks like success in the call counts.

At 0.6 the Sapovirus call rates did exactly what was wanted -- contigs carrying
all seven ORF1 products went from 133 of 218 to 214 -- and the boundaries fell
apart underneath:

    NS1 median length         60 -> 36     (Mc10 69)
    NS1 starts at POLY start  96.2% -> 5.6%
    NS1/NS2 junction overlap  0/133 -> 198/214
    sum(products) - POLY      -1 -> -32 aa

A 69-residue profile at 16 bits matches almost anything, so a spurious internal
match outscored the true N-terminal one and the feature was called in the wrong
place on nearly every genome. The p05 of every other product degraded with it:
NS4 274 -> 193, NS5 114 -> 81. **Counting calls cannot see this. Only the
tiling test can**, which is why `check_cleaved_ends.py` and the sum-minus-
precursor check belong in the same loop as this script.

So: sit just below the weakest member the profile must still find, not far
below it. Lower the margin only for a feature that is genuinely absent from
most genomes and whose boundaries do not matter, and re-run the tiling test
afterwards either way.
"""
import argparse, glob, json, os, statistics, subprocess, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--module", required=True)
    ap.add_argument("--json", help="module JSON (default <workdir>/<M>_Viral_PSSM.json)")
    ap.add_argument("--features", help="features_all.json to update with --write")
    ap.add_argument("--min-cov", type=float, default=0.6,
                    help="hit must span this fraction of the subject (default 0.6)")
    ap.add_argument("--margin", type=float, default=0.9,
                    help="cutoff = margin x weakest genuine member (default 0.9; 0.6 broke Sapovirus NS1 -- see the docstring)")
    ap.add_argument("--max-profiles", type=int, default=4,
                    help="profiles sampled per feature (default 4)")
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()

    W, M = a.workdir, a.module
    jpath = a.json or os.path.join(W, "%s_Viral_PSSM.json" % M)
    if not os.path.exists(jpath):
        sys.exit("no module JSON at %s" % jpath)
    feats = json.load(open(jpath))[M]["features"]

    print("  %-8s %6s | %9s %9s | %9s | %s"
          % ("feature", "now", "raw min", "genuine", "gen p05", "proposed"))
    out, clipping = {}, []
    for key in sorted(feats):
        pssms = sorted(glob.glob(os.path.join(W, "Alignments", M, key, "pssms", "*.pssm")))
        coll = os.path.join(W, "collections", M, key + ".fasta")
        if not pssms or not os.path.exists(coll):
            continue
        rows = []
        for p in pssms[:a.max_profiles]:
            r = subprocess.run(["psiblast", "-in_pssm", p, "-subject", coll,
                                "-comp_based_stats", "0",
                                "-outfmt", "6 bitscore length slen",
                                "-max_target_seqs", "5000"],
                               capture_output=True, text=True)
            for line in r.stdout.splitlines():
                q = line.split("\t")
                if len(q) >= 3:
                    try:
                        rows.append((float(q[0]), int(q[1]), int(q[2])))
                    except ValueError:
                        pass
        if not rows:
            continue
        raw = sorted(x[0] for x in rows)
        genuine = sorted(b for b, ln, sl in rows if sl and ln >= a.min_cov * sl)
        if not genuine:
            print("  %-8s  no hit spans %.0f%% of its subject -- skipped"
                  % (key, 100 * a.min_cov))
            continue
        now = feats[key].get("bit_cutoff")
        proposed = max(10, int(genuine[0] * a.margin))
        out[key] = proposed
        flag = ""
        if now is not None and now >= genuine[0]:
            flag = "   <- CLIPS its own members"
            clipping.append(key)
        print("  %-8s %6s | %9.1f %9.1f | %9.1f | %d%s"
              % (key, now, raw[0], genuine[0],
                 genuine[int(.05 * len(genuine))], proposed, flag))

    if clipping:
        print("\n  %d feature(s) cannot call sequences they were built from: %s"
              % (len(clipping), ", ".join(clipping)))
    else:
        print("\n  no feature clips its own training data")

    if a.write:
        if not a.features:
            sys.exit("--write needs --features")
        spec = json.load(open(a.features))
        n = 0
        for k, v in out.items():
            if k in spec:
                spec[k]["bit_cutoff"] = v
                n += 1
        json.dump(spec, open(a.features, "w"), indent=1)
        print("  wrote bit_cutoff on %d feature(s) -> %s" % (n, a.features))
        print("  now regenerate the module JSON and re-measure coverage")
    else:
        print("  (report only -- rerun with --write --features features_all.json)")


main()
