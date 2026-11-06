#!/usr/bin/env python3
"""Re-derive bit_cutoff from what the ANNOTATOR saw, not from the collections.

derive_cutoffs_separation.py measures `foreign` by running each feature's
profiles against other features' curated protein collections. The annotator
runs them against a GENOME in six frames. Those are not the same subject: a
profile reaches a higher bitscore on a genomic translation than on any curated
member, so a cutoff that cleanly separates the collections can leave the
genome unseparated, and the module ships calling two features on one locus.

Measured on Mastadenovirus, 5 Nov 2026, at loci where two profiles fired on
the SAME EXACT SPAN:

    pair              loser   bits at locus   installed cutoff
    E1B19K / UNC34    UNC34        66.3            31
    E1B55K / UNC23    UNC23        75.1            62
    E4ORF6 / UNC19    UNC19       103.0            83
    FIBER  / UNC6     UNC6        202.0           141

Every false call cleared its cutoff. This script closes that by taking the
foreign evidence from the panel: for each co-located locus it decides which
feature OWNS the locus -- by blastp of the called translation against each
feature's own collection, which is `whichwins`' test and does not depend on
the profile scores being compared -- and then records, for the loser, the
bitscore its profile actually reached there.

    foreign_obs(F) = max bits F reached at a locus F does NOT own
    self_obs(F)    = min bits F reached at a locus F DOES own

A proposal is only emitted when foreign_obs < self_obs; otherwise the pair is
genuinely inseparable by score and is reported for a curator to look at,
because raising the cutoff would cost real calls. Nothing is written unless
--write is passed.
"""
import argparse, collections, glob, json, os, subprocess, sys, tempfile


def span(f):
    L = f.get("location") or []
    if not L:
        return None
    c, s, st, d = L[0][0], int(L[0][1]), L[0][2], int(L[0][3])
    return (c, s - d + 1, s, st) if st == "-" else (c, s, s + d - 1, st)


def sym(ft, mod):
    for fa in ft.get("family_assignments") or []:
        if len(fa) > 1 and fa[1].startswith(mod + "."):
            return fa[1].split(".")[1]
    return None


def best_bits(pssm_dir, feat, subject):
    best = 0.0
    for p in sorted(glob.glob(os.path.join(pssm_dir, feat, "*.pssm"))):
        r = subprocess.run(["tblastn", "-in_pssm", p, "-subject", subject,
                            "-outfmt", "6 bitscore", "-max_hsps", "1"],
                           capture_output=True, text=True)
        for ln in r.stdout.split():
            try:
                best = max(best, float(ln))
            except ValueError:
                pass
    return best


def owns(prot, feat, coll_dir, td):
    """Best blastp bits of the called protein against this feature's own set."""
    c = os.path.join(coll_dir, feat + ".fasta")
    if not os.path.exists(c):
        return 0.0
    q = os.path.join(td, "q.fa")
    with open(q, "w") as fh:
        fh.write(">q\n%s\n" % prot)
    r = subprocess.run(["blastp", "-query", q, "-subject", c, "-outfmt",
                        "6 bitscore", "-max_hsps", "1", "-evalue", "10"],
                       capture_output=True, text=True)
    vals = [float(x) for x in r.stdout.split() if x.replace(".", "", 1).isdigit()]
    return max(vals) if vals else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--module", required=True)
    ap.add_argument("--anno", required=True, help="panel <M>.anno directory")
    ap.add_argument("--pssms", required=True, help="<repo>/Viral-PSSMs/<M>.pssms")
    ap.add_argument("--collections", required=True)
    ap.add_argument("--json", required=True, help="module JSON to read/patch")
    ap.add_argument("--margin", type=float, default=1.15,
                    help="a locus is OWNED when its blastp bits beat the rival "
                         "by this factor; below it the pair is called a tie and "
                         "left alone (D=1.15)")
    ap.add_argument("--self-sample", type=int, default=8,
                    help="how many UNCONTESTED loci per candidate feature to "
                         "score for self evidence (D=8). self_obs taken only "
                         "from contested loci is not the feature's true floor: "
                         "most real calls are uncontested and were never "
                         "looked at, so a cutoff raised on it can silently "
                         "delete them.")
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()

    mod = a.module
    J = json.load(open(a.json), object_pairs_hook=collections.OrderedDict)
    feats = J[mod]["features"]

    foreign = collections.defaultdict(float)   # feature -> worst false bits
    selfmin = {}                               # feature -> lowest true bits
    ties, loci = [], 0
    blame = collections.Counter()
    solo = collections.defaultdict(list)       # feature -> uncontested loci

    with tempfile.TemporaryDirectory() as td:
        for gf in sorted(glob.glob(os.path.join(a.anno, "*.ann.gto"))):
            G = json.load(open(gf))
            ctg = {c["id"]: c["dna"] for c in G.get("contigs", [])}
            by = collections.defaultdict(list)
            for ft in G.get("features", []):
                s, k = span(ft), sym(ft, mod)
                if s and k:
                    by[s].append((k, ft.get("protein_translation") or ""))
            for s, rows in by.items():
                if len(rows) < 2:
                    if len(solo[rows[0][0]]) < 400:
                        solo[rows[0][0]].append((gf, s))
                    continue
                loci += 1
                c, b, e, st = s
                sub = os.path.join(td, "s.fna")
                lo, hi = max(0, b - 300), min(len(ctg[c]), e + 300)
                with open(sub, "w") as fh:
                    fh.write(">locus\n%s\n" % ctg[c][lo:hi])
                prot = next((p for _, p in rows if p), "")
                scored = [(k, owns(prot, k, a.collections, td)) for k, _ in rows]
                scored.sort(key=lambda x: -x[1])
                win, lose = scored[0], scored[1]
                if not win[1] or win[1] < a.margin * lose[1]:
                    ties.append((win[0], lose[0]))
                    continue
                bits = best_bits(a.pssms, lose[0], sub)
                foreign[lose[0]] = max(foreign[lose[0]], bits)
                blame[(win[0], lose[0])] += 1
                wb = best_bits(a.pssms, win[0], sub)
                if win[0] not in selfmin or wb < selfmin[win[0]]:
                    selfmin[win[0]] = wb

    #  SELF EVIDENCE FROM UNCONTESTED LOCI.
    #  selfmin above is taken only where two profiles fired. That is a biased
    #  sample: it is the loci where this feature was in a fight, not the loci
    #  where it is simply right. A cutoff set above a feature's true floor
    #  deletes calls, and nothing downstream reports it -- the call rate just
    #  drops. So sample loci where this feature was the ONLY caller and fold
    #  the lowest score in.
    import random
    random.seed(0)
    with tempfile.TemporaryDirectory() as td:
        for k in sorted(set(foreign)):
            rows = solo.get(k) or []
            if not rows:
                continue
            for gf, s in random.sample(rows, min(a.self_sample, len(rows))):
                G = json.load(open(gf))
                ctg = {c["id"]: c["dna"] for c in G.get("contigs", [])}
                c, b, e, st = s
                sub = os.path.join(td, "s.fna")
                lo, hi = max(0, b - 300), min(len(ctg[c]), e + 300)
                with open(sub, "w") as fh:
                    fh.write(">locus\n%s\n" % ctg[c][lo:hi])
                bits = best_bits(a.pssms, k, sub)
                if bits and (k not in selfmin or bits < selfmin[k]):
                    selfmin[k] = bits

    print("## %s -- %d co-located locus(es), %d resolved, %d tie(s)"
          % (mod, loci, sum(blame.values()), len(ties)))
    print("  %-14s %10s %10s %10s %10s  %s"
          % ("feature", "cur_cut", "foreign", "self_obs", "proposed", "false calls for"))
    prop = {}
    for k in sorted(foreign):
        cur = (feats.get(k) or {}).get("bit_cutoff")
        sm = selfmin.get(k)
        new = int(foreign[k]) + 1
        #  No self evidence at all means no floor to test against, and a
        #  cutoff raised blind is a call-rate cut waiting to happen.
        ok = sm is not None and new < sm
        who = ", ".join("%s(%d)" % (w, n) for (w, l), n in blame.items() if l == k)
        print("  %-14s %10s %10.1f %10s %10s  %s"
              % (k, cur, foreign[k], ("%.1f" % sm) if sm else "-",
                 new if ok else ("INSEPARABLE" if sm is not None
                                 else "NO-SELF-EVIDENCE"), who))
        if ok and (cur is None or new > cur):
            prop[k] = new
    if a.write and prop:
        for k, v in prop.items():
            feats[k]["bit_cutoff"] = v
        tmp = a.json + ".tmp.%d" % os.getpid()
        with open(tmp, "w") as fh:
            json.dump(J, fh, indent=3)
            fh.write("\n")
        os.replace(tmp, a.json)
        print("\n  wrote %d raised cutoff(s) -> %s" % (len(prop), a.json))
    elif prop:
        print("\n  %d cutoff(s) would be raised; --write to apply" % len(prop))
    return 0


if __name__ == "__main__":
    sys.exit(main())
