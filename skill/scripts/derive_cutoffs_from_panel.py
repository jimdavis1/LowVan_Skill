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
    ap.add_argument("--self-sample", type=int, default=0,
                    help="how many UNCONTESTED loci per candidate feature to "
                         "score for self evidence. 0 (the default) scores "
                         "EVERY one. self_obs taken only from contested loci "
                         "is not the feature's true floor -- most real calls "
                         "are uncontested and were never looked at -- and "
                         "neither is a small random sample of the "
                         "uncontested ones. This defaulted to 8, and on "
                         "Mastadenovirus UNC8 those 8 put the floor at 119 "
                         "against a foreign 117, so a cutoff of 118 looked "
                         "separable; scoring all 37 loci put the real floor "
                         "at 116, BELOW the foreign score. The pair is "
                         "inseparable and the 118 would have deleted real "
                         "calls with nothing but the call rate to show it. A "
                         "partial sample is still allowed, but it can only "
                         "report -- no cutoff is proposed from one.\n")
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--cutoffs",
                    help="build/<M>/bit_cutoffs.json -- the SIDECAR that "
                         "setup_build.py reads back. Required with --write, "
                         "because the module JSON is regenerated from "
                         "features.json and a cutoff written only there is "
                         "erased by the next build.")
    a = ap.parse_args()
    #  Checked HERE, not at the write: this script blastp's every contested
    #  locus and every self sample first, which is minutes of work, and a
    #  refusal that arrives after all of it is a refusal nobody waits for.
    if a.write and not a.cutoffs:
        sys.exit(
            "REFUSED: --write without --cutoffs.\n"
            "  A cutoff written only to the module JSON is erased the next\n"
            "  time setup_build.py runs, because it reads cutoffs back from\n"
            "  build/<module>/bit_cutoffs.json. Pass that path as --cutoffs\n"
            "  so the measurement survives a rebuild.")

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
    seen = {}                                  # feature -> (scored, available)
    with tempfile.TemporaryDirectory() as td:
        for k in sorted(set(foreign)):
            rows = solo.get(k) or []
            if not rows:
                continue
            n = len(rows) if a.self_sample <= 0 else min(a.self_sample, len(rows))
            seen[k] = (n, len(rows))
            for gf, s in random.sample(rows, n):
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
    print("  %-14s %10s %10s %10s %8s %10s  %s"
          % ("feature", "cur_cut", "foreign", "self_obs", "sampled",
             "proposed", "false calls for"))
    prop = {}
    partial = []
    for k in sorted(foreign):
        cur = (feats.get(k) or {}).get("bit_cutoff")
        sm = selfmin.get(k)
        new = int(foreign[k]) + 1
        n, tot = seen.get(k, (0, 0))
        #  No self evidence at all means no floor to test against, and a
        #  cutoff raised blind is a call-rate cut waiting to happen. Neither
        #  is a PARTIAL floor evidence: the lowest of 8 sampled loci is an
        #  upper bound on the true floor, never the floor itself, so a
        #  cutoff that clears it may still sit above real calls that were
        #  not looked at. Report those, do not propose from them.
        ok = sm is not None and new < sm
        if ok and tot and n < tot:
            partial.append(k)
            ok = False
        who = ", ".join("%s(%d)" % (w, n2) for (w, l), n2 in blame.items() if l == k)
        print("  %-14s %10s %10.1f %10s %8s %10s  %s"
              % (k, cur, foreign[k], ("%.1f" % sm) if sm else "-",
                 ("%d/%d" % (n, tot)) if tot else "-",
                 new if ok else ("PARTIAL-SAMPLE" if k in partial else
                                 "INSEPARABLE" if sm is not None
                                 else "NO-SELF-EVIDENCE"), who))
        if ok and (cur is None or new > cur):
            prop[k] = new
    if partial:
        print("\n  %d feature(s) looked separable but were scored on only part\n"
              "  of their uncontested loci, so the floor is an upper bound and\n"
              "  no cutoff is proposed: %s\n"
              "  Re-run with --self-sample 0 to score every locus."
              % (len(partial), ", ".join(partial)))
    if a.write and prop:
        #  THE MODULE JSON IS AN OUTPUT, NOT A STORE.
        #
        #  setup_build.py rebuilds features.json from the collections and
        #  re-attaches each bit_cutoff from build/<M>/bit_cutoffs.json;
        #  gen_module_json.py then writes the module JSON from features.json,
        #  falling back to a LENGTH-BASED default -- 60 bits at 250 residues --
        #  for anything the spec does not carry. That default is a floor for
        #  "callable at all", not a threshold that discriminates, so a cutoff
        #  measured here and written only to the JSON does not come back as an
        #  error on the next regeneration: it silently reverts, and the module
        #  ships calling two features on one locus again. This is the same
        #  shape as the `special` flag that four generators used to delete --
        #  a curation decision living only in a generated file -- and it is
        #  why --cutoffs is required rather than optional.
        if not a.cutoffs:
            sys.exit(
                "REFUSED: --write without --cutoffs.\n"
                "  A cutoff written only to %s is erased the next time\n"
                "  setup_build.py runs, because it reads cutoffs back from\n"
                "  build/<module>/bit_cutoffs.json. Pass that path as\n"
                "  --cutoffs so the measurement survives a rebuild." % a.json)
        side = {}
        if os.path.exists(a.cutoffs):
            try:
                side = json.load(open(a.cutoffs))
            except Exception:
                side = {}
        for k, v in prop.items():
            feats[k]["bit_cutoff"] = v
            side[k] = v
        tmp = a.json + ".tmp.%d" % os.getpid()
        with open(tmp, "w") as fh:
            json.dump(J, fh, indent=3)
            fh.write("\n")
        os.replace(tmp, a.json)
        tmp = a.cutoffs + ".tmp.%d" % os.getpid()
        with open(tmp, "w") as fh:
            json.dump(side, fh, indent=1, sort_keys=True)
            fh.write("\n")
        os.replace(tmp, a.cutoffs)
        print("\n  wrote %d raised cutoff(s) -> %s\n  and -> %s"
              % (len(prop), a.json, a.cutoffs))
    elif prop:
        print("\n  %d cutoff(s) would be raised; --write to apply" % len(prop))
    return 0


if __name__ == "__main__":
    sys.exit(main())
