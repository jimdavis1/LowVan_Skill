#!/usr/bin/env python3
"""Set bit_cutoff to the lowest value that admits no OTHER feature's protein.

`derive_cutoffs.py` sets the cutoff to `--margin` x the weakest genuine member
of the feature's own collection. Its own docstring says why that is only a
ceiling: it is derived from the training collection, so it says what the
profiles need to find THEMSELVES, and wild sequences are routinely more
divergent than anything in the collection. Applying it to Adenoviridae cost
exactly that -- duplicate co-located loci fell 4,103 -> 784, and with them

    Aviadenovirus   6,601 -> 3,800 genome-level calls   (-42%)
    Mastadenovirus  7,299 -> 5,346                      (-27%)
    LIPASE            71  -> 15 genomes      I_LEADER  84 -> 17
    Mastadenovirus UNC4, UNC6, UNC53 -> called on NOTHING

The docstring also names the rule that works, and this script is that rule
applied WITHIN a module rather than across taxa. The question a cutoff has to
answer is not "how low can I go before I stop recognising my own training
data" but "how low can I go before I start calling something that belongs to
another feature". For a within-module duplicate the competitor is not a
foreign genus -- it is the feature next door.

So, per feature:

    self    = min over members of (max over its own profiles), counting only
              hits that span --min-cov of the subject.   [same as derive_cutoffs]
    foreign = the highest score this feature's profiles reach against any
              OTHER feature's collection in the same module.

    cutoff  = foreign + 1            when foreign + 1 <= self
            = int(margin x self)     when foreign >= self   -- INSEPARABLE
                                     (and the pair is reported: no threshold
                                      divides them, so they must be merged or
                                      documented as paralogues)

The first case is the point. A feature with no within-module competitor gets a
LOW cutoff and keeps its divergent wild calls; only a feature that actually
collides with a neighbour pays for it. Under the margin rule every feature
paid, whether it collided or not.

    python3 derive_cutoffs_separation.py --workdir build/<M> --module <M> \
            --json <M>_Viral_PSSM.json
    ... --write-json          patch the runtime JSON in place (preserves
                              `special` declarations and PMID lists, which
                              regenerating the JSON would erase)

Requires psiblast. Profiles live at
<workdir>/Alignments/<M>/<FEAT>/pssms/*.pssm and collections at
<workdir>/collections/<M>/<FEAT>.fasta, the same layout derive_cutoffs reads.
"""
import argparse, collections, glob, json, os, subprocess, sys, tempfile


def psi(pssm, subject, maxt=5000):
    r = subprocess.run(["psiblast", "-in_pssm", pssm, "-subject", subject,
                        "-comp_based_stats", "0",
                        "-outfmt", "6 sseqid bitscore length slen",
                        "-max_target_seqs", str(maxt)],
                       capture_output=True, text=True)
    out = []
    for line in r.stdout.splitlines():
        q = line.split("\t")
        if len(q) >= 4:
            try:
                out.append((q[0], float(q[1]), int(q[2]), int(q[3])))
            except ValueError:
                pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--module", required=True)
    ap.add_argument("--json", required=True)
    ap.add_argument("--min-cov", type=float, default=0.6)
    ap.add_argument("--margin", type=float, default=0.9,
                    help="fallback multiplier when no threshold separates (default 0.9)")
    ap.add_argument("--max-profiles", type=int, default=4)
    ap.add_argument("--write-json", action="store_true")
    ap.add_argument("--sidecar", help="also write {feature: cutoff} here, so "
                    "regenerating the module JSON keeps the measured values")
    a = ap.parse_args()

    W, M = a.workdir, a.module
    feats = json.load(open(a.json))[M]["features"]

    colls, pssms = {}, {}
    for k in sorted(feats):
        c = os.path.join(W, "collections", M, k + ".fasta")
        p = sorted(glob.glob(os.path.join(W, "Alignments", M, k, "pssms", "*.pssm")))
        if os.path.exists(c) and p:
            colls[k], pssms[k] = c, p[:a.max_profiles]

    print("  %-10s %8s %8s | %8s | %s" % ("feature", "self", "foreign", "cutoff", "collides with"))
    out, inseparable = {}, []
    for k in sorted(colls):
        #  self: every member must be callable by SOME profile of this feature
        bysub = {}
        for p in pssms[k]:
            for sid, b, ln, sl in psi(p, colls[k]):
                if sl and ln >= a.min_cov * sl and b > bysub.get(sid, -1):
                    bysub[sid] = b
        if not bysub:
            print("  %-10s  no hit spans %.0f%% of its subject -- skipped"
                  % (k, 100 * a.min_cov))
            continue
        self_ = min(bysub.values())

        #  foreign: the best this feature's profiles do on anyone else's
        #  collection. One concatenated subject keeps it to one psiblast per
        #  profile instead of one per (profile, feature) pair.
        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as t:
            for j, c in colls.items():
                if j == k:
                    continue
                for line in open(c, errors="replace"):
                    t.write(line.replace(">", ">%s|" % j, 1) if line[0] == ">" else line)
            other = t.name
        foreign, worst = 0.0, None
        for p in pssms[k]:
            for sid, b, ln, sl in psi(p, other):
                if b > foreign:
                    foreign, worst = b, sid.split("|")[0]
        os.unlink(other)

        #  `foreign + 1` alone can land at 1-25 bits for a feature with no
        #  within-module competitor, and the annotator searches a ~35 kb
        #  GENOME in six frames, not a protein set. gen_module_json.py's
        #  length default is the floor for "callable at all"; nothing below
        #  it is a threshold.
        med = feats[k].get("median_aa") or 0
        floor = 60 if med >= 250 else 30 if med >= 80 else 25
        if max(foreign + 1, floor) <= self_:
            cut = max(int(foreign) + 1, floor)
            note = (worst or "") + ("" if int(foreign) + 1 >= floor else " [floor]")
        elif foreign >= self_:
            #  No threshold divides them: the competitor outscores this
            #  feature's own weakest member.
            cut = max(10, int(self_ * a.margin))
            note = "%s INSEPARABLE" % (worst or "")
            inseparable.append((k, worst, round(self_, 1), round(foreign, 1)))
        else:
            #  Separable, but the length floor sits above the weakest member --
            #  the collection holds something far more divergent than its
            #  median length implies. Sit below that member, not on the floor.
            cut = max(10, int(self_ * a.margin))
            note = (worst or "") + " [floor %d > self]" % floor
        out[k] = cut
        print("  %-10s %8.1f %8.1f | %8d | %s" % (k, self_, foreign, cut, note))

    if inseparable:
        print("\n  %d feature(s) that no threshold separates -- merge them or"
              " document them as paralogues:" % len(inseparable))
        for k, w, s, f in inseparable:
            print("    %-12s self=%-8s foreign=%-8s (%s)" % (k, s, f, w))

    if a.write_json:
        J = json.load(open(a.json), object_pairs_hook=collections.OrderedDict)
        F = J[M]["features"]
        n = 0
        for k, v in out.items():
            if F.get(k, {}).get("special"):
                continue
            F[k]["bit_cutoff"] = v
            n += 1
        json.dump(J, open(a.json, "w"), indent=3)
        open(a.json, "a").write("\n")
        print("\n  wrote bit_cutoff on %d feature(s) -> %s" % (n, a.json))

    if a.sidecar:
        json.dump(out, open(a.sidecar, "w"), indent=1, sort_keys=True)
        print("  wrote %d cutoff(s) -> %s" % (len(out), a.sidecar))


if __name__ == "__main__":
    main()
