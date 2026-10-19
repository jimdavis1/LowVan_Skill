#!/usr/bin/env python3
"""Find the lowest bit_cutoff that admits no protein from another taxon.

`derive_cutoffs.py` answers "how low can this cutoff go before the profiles
stop finding their own training data". That is a ceiling, not an answer. Wild
sequences are routinely more divergent than anything in the collection, so a
cutoff set from the training data alone can sit above what real genomes score
and quietly cost coverage -- re-cutting the five Caliciviridae modules from
that statistic alone dropped Sapovirus NS2 and NS3 by 26 points each and
Vesivirus VP2 from 100% to 45%, while gaining nothing on POLY in either
Sapovirus or Lagovirus because nothing foreign was anywhere near it.

The question a cutoff has to answer is the other one: what is the lowest value
that still admits no protein from a different taxon? That is measurable. Score
the feature's own profiles against every protein in the family, keep the hits
that fall inside the feature's declared min_len/max_len -- because the
annotator applies that window too, and a cutoff only has to exclude what the
window lets through -- and take the highest score reached by anything from
another group.

    cutoff = the curated value           if it already exceeds that score
           = highest foreign score + 1   otherwise
           = margin x weakest member     when foreign outscores the weakest
                                         member, so no cutoff separates them

The third case is not a failure. Sapovirus NS5, VP1 and VP2 and Vesivirus VP2
all land there: a short or well-conserved protein simply is not separable from
its relatives by score. Genome routing is what keeps a norovirus out of the
sapovirus module; the cutoff is the second line of defence, not the first.

    python3 foreign_leak.py --workdir sapo --module Sapovirus \
            --db /tmp/calidb --group-field genus
    python3 foreign_leak.py ... --apply features_all.json

`--db` is a protein BLAST database and `--map` a TSV of <subject-id> <TAB>
<group> saying which taxon each subject belongs to; any subject whose group
differs from --module (and is not blank/unknown) counts as foreign. Build both
from the dump -- see references/inputs.md.

**Put only FOREIGN sequences in --db.** Member sequences are scored and then
thrown away: `maxf` keeps only hits whose group differs from the module, and
the weakest-member figure comes from the collection FASTA, not from the
database. So every member in --db is pure cost, paid once per profile.

That is cheap on a small module and ruinous on a large one. HIV-1 has 795
profiles, and a whole-family database of 92,248 sequences had the run at
2 hours and a third of the way through ENV; the foreign-only database is
5,002 sequences -- 18x smaller -- and returns bit-identical cutoffs.

Keep the full membership in --map even so. The map is what tells the run
which hits are kin, and the kin check below reads it.
"""
import argparse, glob, json, os, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def best_hits(pssms, db, max_targets):
    """{subject: (best bitscore, aligned length)} over all of a feature's profiles."""
    best = {}
    for p in pssms:
        r = subprocess.run(
            ["psiblast", "-in_pssm", p, "-db", db, "-comp_based_stats", "0",
             "-evalue", "10", "-outfmt", "6 sseqid bitscore length",
             "-num_alignments", str(max_targets)],
            capture_output=True, text=True)
        for line in r.stdout.splitlines():
            f = line.split()
            if len(f) < 3:
                continue
            try:
                b, ln = float(f[1]), int(f[2])
            except ValueError:
                continue
            if f[0] not in best or b > best[f[0]][0]:
                best[f[0]] = (b, ln)
    return best


def weakest_member(pssms, collection, min_cov):
    """Lowest score among members, each scored by its BEST covering profile.

    Same statistic as derive_cutoffs: min over members of max over profiles.
    Taking the global minimum instead would pick up a cross-clade hit from a
    profile that is not responsible for that member.
    """
    if not os.path.exists(collection):
        return 0.0
    member = {}
    for p in pssms:
        r = subprocess.run(
            ["psiblast", "-in_pssm", p, "-subject", collection,
             "-comp_based_stats", "0", "-outfmt", "6 sseqid bitscore length slen",
             "-max_target_seqs", "5000"], capture_output=True, text=True)
        for line in r.stdout.splitlines():
            f = line.split()
            if len(f) < 4:
                continue
            try:
                b, ln, sl = float(f[1]), int(f[2]), int(f[3])
            except ValueError:
                continue
            if sl and ln >= min_cov * sl and b > member.get(f[0], -1.0):
                member[f[0]] = b
    return min(member.values()) if member else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--module", required=True)
    ap.add_argument("--json", help="module JSON (default <workdir>/<module>_Viral_PSSM.json)")
    ap.add_argument("--db", required=True,
                help="protein BLAST db of FOREIGN sequences only -- members are "
                     "scored then discarded, so including them is pure cost")
    ap.add_argument("--map", required=True, help="TSV: subject-id <TAB> group")
    ap.add_argument("--margin", type=float, default=0.9,
                    help="fallback multiplier on the weakest member when foreign "
                         "cannot be separated (default 0.9)")
    ap.add_argument("--min-cov", type=float, default=0.6)
    ap.add_argument("--max-targets", type=int, default=6000)
    ap.add_argument("--apply", help="feature spec JSON to rewrite with the result")
    a = ap.parse_args()

    W, M = a.workdir, a.module
    jpath = a.json or os.path.join(W, "%s_Viral_PSSM.json" % M)
    feats = json.load(open(jpath))[M]["features"]

    group = {}
    for line in open(a.map):
        c = line.rstrip("\n").split("\t")
        if len(c) >= 2:
            group[c[0]] = c[1].strip()

    print("  %-10s %7s %10s %9s %9s  %s"
          % ("feature", "curated", "maxForeign", "weakest", "-> cutoff", "basis"))
    out, unsep, kin = {}, [], []
    for k in sorted(feats):
        pssms = sorted(glob.glob(os.path.join(W, "Alignments", M, k, "pssms", "*.pssm")))
        if not pssms:
            continue
        cur = feats[k].get("bit_cutoff")
        lo, hi = feats[k].get("min_len", 0), feats[k].get("max_len", 10 ** 9)
        best = best_hits(pssms, a.db, a.max_targets)
        foreign = [b for s, (b, ln) in best.items()
                   if lo <= ln <= hi
                   and group.get(s, "") not in ("", M, "(none)", "unknown")]
        maxf = max(foreign) if foreign else 0.0
        weak = weakest_member(pssms, os.path.join(W, "collections", M, k + ".fasta"),
                              a.min_cov)
        if cur is not None and cur > maxf:
            cut, basis = cur, "curated value already excludes all foreign"
        elif weak and maxf + 1 > weak:
            cut, basis = max(10, int(weak * a.margin)), "NOT SEPARABLE -- routing is the defence"
            unsep.append(k)
            #  Before believing that, ask whether the top "foreign" hit is
            #  really foreign. Every NOT SEPARABLE verdict seen while building
            #  the Spinareoviridae modules was a bad --map, not a property of
            #  the protein: sibling dumps hold records of the module's OWN
            #  taxon under other names, and scoring them as foreign puts the
            #  bar above the members. Pteropine orthoreovirus is deposited as
            #  Nelson Bay, Melaka, Kampar, Pulau, Paguyaman, Xi River, Kasama,
            #  HK23629 and Lopburi, scattered across three dumps; a name-based
            #  map missed Lopburi01 and Lopburi02 at 95-97% identity and made
            #  sigma B look inseparable when its real cutoff is 277 against a
            #  weakest member of 696. A hit this close is kin, so say so.
            top = max((b, s) for s, (b, ln) in best.items()
                      if lo <= ln <= hi
                      and group.get(s, "") not in ("", M, "(none)", "unknown"))
            kin.append((k, top[1], group.get(top[1], "?")))
        else:
            cut, basis = int(maxf) + 1, "lowest non-leaking value"
        out[k] = max(10, int(cut))
        print("  %-10s %7s %10.0f %9.0f %9d  %s"
              % (k, cur, maxf, weak, out[k], basis))

    if kin:
        print("\n  Before accepting those, CHECK THE MAP. The sequence that set the bar")
        print("  is named below with the group --map assigns it. If it is really this")
        print("  module's own taxon under another name -- sibling dumps routinely hold")
        print("  such records -- it is kin, not foreign, and the verdict is an artefact:")
        for k, sid, grp in kin:
            print("     %-10s top 'foreign' hit %-44s as %s" % (k, sid[:44], grp))
    if unsep:
        print("\n  %d feature(s) cannot be separated from other taxa by score: %s"
              % (len(unsep), ", ".join(unsep)))
        print("  That is a statement about the protein, not a defect. Check that "
              "genome routing\n  sends those taxa elsewhere, and say so in the "
              "coverage audit.")

    if a.apply:
        spec = json.load(open(a.apply))
        n = 0
        for k, v in spec.items():
            if k in out and v.get("bit_cutoff") != out[k]:
                v["bit_cutoff"] = out[k]
                n += 1
        json.dump(spec, open(a.apply, "w"), indent=1)
        print("\n  rewrote bit_cutoff on %d feature(s) -> %s" % (n, a.apply))
        print("  regenerate the module JSON, reinstall, and RE-RUN THE PANEL: a "
              "cutoff change\n  moves call rates in both directions.")
    else:
        print("\n  (report only -- rerun with --apply <feature spec> to write)")
    print(json.dumps(out))


main()
