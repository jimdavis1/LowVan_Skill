#!/usr/bin/env python3
"""Add copy_num to the features that earned it, after the module has run.

The safe build order is: ship with no copy_num, measure called/routed on real
genomes, then mark essential only what clears the threshold. copy_num is what
makes a feature essential, so a feature carrying it flags EVERY genome that
lacks it -- three modules scored 0.0% clean that way while their core
features called at 91-100%.

This reads the scored genomes an evaluation produced and rewrites the module
JSON, setting copy_num on features called on >= --threshold of routed
genomes and removing it from those below. It reports both directions and
writes nothing without --write.

    python3 apply_copy_num.py --eval Crinivirus/eval \
            --json Crinivirus/Crinivirus_Viral_PSSM.json --module Crinivirus
    python3 apply_copy_num.py ... --write
"""
import argparse, collections, glob, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from json_canon import write_json


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eval", required=True, help="directory of *.qual.gto")
    ap.add_argument("--json", required=True)
    ap.add_argument("--module", required=True)
    ap.add_argument("--threshold", type=float, default=0.85)
    ap.add_argument("--contain", type=float, default=0.80,
                    help="a call at least this contained in a LONGER feature's "
                         "span does not count toward the call rate")
    ap.add_argument("--standalone", action="store_true",
                    help="decide on the STANDALONE rate, discounting calls "
                         "nested in a longer feature. Off by default because "
                         "genuine overlapping-ORF architecture is common; use "
                         "with --expected. The warning prints either way.")
    ap.add_argument("--expected", default="",
                    help="CHILD:PARENT overlaps that are real genome "
                         "architecture, e.g. CP:POLY,MP:POLY. Not discounted.")
    ap.add_argument("--require", default="",
                    help="comma-separated features the curator holds essential "
                         "on biology, e.g. ENV. Marked whatever their rate and "
                         "never dropped; the rate is still reported.")
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()

    blk = json.load(open(a.json))
    feats = blk[a.module]["features"]

    #  CURATOR-REQUIRED FEATURES BYPASS THE THRESHOLD.
    #
    #  The 85% rule guards against FALSE flags: a genus-specific accessory
    #  marked essential fails good genomes that never had it. It is the wrong
    #  test when the genomes lacking a feature are genuinely defective. Per the
    #  curator, 28 Oct 2026: every replication-competent retrovirus encodes
    #  env, so ENV keeps copy_num 1 in all four retro modules. Betaretrovirus
    #  ENV calls on 50% of exemplars, but 18 of its 26 non-callers are deposits
    #  their own submitters label "truncated env" or "nonfunctional env" --
    #  endogenous proviruses that SHOULD score poor. Without this switch every
    #  later --write silently dropped it again.
    required = {x.strip() for x in a.require.split(",") if x.strip()}
    for k in sorted(required):
        if k not in feats:
            sys.exit("--require %s: no such feature in %s" % (k, a.module))
        if feats[k].get("special") == "splice":
            sys.exit("--require %s: a spliced feature is never essential" % k)

    #  Two rates per feature, per genome:
    #    raw        the feature was called at all
    #    standalone it was called somewhere NOT contained in a longer
    #               feature's span
    #  They differ when a profile fires inside another gene. Closterovirus
    #  P20 is 89.3% raw and 53.3% standalone. But genuine overlapping ORFs
    #  are common -- Capillovirus CP sits inside its own polyprotein and is
    #  called on 100% of genomes with a 2% standalone rate -- so the two
    #  rates are used asymmetrically, see below.
    #  A feature's span is its FULL extent, first exon start to last exon end.
    #  Reading location[0] alone measured a spliced feature by its first exon,
    #  and for HIV2 REV that is a 61-70 nt fragment of a 324 nt, 107 aa
    #  protein. Scored that way REV sits inside TAT's first exon in 35 of 35
    #  genomes and reported 0.0% standalone -- a number about one fifth of the
    #  feature. EIAV REV read 3.3% standalone by exon 1 and 23.3% by full
    #  extent. Neither verdict changed here, but the figure was wrong and it
    #  was the figure that got reported to the curator.
    def exons(f):
        """the feature's exon intervals, which is what overlap must use.

        The SPAN of a spliced feature covers its intron, and the intron is not
        part of the feature. Deltaretrovirus G4 is join(431..503, 7071..7315),
        so its span is 6.9 kb of an 8.8 kb provirus and GAG, PRO, POL and ENV
        all fall inside it -- while sharing not one base with either exon.
        Scored that way ENV read 12.1% standalone against a 93.9% call rate
        and was refused copy_num it had plainly earned, and PRO was recorded
        as a separate gene "inside G4" that it never touches. Measuring by
        exon, ENV is 93.9% standalone. The earlier fix here, from location[0]
        to the full extent, was correcting the opposite error on HIV2 REV and
        went one step too far: neither endpoint pair is the feature, the exons
        are.
        """
        locs = f.get("location") or []
        ivs = []
        for l in locs:
            try:
                c, st0, strand, d0 = l[0], int(l[1]), l[2], int(l[3])
            except Exception:
                continue
            ivs.append((c, st0 - d0 + 1, st0) if strand == "-"
                       else (c, st0, st0 + d0 - 1))
        if not ivs or len({i[0] for i in ivs}) != 1:
            return None
        return ivs

    def span(f):
        locs = f.get("location") or []
        ivs = []
        for l in locs:
            try:
                c, st0, strand, d0 = l[0], int(l[1]), l[2], int(l[3])
            except Exception:
                continue
            ivs.append((c, st0 - d0 + 1, st0) if strand == "-"
                       else (c, st0, st0 + d0 - 1))
        if not ivs:
            return None
        #  Exons of one feature are on one contig; if they are not, the call is
        #  broken and scoring its span is meaningless.
        if len({i[0] for i in ivs}) != 1:
            return None
        return (ivs[0][0], min(i[1] for i in ivs), max(i[2] for i in ivs))

    expected = set()
    for pair in a.expected.split(","):
        if ":" in pair:
            c, par = pair.split(":", 1)
            expected.add((c.strip(), par.strip()))

    #  Frame, for telling a separate gene from a cleavage product.
    #
    #  Per the curator, 26 Oct: a feature that is present in every genome, is
    #  not a mat_peptide of a larger ORF, and has its own start or lies in a
    #  different reading frame, is a real gene and should be marked and
    #  checked for -- even though it overlaps something longer.
    #
    #  That distinction is measurable. A cleavage product is a sub-segment of
    #  its parent and therefore IN FRAME with it. A distinct overlapping gene
    #  is translated from its own start codon, usually in another frame.
    #  EIAV S2 is a median 87.4% inside env and is offset from env's start by
    #  a non-multiple of three in 29 of 29 genomes -- it is never in env's
    #  frame, so it cannot be a piece of env being re-detected.
    def frame_offset(child, parent, strand):
        """0 when the two are in the same reading frame."""
        if strand == "-":
            return (parent[2] - child[2]) % 3
        return (child[1] - parent[1]) % 3

    def strand_of(f):
        locs = f.get("location") or []
        return locs[0][2] if locs and len(locs[0]) > 2 else "+"

    def _cod(ivs):
        """coding length of an exon list -- the feature, not its footprint."""
        return sum(i[2] - i[1] + 1 for i in ivs) if ivs else 0

    anno2key = {v["anno"]: k for k, v in feats.items()}

    #  WHICH MODULE DID THIS GENOME ROUTE TO?
    #
    #  `anno` is not unique across modules. FIV, SRLV and EIAV all call their
    #  features "Gag polyprotein", "Pol polyprotein", "Envelope glycoprotein"
    #  and "Virion infectivity factor Vif", so matching on `anno` alone counts
    #  every genome in the evaluation directory toward every module -- in both
    #  the numerator and the denominator. Run against one Quality_GTOs holding
    #  all three lentivirus genera on 26 Oct, that reported FIV ENV at 99.0%
    #  (really 100/101 = 51 SRLV + 39 FIV + 10 EIAV), EIAV S2 at 10.9% when it
    #  is 11/11, and SRLV TAT at 48.5% when it is 49/51. The S2 figure would
    #  have dropped a copy_num the curator had just decided to keep.
    #
    #  family_assignments carries the module: "SRLV.GAG.2" -> ("SRLV", "GAG").
    #  Fall back to `anno` when a feature has no assignment, which keeps the
    #  old behaviour for a single-module directory.
    def mod_and_key(f):
        fam = (f.get("family_assignments") or [[None, "", "", ""]])[0]
        key = (fam[1] or "") if len(fam) > 1 else ""
        if key:
            bits = key.split(".")
            if len(bits) >= 2 and bits[0] and bits[1]:
                return bits[0], bits[1]
        return None, anno2key.get(f.get("function") or "")

    raw = collections.Counter(); alone = collections.Counter(); routed = 0
    skipped_other = 0
    #  per feature: how many contained calls were out of frame with the
    #  containing feature, and which features did the containing
    oof = collections.Counter(); contained = collections.Counter()
    parents = collections.defaultdict(collections.Counter)
    for q in sorted(glob.glob(os.path.join(a.eval, "*.qual.gto"))):
        try:
            d = json.load(open(q))
        except Exception:
            continue
        if not d.get("features"):
            continue
        parsed = [(mod_and_key(f), span(f), strand_of(f), exons(f))
                  for f in d["features"]]
        seen_mods = collections.Counter(m for (m, _k), _s, _d, _e in parsed if m)
        gmod = seen_mods.most_common(1)[0][0] if seen_mods else None
        if gmod and gmod != a.module:
            skipped_other += 1
            continue
        routed += 1
        fs = [(k, sp, sd, ex) for (m, k), sp, sd, ex in parsed
              if k and sp and ex and (m is None or m == a.module)]
        got, got_alone = set(), set()
        for i, (k, sp, sd, ex) in enumerate(fs):
            got.add(k)
            nested = False
            for j, (k2, sp2, sd2, ex2) in enumerate(fs):
                if i == j or k2 == k or sp[0] != sp2[0]:
                    continue
                if _cod(ex2) <= _cod(ex):
                    continue
                ov = sum(max(0, min(u[2], v[2]) - max(u[1], v[1]) + 1)
                         for u in ex for v in ex2)
                mine = _cod(ex)
                if mine > 0 and ov / mine >= a.contain and (k, k2) not in expected:
                    nested = True
                    contained[k] += 1
                    parents[k][k2] += 1
                    #  Opposite strands are never the same frame.
                    if sd != sd2 or frame_offset(sp, sp2, sd) != 0:
                        oof[k] += 1
                    break
            if not nested:
                got_alone.add(k)
        for k in got:
            raw[k] += 1
        for k in got_alone:
            alone[k] += 1
    called = raw

    if not routed:
        sys.exit("no scored genome in %s carried any feature" % a.eval)

    #  ASYMMETRIC, deliberately.
    #    ADD    needs BOTH rates over the threshold. A feature whose raw rate
    #           clears it only because the profile fires inside a longer gene
    #           would otherwise earn essentiality on another feature's locus.
    #    DROP   uses the RAW rate only. A low standalone rate is just as often
    #           real overlapping-ORF architecture, and silently demoting a
    #           correct feature is worse than leaving one un-demoted.
    #  SPLICED FEATURES ARE NEVER ESSENTIAL.
    #
    #  copy_num makes a feature essential, so every genome lacking it is
    #  flagged. A spliced feature is not called by a PSSM at all -- it is
    #  called by get_splice_variant_features.pl from an exact donor/acceptor
    #  motif lookup, and that lookup misses a real fraction of genomes for
    #  reasons that are about the vocabulary, not about the genome. HIV-1 tat
    #  and rev reproduce the deposited protein on 75% and 71% of held-out
    #  genomes; marking them essential would condemn a quarter of HIV-1 as
    #  poor on a splice-reference gap.
    #
    #  Per the curator, 26 Oct: essentiality is for the principal ORFs that
    #  should occur in 100% of a module's genomes. Every splice feature in the
    #  repo was unmarked except HIV2 TAT, which this tool had added.
    #
    #  transcript_edit is deliberately NOT covered here: Closterovirus RDRP
    #  and Orthoebolavirus GPC are principal ORFs that happen to be reached by
    #  a frameshift, and they are correctly essential.
    add, drop, keep, blocked, spliced = [], [], [], [], []
    outofframe = []
    for k, v in sorted(feats.items()):
        frac = raw[k] / routed
        sa = alone[k] / routed
        has = bool(v.get("copy_num"))
        if v.get("special") == "splice":
            spliced.append((k, frac, has))
            continue
        if k in required:
            (keep if has else add).append((k, frac) if has else (k, frac, sa))
            continue
        if frac >= a.threshold and not has:
            if sa >= a.threshold:
                add.append((k, frac, sa))
            #  The curator's rule, applied mechanically: a contained feature
            #  whose calls are consistently OUT OF FRAME with whatever
            #  contains them is a separate gene, not a cleavage product, so
            #  the containment is real architecture and must not block it.
            #  mat_peptide is excluded outright -- that type IS a piece of a
            #  larger ORF, which is the one case this must never wave through.
            elif (v.get("feature_type") != "mat_peptide"
                  and v.get("type") != "mat_peptide"
                  and contained[k]
                  and oof[k] == contained[k]):
                add.append((k, frac, sa))
                outofframe.append((k, frac, sa, parents[k].most_common(1)[0][0]))
            else:
                blocked.append((k, frac, sa))
        elif frac < a.threshold and has:
            drop.append((k, frac))
        elif has:
            keep.append((k, frac))

    if skipped_other:
        print("  %d genome(s) routed to another module, not counted"
              % skipped_other)
    print("  %s: %d routed genome(s), threshold %.0f%%%s"
          % (a.module, routed, 100 * a.threshold,
             ""))
    #  Warn either way. A feature whose raw rate clears the threshold on
    #  another feature's locus is the failure this exists to prevent --
    #  Closterovirus P20 is 89.3% raw and 53.3% standalone -- but genuine
    #  overlapping ORFs are common (Capillovirus CP sits inside its own
    #  polyprotein), so the default does NOT silently demote them.
    print("\n  %-12s %7s  %s" % ("feature", "called", "copy_num"))
    for k in sorted(feats, key=lambda x: -called[x]):
        frac = called[k] / routed
        sa = alone[k] / routed
        mark = ("splice/never" if feats[k].get("special") == "splice" else
                ("ADD" if not feats[k].get("copy_num") else "keep") + " (required)"
                if k in required else
                "ADD" if any(x[0] == k for x in add) else
                "DROP" if (k, frac) in drop else
                "BLOCKED" if any(x[0] == k for x in blocked) else
                "keep" if feats[k].get("copy_num") else "-")
        print("  %-12s %6.1f%%  %s" % (k, 100 * frac, mark))

    if outofframe:
        print("\n  marked essential DESPITE containment -- every call is out of")
        print("  frame with the feature containing it, so it is a separate gene")
        print("  and not a cleavage product of one:")
        for k, f, sa, par in outofframe:
            print("    %-11s %5.1f%% raw, %5.1f%% standalone, inside %s, 0 in frame"
                  % (k, 100 * f, 100 * sa, par))

    if blocked:
        print("\n  NOT marked essential -- raw rate clears the threshold but the")
        print("  standalone rate does not, so the calls are on another feature's")
        print("  locus. Resolve with qc_duplicate_calls.py, or pass --expected")
        print("  CHILD:PARENT if the overlap is real genome architecture:")
        for k, f, sa in blocked:
            print("    %-11s %5.1f%% raw  ->  %5.1f%% standalone" % (k, 100 * f, 100 * sa))

    #  A splice feature that already carries copy_num is stripped, not left
    #  alone: it was marked in error and leaving it would keep flagging
    #  genomes the splice vocabulary simply did not reach.
    unmark = [k for k, _f, has in spliced if has]
    if spliced:
        print("\n  spliced -- never marked essential (called by")
        print("  get_splice_variant_features.pl, not by a PSSM):")
        for k, f, has in spliced:
            print("    %-11s %5.1f%% called%s"
                  % (k, 100 * f, "   STRIPPING copy_num" if has else ""))

    if a.write:
        for k, _f, _s in add:
            feats[k]["copy_num"] = 1
        for k, _f in drop:
            feats[k].pop("copy_num", None)
        for k in unmark:
            feats[k].pop("copy_num", None)
        write_json(blk, a.json)
        print("\n  wrote %s: +%d copy_num, -%d, %d blocked, %d splice stripped"
              % (a.json, len(add), len(drop), len(blocked), len(unmark)))
    else:
        print("\n  dry run -- rerun with --write (+%d, -%d, %d blocked)"
              % (len(add), len(drop), len(blocked)))


main()
