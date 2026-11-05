#!/usr/bin/env python3
"""Length-audit EVERY collection, against the literature, BEFORE profiles.

`capsid_length_audit.py` does this for nested alternative-start forms, which
is where the defect was first met. It is not where the defect lives. Any
collection can hold a fragment, a fusion, or two genes at once, and none of
the downstream checks can see it:

  * a PSSM built from a mixed bin recalls its own training data perfectly,
  * the call rate is unaffected or improves,
  * the annotation string is correct,
  * and the only symptom is duplicate calls -- which look like a cutoff
    problem, and get "fixed" by spending a coverage run on cutoffs.

Adenoviridae, 4-5 Nov 2026, and it is the SECOND bin contamination in that
one family. `FIBER` held a **5 aa** "fiber protein" (28285.923, HAdV-5)
alongside a 1500 aa member; `FIBER1` and `FIBER2` each held the other's gene.
That is what drove Mastadenovirus FIBER's weakest-member score to 22.8 bits
and put FIBER/FIBER1/FIBER2 at 95-99% identity to each other -- read, wrongly
and expensively, as "one protein split three ways" only after two full
coverage panels had been run on it.

So the audit is not a diagnosis you reach for when something looks odd. It
runs at bin time, and `install_module.py` refuses a module without a clean
and current report.

WHAT IT TESTS

  fragment   a member far below its own species' median. The partial-deposit
             filter cannot catch these: the deposit string is a plain
             "fiber protein" with no `fragment` or `truncated` marker, so
             nothing upstream objects.
  giant      a member far above it -- a fusion, a mis-join, or two ORFs.
  split      one species showing two separated length classes inside ONE
             collection. That is the signature of two genes in one bin, and
             it is the finding that matters most, because homology cannot
             undo it later.
  literature a member outside the range a paper gives for that feature.
             Every range needs a PMID; a range without one is not evidence
             and the tool refuses it.
  unverified a feature with no literature entry at all. Reported, not fatal,
             so that what has been checked is distinguishable from what has
             merely not failed yet.

SPECIES KEY

The genome id prefix IS the NCBI taxon id, so `10524.1311` -> `10524`, exact
and free. Do not key on the genome NAME: "Human adenovirus 41 Hannover_2022"
and "Human adenovirus 41 DOM03186" are one species under two strings, and
keying on the name reported 53 within-species matches as cross-species while
this very defect was being diagnosed.

Compare per species, never per module or genus -- the whole point of
`capsid_length_audit.py`'s docstring, and it holds here. Fowl adenovirus
makes the stronger version of the argument: fiber-1 is the LONG fiber in
FAdV-A (793 aa vs 410) and the SHORT one in FAdV-C (431 vs 479), so a
module-wide length rule does not merely blur the two genes, it swaps them
(PMID 35070851; CELO genome PMID 8627769).

    python3 audit_collection_lengths.py --workdir . --module Mastadenovirus \\
            --refs length_refs.json --write-md LENGTH_AUDIT.md

Exit status is 1 if anything unresolved is reported, so it gates a build.
"""
import argparse, collections, glob, json, os, re, sys


def read_fasta(path):
    out, h, s = [], None, []
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if h is not None:
                    out.append((h, "".join(s)))
                h, s = line[1:].strip(), []
            else:
                s.append(line.strip())
    if h is not None:
        out.append((h, "".join(s)))
    return out


def genome_of(h):
    return h.split()[0].split("|")[0]


def taxon_of(h):
    """'10524.1311 long fiber protein' -> '10524'. The exact species key."""
    return genome_of(h).split(".")[0]


def median(xs):
    xs = sorted(xs)
    n = len(xs)
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2.0


def classes(lengths, gap=0.25):
    """Split sorted lengths where a neighbour jumps more than `gap`.

    Two genes in one bin show as two clumps with a wide space between, not
    as a smooth spread. A relative gap is what detects both the 227 aa
    VP1/VP2 offset and the 42 aa Amdoparvovirus one -- no absolute
    threshold separates those two at once.
    """
    out, cur = [], []
    for x in sorted(lengths):
        if cur and (x - cur[-1]) > gap * max(cur[-1], 1):
            out.append(cur)
            cur = []
        cur.append(x)
    if cur:
        out.append(cur)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--module", required=True)
    ap.add_argument("--collections",
                    help="collections dir. Default: <workdir>/<module>/"
                         "collections, else <workdir>/collections/<module>")
    ap.add_argument("--refs",
                    help="JSON: {feature: {min, max, pmid, note, "
                         "per_taxon:{taxid:[min,max]}}}. Every entry needs a "
                         "pmid -- a range with no paper behind it is a guess "
                         "wearing a citation's clothes.")
    ap.add_argument("--frag", type=float, default=0.60,
                    help="member below this fraction of its species median "
                         "is a fragment (default 0.60)")
    ap.add_argument("--giant", type=float, default=1.60,
                    help="member above this multiple of it is a giant")
    ap.add_argument("--min-species", type=int, default=3,
                    help="species with fewer members than this are judged "
                         "against the whole collection instead, since a "
                         "median of one or two says nothing")
    ap.add_argument("--accept",
                    help="TSV of findings a curator has signed off, one per "
                         "line: KIND <tab> FEATURE <tab> WHO <tab> REASON. "
                         "WHO is the member id or taxon the finding names, or "
                         "'*' for every finding of that kind in that feature. "
                         "An accepted finding stays in the report and stops "
                         "counting as unresolved.")
    ap.add_argument("--write-md")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()

    W, M = a.workdir, a.module
    cdir = a.collections
    for cand in (os.path.join(W, M, "collections"),
                 os.path.join(W, "collections", M)):
        if cdir:
            break
        if os.path.isdir(cand):
            cdir = cand
    if not cdir or not os.path.isdir(cdir):
        sys.exit("no collections dir for %s under %s" % (M, W))

    refs = {}
    if a.refs and os.path.exists(a.refs):
        refs = json.load(open(a.refs))
        bad = [k for k, v in refs.items() if not str(v.get("pmid", "")).strip()]
        if bad:
            sys.exit("refs entries with no pmid: %s\n"
                     "  Every range must name the paper it came from. See the "
                     "literature-authority rule in the skill." % ", ".join(sorted(bad)))

    #  A signed-off finding is not a silenced one. It stays in the report
    #  with the reason beside it, because the whole failure this tool exists
    #  for was a defect that nothing printed -- and an exception list that
    #  deletes its own entries is the same failure wearing a process.
    accepted = {}
    if a.accept and os.path.exists(a.accept):
        for line in open(a.accept, errors="replace"):
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                sys.exit("--accept needs 4 columns: %s" % line.strip())
            accepted[(f[0].strip(), f[1].strip(), f[2].strip())] = f[3].strip()

    def signed(kind, feat, who):
        return (accepted.get((kind, feat, who))
                or accepted.get((kind, feat, "*"))
                or accepted.get((kind, feat, who.split()[0] if who else "")))

    def report(kind, feat, who, detail):
        why = signed(kind, feat, who)
        (okfindings if why else findings).append((kind, feat, who, detail, why))

    findings = []   # (kind, feature, who, detail, accepted_reason)
    okfindings = []
    unverified = []
    summary = []

    for path in sorted(glob.glob(os.path.join(cdir, "*.fasta"))):
        feat = os.path.basename(path)[:-6]
        rows = read_fasta(path)
        if not rows:
            continue
        alllen = [len(s) for _, s in rows]
        allmed = median(alllen)
        ref = refs.get(feat)
        if ref is None:
            unverified.append(feat)

        #  A length that several unrelated members share is a FORM of the
        #  protein, not a defect in one deposit -- the same rule the binner
        #  uses to decide what is an unmarked partial, and it has to be the
        #  same rule or the two disagree about the same sequence. 197
        #  E4ORF67 deposits at exactly 83 aa are a second form; a lone 5 aa
        #  "fiber protein" is a fragment. Reporting the form once as a SPLIT
        #  puts the decision where it belongs; reporting it as 197 FRAGMENTs
        #  buries it.
        nlen = collections.Counter(len(s) for _, s in rows)
        FORM = 3

        bytax = collections.defaultdict(list)
        for h, s in rows:
            bytax[taxon_of(h)].append((h, len(s)))

        nfrag = ngiant = nlit = 0
        splits = []
        for tax, mem in sorted(bytax.items()):
            L = [n for _, n in mem]
            #  A one- or two-member species has no usable median of its own;
            #  judging it against itself would clear every fragment that
            #  happens to be alone in its species.
            med = median(L) if len(L) >= a.min_species else allmed
            for h, no in mem:
                if nlen[no] >= FORM:
                    continue
                if no < a.frag * med:
                    report("FRAGMENT", feat, h,
                           "%d aa vs %s median %d" % (no, tax, med))
                    nfrag += 1
                elif no > a.giant * med:
                    report("GIANT", feat, h,
                           "%d aa vs %s median %d" % (no, tax, med))
                    ngiant += 1
            if len(L) >= 2 * a.min_species:
                cl = classes(L)
                #  Ignore a clump of one -- that is the fragment/giant test's
                #  job and it has already run. Two populated classes is the
                #  two-genes-in-one-bin signature.
                cl = [c for c in cl if len(c) >= 2]
                if len(cl) > 1:
                    splits.append((tax, [(len(c), min(c), max(c)) for c in cl]))
                    report("SPLIT", feat, tax,
                           " + ".join("%dx %d-%d aa" % (len(c), min(c), max(c))
                                      for c in cl))

            if ref:
                lo, hi = ref.get("min"), ref.get("max")
                pt = (ref.get("per_taxon") or {}).get(tax)
                if pt:
                    lo, hi = pt
                if lo is not None and hi is not None:
                    for h, no in mem:
                        if no < lo or no > hi:
                            report("LITERATURE", feat, h,
                                   "%d aa outside %d-%d aa (PMID %s)"
                                   % (no, lo, hi, ref["pmid"]))
                            nlit += 1

        summary.append((feat, len(rows), min(alllen), int(allmed), max(alllen),
                        nfrag, ngiant, len(splits), nlit,
                        "" if ref is None else "PMID %s" % ref["pmid"]))

    #  ---------------------------------------------------------------- out
    lines = []
    w = lines.append
    w("# Length audit -- %s" % M)
    w("")
    w("Collections: `%s`" % cdir)
    w("")
    w("| feature | n | min | median | max | frag | giant | split | lit | reference |")
    w("|---|---|---|---|---|---|---|---|---|---|")
    for r in summary:
        w("| `%s` | %d | %d | %d | %d | %s | %s | %s | %s | %s |"
          % (r[0], r[1], r[2], r[3], r[4],
             r[5] or "", r[6] or "", r[7] or "", r[8] or "", r[9] or "—"))
    w("")
    if unverified:
        w("## Features with no literature range (%d)" % len(unverified))
        w("")
        w("Not a defect. It records that these have been checked for internal")
        w("consistency only, so that what is verified stays distinguishable")
        w("from what has merely not failed yet.")
        w("")
        w("`" + "`, `".join(unverified) + "`")
        w("")
    if findings:
        w("## Findings (%d)" % len(findings))
        w("")
        bykind = collections.defaultdict(list)
        for kind, feat, who, detail, _why in findings:
            bykind[kind].append((feat, who, detail))
        for kind in ("SPLIT", "LITERATURE", "FRAGMENT", "GIANT"):
            if not bykind[kind]:
                continue
            w("### %s (%d)" % (kind, len(bykind[kind])))
            w("")
            for feat, who, detail in bykind[kind][:200]:
                w("- `%s`  %s  — %s" % (feat, who, detail))
            if len(bykind[kind]) > 200:
                w("- ... %d more" % (len(bykind[kind]) - 200))
            w("")
        w("**Unresolved: %d.** Resolve each one by re-binning, by dropping the"
          % len(findings))
        w("member, or by widening the literature range with a PMID that")
        w("supports it -- then re-run. Do not build profiles over this.")
    else:
        w("## Findings")
        w("")
        w("**Unresolved: 0.**")
    if okfindings:
        w("")
        w("## Signed off (%d)" % len(okfindings))
        w("")
        w("Checked and accepted, with the reason. These still describe the")
        w("collections -- they are not suppressed, only answered.")
        w("")
        w("| kind | feature | member | measurement | why it is correct |")
        w("|---|---|---|---|---|")
        for kind, feat, who, detail, why in okfindings[:400]:
            w("| %s | `%s` | %s | %s | %s |" % (kind, feat, who, detail, why))
        if len(okfindings) > 400:
            w("")
            w("... %d more" % (len(okfindings) - 400))
    out = "\n".join(lines) + "\n"

    if a.write_md:
        open(a.write_md, "w").write(out)
    if not a.quiet:
        sys.stdout.write(out)
    return 1 if findings else 0


if __name__ == "__main__":
    sys.exit(main())
