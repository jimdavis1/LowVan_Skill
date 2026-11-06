#!/usr/bin/env python3
"""Every member of every collection must have a sibling at the identity floor.

    python3 audit_floor.py --workdir . --module Mastadenovirus
    python3 audit_floor.py ... --write          # then REBUILD

The minimum clustering identity is settled at 0.80. A member whose best
within-bin hit falls below it was never a member of that protein family: it
violates the binning rule and should not have clustered there in the first
place. No genome check and no literature check is owed to such a sequence --
the curator's rule, 5 Nov 2026:

    "In every case the proteins with bad homology within the set should be
     discarded. They violate our binning rules and should not have been
     clustered in the first place."

This is separate from audit_bins.py, which judges whether a BIN should exist.
A bin can be entirely legitimate and still carry a handful of members that do
not belong to it, and those members do three kinds of damage:

  - they widen the profile until it matches the neighbours (Aviadenovirus
    FIBER held one sequence 99.8% identical to a FIBER2 member, and that one
    sequence is why the FIBER profile scored 352 bits on FIBER2's locus and
    93.6 on its own);
  - they make a coherent bin look incoherent, so audit_bins proposes
    retiring a real protein;
  - they survive every other check, because a profile recalls its own
    training data however mixed that data is.

Applied to four flagged Adenoviridae bins it resolved one outright --
Barthadenovirus RH6's three members have best siblings at 74%, 74% and 40%,
so nothing survives the floor -- and shrank the others to a core worth
asking about: UNC58 4 -> 3, UNC29 9 -> 4, UNC5 5 -> 2, UNC19 18 -> 5.

WHAT THIS DELIBERATELY DOES NOT DO. It never retires a whole bin, even when
no member clears the floor; it reports that as EMPTIED and leaves the bin for
audit_bins and the curator. And --min-keep refuses to strip a bin below the
member floor, because a profile built from two survivors is its own problem.

THE FLOOR APPLIES TO BINS WHOSE IDENTITY IS IN QUESTION, NOT TO NAMED
ORTHOLOGUE SETS. This was the scope error, caught in pre-flight on 6 Nov 2026
before anything was written: run blind on Siadenovirus it flagged 24 of 26
bins and proposed dropping 65 sequences, among them members of HEXON, POL,
PENTON, DBP, PROTEASE and FIBER. Those are the core structural genes and they
are correct. Siadenovirus spans frog, turkey, raptor and tortoise
adenoviruses, so a genuine hexon orthologue whose nearest sibling is 73.7%
identical is exactly what belongs there.

0.80 is a CLUSTERING threshold. It governs how the members of a feature are
grouped into sub-clusters, each of which becomes its own PSSM. It is not a
membership test for the collection, and a collection is an orthologue set
across the whole taxon.

So the evidence of identity decides the scope:

  named feature   identity is established by the annotation string AND the
                  ORF's position in the genome. A divergent member is a
                  divergent orthologue. Leave it alone.
  UNC bin         mutual similarity is the ONLY evidence that the members are
                  one protein. A member that does not resemble the others has
                  no claim to membership at all, and no genome or literature
                  check is owed to it.

Hence --unnamed-only is the default. Auditing a named feature takes
--include-named and a reason, and --exempt still removes individual sets.
"""
import argparse, collections, glob, os, shutil, statistics, subprocess
import sys, tempfile

BV = "/Applications/BV-BRC.app/runtime/bin"
if os.path.isdir(BV):
    os.environ["PATH"] = BV + ":" + os.environ.get("PATH", "")


def read_fasta(path):
    out, name, cur = [], None, []
    for line in open(path, errors="replace"):
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(cur)))
            name, cur = line[1:].rstrip("\n"), []
        elif name is not None:
            cur.append(line.strip())
    if name is not None:
        out.append((name, "".join(cur)))
    return out


def collections_dir(workdir, module):
    for d in (os.path.join(workdir, module, "collections"),
              os.path.join(workdir, "collections", module)):
        if os.path.isdir(d):
            return d
    sys.exit("no collections directory for %s under %s" % (module, workdir))


def best_siblings(path, tmp):
    """Best within-bin percent identity for each member, by bitscore.

    A collection can hold the same id twice (one genome depositing two
    copies), so a name match is not proof of self -- but a full-length
    100% hit is, and excluding it is the conservative error.
    """
    seqs = read_fasta(path)
    if len(seqs) < 2:
        return {h.split()[0]: None for h, _ in seqs}
    r = subprocess.run(
        ["blastp", "-query", path, "-subject", path,
         "-outfmt", "6 qseqid sseqid pident length qlen bitscore",
         "-evalue", "10", "-num_threads", "4"],
        capture_output=True, text=True)
    best = {}
    for line in r.stdout.splitlines():
        q, s, pid, aln, qlen, bits = line.split("\t")
        if q == s and float(pid) > 99.5 and int(aln) >= int(qlen):
            continue
        pid, bits = float(pid), float(bits)
        if q not in best or bits > best[q][1]:
            best[q] = (pid, bits)
    return {h.split()[0]: (best[h.split()[0]][0] if h.split()[0] in best else 0.0)
            for h, _ in seqs}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--module", required=True)
    ap.add_argument("--floor", type=float, default=80.0,
                    help="minimum within-bin percent identity (default 80, "
                         "the settled clustering floor)")
    ap.add_argument("--min-keep", type=int, default=3,
                    help="never strip a bin below this many members; report "
                         "it instead (default 3)")
    ap.add_argument("--include-named", action="store_true",
                    help="also audit named features. OFF by default: for a "
                         "named feature, identity is established by its name "
                         "and its genome position, so a divergent member is a "
                         "divergent orthologue, not a binning error. Running "
                         "blind on Siadenovirus flagged HEXON, POL, PENTON, "
                         "DBP and PROTEASE.")
    ap.add_argument("--unnamed-prefix", default="UNC",
                    help="prefix marking a bin whose identity rests only on "
                         "mutual similarity (default UNC)")
    ap.add_argument("--exempt", action="append", default=[],
                    help="a named orthologue set that is genuinely divergent "
                         "across host species. Repeatable.")
    ap.add_argument("--write-md")
    ap.add_argument("--write", action="store_true",
                    help="rewrite each collection without its sub-floor "
                         "members, keeping the original as .prefloor")
    a = ap.parse_args()

    cdir = collections_dir(a.workdir, a.module)
    feats = sorted(os.path.basename(p)[:-6]
                   for p in glob.glob(os.path.join(cdir, "*.fasta"))
                   if not p.endswith(".outliers.fasta"))
    exempt = set(a.exempt)
    tmp = tempfile.mkdtemp(prefix="floor_")
    rows, total_drop = [], 0
    try:
        for f in feats:
            if f in exempt:
                continue
            if not a.include_named and not f.startswith(a.unnamed_prefix):
                continue
            path = os.path.join(cdir, f + ".fasta")
            sib = best_siblings(path, tmp)
            bad = sorted(k for k, v in sib.items()
                         if v is not None and v < a.floor)
            if not bad:
                continue
            n = len(sib)
            keep = n - len(bad)
            rows.append(dict(feat=f, n=n, bad=bad, keep=keep,
                             worst=min(sib[k] for k in bad),
                             status=("EMPTIED" if keep == 0 else
                                     "BELOW-MIN-KEEP" if keep < a.min_keep
                                     else "strip")))
            if rows[-1]["status"] == "strip":
                total_drop += len(bad)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    lines = ["# Identity-floor audit -- %s" % a.module, "",
             "Every member must have a sibling at >= %.0f%% identity." % a.floor,
             "", "Scope: %s." % ("every collection (--include-named)"
                                 if a.include_named else
                                 "%s* bins only -- a named feature's identity "
                                 "rests on its name and genome position, not "
                                 "on mutual similarity" % a.unnamed_prefix),
             "", "**Bins with sub-floor members: %d of %d audited.  Sequences "
             "to drop: %d.**" % (len(rows), sum(
                 1 for f in feats if f not in exempt and
                 (a.include_named or f.startswith(a.unnamed_prefix))),
                 total_drop), ""]
    if exempt:
        lines += ["Exempt (divergent orthologue sets): %s"
                  % ", ".join(sorted(exempt)), ""]
    if rows:
        lines += ["| bin | members | sub-floor | would keep | worst | action |",
                  "|---|---|---|---|---|---|"]
        for r in sorted(rows, key=lambda x: -len(x["bad"])):
            lines.append("| `%s` | %d | %d | %d | %.1f%% | %s |"
                         % (r["feat"], r["n"], len(r["bad"]), r["keep"],
                            r["worst"], r["status"]))
        lines.append("")
        for r in sorted(rows, key=lambda x: -len(x["bad"])):
            if r["status"] != "strip":
                lines += ["**`%s` -> %s**: %d of %d members are sub-floor. "
                          "Not stripped automatically; audit_bins and the "
                          "curator decide this one."
                          % (r["feat"], r["status"], len(r["bad"]), r["n"]), ""]
    else:
        lines += ["Every member of every collection clears the floor.", ""]

    report = "\n".join(lines)
    print(report)
    if a.write_md:
        open(a.write_md, "w").write(report + "\n")

    if a.write:
        moved = 0
        for r in rows:
            if r["status"] != "strip":
                continue
            path = os.path.join(cdir, r["feat"] + ".fasta")
            shutil.copy2(path, path + ".prefloor")
            drop = set(r["bad"])
            with open(path, "w") as o:
                for h, s in read_fasta(path + ".prefloor"):
                    if h.split()[0] in drop:
                        continue
                    o.write(">%s\n%s\n" % (h, s))
            moved += len(r["bad"])
        print("dropped %d sequence(s); originals kept as .prefloor -- "
              "REBUILD profiles and the JSON" % moved)


if __name__ == "__main__":
    main()
