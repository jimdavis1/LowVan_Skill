#!/usr/bin/env python3
"""Every CDS profile must begin at the protein's start. Flag the ones that do not.

A mature peptide begins at a protease cut, so a non-Met first residue is correct
and expected there. A CDS has no such excuse: a translated protein starts at Met,
so a CDS alignment whose first column is not Met is wrong until shown otherwise.

The upstream pipeline already half-knows this. `process_alignment` trims ragged
N-terminal columns, records where the trim stopped, and when it stopped on a
non-Met residue writes the cluster to FLAGGED_NTERM -- "for review, not trimmed
further". Nothing consumes that file. Norovirus POLY shipped with eight flagged
alignments, cluster 1 of 2598 sequences trimmed at column 349, and the annotator
then called the polyprotein a median 127 codons late on 409 genomes while every
mature peptide inside it tiled correctly.

FLAGGED_NTERM is also incomplete. It is written during the first pass, so the
alignments produced by N-terminal re-clustering and by the leftover pass are
never tested: Norovirus POLY cluster 4 was flagged, re-clustered into 4_1, 4_2
and 4_3, and all three inherited an 80-residue trim with no flag anywhere.

This checks the alignments that actually became installed profiles, and separates
the two failure modes, because they have different fixes:

  TRIMMED   the alignment starts N residues into its own members' sequences.
            The profile cannot match the first N residues, so the annotator's
            hit begins late and scan_to_met_start then stops at whatever
            internal Met it meets first. Fix by splitting the cluster on its
            N-terminus, not by trimming further.

  NON_MET   the alignment starts at offset 0 but the members themselves do not
            begin with Met. The sequences are N-terminally incomplete in the
            source annotation. Fix in the collection, by dropping sequences that
            another sequence extends.

    python3 check_met_start.py --workdir . --module Norovirus
    python3 check_met_start.py --workdir . --module Norovirus --key POLY --repo ~/VA

With --repo the check is restricted to profiles that are actually installed;
without it, every alignment on disk is tested.
"""
import argparse, collections, glob, json, os, statistics, sys


def read_fasta(path):
    out, h, s = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if h is not None:
                out.append((h, "".join(s)))
            h, s = line.split()[0][1:], []
        else:
            s.append(line.strip())
    if h is not None:
        out.append((h, "".join(s)))
    return out


def modal_first_residue(seqs):
    """First column with >=50% occupancy, and its modal residue."""
    if not seqs:
        return None, 0.0
    width = max(len(x) for x in seqs)
    for c in range(width):
        col = [x[c] for x in seqs if c < len(x) and x[c] != "-"]
        if len(col) / len(seqs) >= 0.5:
            res, n = collections.Counter(col).most_common(1)[0]
            return res, n / len(col)
    return None, 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--module", required=True)
    ap.add_argument("--key", default="", help="one feature; default every CDS")
    ap.add_argument("--repo", default="", help="only test installed profiles")
    ap.add_argument("--json", default="", help="module JSON; default <workdir>/<M>_Viral_PSSM.json")
    ap.add_argument("--max-offset", type=int, default=3,
                    help="offset at or below this is not a trim (default 3)")
    a = ap.parse_args()

    W, M = a.workdir, a.module
    jpath = a.json or os.path.join(W, "%s_Viral_PSSM.json" % M)
    if not os.path.exists(jpath):
        sys.exit("no module JSON at %s -- pass --json" % jpath)
    feats = json.load(open(jpath))[M]["features"]

    keys = [a.key] if a.key else sorted(feats)
    findings, tested = [], 0

    for key in keys:
        spec = feats.get(key) or {}
        ftype = spec.get("feature_type") or spec.get("type") or "CDS"
        if ftype == "mat_peptide":
            continue                      # a cut site is not expected to be Met

        cpath = os.path.join(W, "collections", M, key + ".fasta")
        if not os.path.exists(cpath):
            continue
        full = dict(read_fasta(cpath))

        installed = None
        if a.repo:
            pat = os.path.join(a.repo, "Viral-PSSMs", "%s.pssms" % M, key, "*.pssm")
            installed = {os.path.basename(p).split(".")[2] for p in glob.glob(pat)}
            if not installed:
                continue

        adir = os.path.join(W, "Alignments", M, key)
        for sub in ("corrected_alis", "reclustered_alis", "leftover_alis"):
            for f in sorted(glob.glob(os.path.join(adir, sub, "*.fa"))):
                tag = os.path.basename(f)[:-3]
                if installed is not None and tag not in installed:
                    continue
                rows = read_fasta(f)
                if not rows:
                    continue
                tested += 1
                offs = []
                for hid, aln in rows:
                    frag = aln.replace("-", "")[:30]
                    fu = full.get(hid)
                    if fu and frag:
                        i = fu.find(frag)
                        if i >= 0:
                            offs.append(i)
                res, frac = modal_first_residue([x for _, x in rows])
                med = statistics.median(offs) if offs else None

                if med is not None and med > a.max_offset:
                    findings.append(("TRIMMED", key, sub + "/" + tag, len(rows), res,
                                     "alignment starts %d residues into its members"
                                     % med))
                elif res != "M":
                    findings.append(("NON_MET", key, sub + "/" + tag, len(rows), res,
                                     "members do not begin with Met (modal first "
                                     "residue %s, %.0f%%)" % (res, 100 * frac)))

    print("  %s: %d %s alignment(s) tested" % (M, tested,
          "installed" if a.repo else "on-disk"))
    if not findings:
        print("  every CDS profile begins at its protein's start")
        return 0

    by = collections.Counter(f[0] for f in findings)
    print("  %d finding(s): %s" % (len(findings),
          ", ".join("%s %d" % (k, v) for k, v in sorted(by.items()))))
    print()
    print("  %-9s %-10s %-26s %6s %5s  %s"
          % ("kind", "feature", "alignment", "nseqs", "first", "detail"))
    for kind, key, tag, n, res, detail in findings:
        print("  %-9s %-10s %-26s %6d %5s  %s" % (kind, key, tag, n, res, detail))
    print()
    print("  TRIMMED  -> split the cluster on its N-terminus; trimming further "
          "deletes real protein")
    print("  NON_MET  -> drop collection sequences that another sequence extends, "
          "then rebuild")
    return 1


sys.exit(main())
