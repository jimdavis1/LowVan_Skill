#!/usr/bin/env python3
"""Validate what the annotator produced for one genome, before you commit.

`evaluate_module.py` asks whether the calls line up with GenBank. This asks a
different and more basic question: is the output *internally* sound? It needs no
reference annotation, so it works on a genome nobody has annotated yet — which
is the case you care about when deciding whether a module is fit to ship.

Per called feature it checks the protein the annotator actually emitted:

  NO_MET        does not begin with M. Sometimes real, usually a start-site miss.
                Skipped for a mature peptide whose N-terminus is a cleavage site
                (upstream_ext 0) -- there is no start codon there to find.
  INTERNAL_STOP a `*` inside the protein. A defect UNLESS the feature declares
                "internal_stop" : 1, which says a single stop is expected there
                and read through -- alphavirus nsP3 carries an opal six codons
                before its C-terminus. More stops than declared is still a
                defect, and is reported as INTERNAL_STOP_EXCESS.
  LENGTH        outside the min_len/max_len the JSON declares for that feature.
  BOUNDS        coordinates reversed, or running off the end of the contig.

And across the genome:

  DUPLICATE     two features on identical coordinates. Always a defect -- some
                cluster is binned under the wrong feature. See qc_cross_feature.py.
  MISSING       a feature the module declares with copy_num that was not called.
  DERIVED_GAP   a CDS was called but its mature peptide / signal peptide was not.
                Usually means the derived feature covers fewer clusters than its
                parent, so some genomes get the CDS and nothing else.
  UNROUTED      no module matched; the genome was rejected before any PSSM ran.

    python3 validate_calls.py --repo ~/Viral_Annotation --fasta genome.fna
    python3 validate_calls.py --repo ~/Viral_Annotation --prefix out   # reuse a run

Exit is non-zero if any hard defect (duplicate, internal stop, bounds) is found.
"""

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict

HARD = {"DUPLICATE", "INTERNAL_STOP", "INTERNAL_STOP_EXCESS", "BOUNDS", "UNROUTED"}

# How many internal stops a feature declaring internal_stop may carry. Must
# match the annotator's -mis default, or this reports as broken exactly the
# features the annotator deliberately kept.
MAX_INTERNAL_STOPS = 1


def read_fasta(path):
    hdr, seq = None, []
    for line in open(path, errors="replace"):
        if line.startswith(">"):
            if hdr is not None:
                yield hdr, "".join(seq)
            hdr, seq = line[1:].rstrip("\n"), []
        else:
            seq.append(line.strip())
    if hdr is not None:
        yield hdr, "".join(seq)


def run_annotator(repo, fasta, prefix, threads):
    env = dict(os.environ, LOWVAN_DATA_DIR=repo)
    r = subprocess.run(
        ["perl", os.path.join(repo, "annotate_by_viral_pssm.pl"),
         "-i", fasta, "-p", prefix, "-threads", str(threads)],
        capture_output=True, text=True, env=env)
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True, help="Viral_Annotation checkout")
    ap.add_argument("--fasta", help="genome to annotate")
    ap.add_argument("--prefix", help="reuse an existing run's output prefix")
    ap.add_argument("--json", help="Viral_PSSM.json (default: <repo>/Viral_PSSM.json)")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    repo = os.path.abspath(os.path.expanduser(args.repo))
    jf = args.json or os.path.join(repo, "Viral_PSSM.json")
    with open(jf) as fh:
        mod = json.load(fh)

    tmp = None
    prefix = args.prefix
    if not prefix:
        if not args.fasta:
            sys.exit("give --fasta or --prefix")
        tmp = tempfile.mkdtemp(prefix="vcalls.")
        prefix = os.path.join(tmp, "run")
        r = run_annotator(repo, os.path.abspath(args.fasta), prefix, args.threads)
        if not os.path.exists(prefix + ".feature.tbl"):
            print("annotator produced no feature table.")
            print((r.stderr or r.stdout or "")[-900:])
            return 1

    tbl, faa = prefix + ".feature.tbl", prefix + ".faa"
    if not os.path.exists(tbl):
        sys.exit("no %s" % tbl)

    rows = []
    for line in open(tbl):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        # Column 7 is the GENE SYMBOL, not the JSON feature key -- "G_mature"
        # where the JSON says "G_MAT". The PSSM filename in the last-but-one
        # column carries the real key as <Module>.<KEY>.<cluster>.pssm, so take
        # it from there and fall back to the symbol only when there is no PSSM
        # (a "special" feature called by an external program).
        pssm = f[16] if len(f) > 16 else ""
        m = re.match(r"^.+?\.(.+)\.[^.]+\.pssm$", os.path.basename(pssm)) if pssm else None
        # A no_features_called row carries the routing decision (module, rep
        # contig, bitscore) with EMPTY coordinates, so that a genome which
        # routed but called nothing still leaves a record. It is not a feature.
        # Without this guard the script dies with ValueError on int('') for
        # exactly the genomes most worth validating -- the ones that called
        # nothing. check_cleaved_ends.py already guards the same row shape.
        try:
            int(f[7]); int(f[8])
        except (ValueError, IndexError):
            continue
        rows.append({"contig": f[2], "type": f[4], "fid": f[5],
                     "key": m.group(1) if m else f[6], "symbol": f[6],
                     "start": int(f[7]), "end": int(f[8]), "strand": f[9],
                     "module": f[11], "anno": f[-1]})
    prot = {}
    if os.path.exists(faa):
        for hdr, s in read_fasta(faa):
            prot[hdr.split()[0]] = s

    if not rows:
        print("UNROUTED -- no features called. The genome matched no representative")
        print("contig above the bitscore floor, so no PSSM ever ran.")
        print("Run check_rep_contigs.py; this needs another rep contig, not a better PSSM.")
        return 1

    module = rows[0]["module"]
    declared = mod.get(module, {}).get("features", {})
    contig_len = {}
    if args.fasta:
        for hdr, s in read_fasta(args.fasta):
            contig_len[hdr.split()[0]] = len(s)

    findings = []

    # duplicates: identical coordinates under two keys
    bycoord = defaultdict(list)
    for r in rows:
        if r["type"] == "CDS":
            bycoord[(r["contig"], r["start"], r["end"])].append(r["key"])
    for coord, keys in bycoord.items():
        if len(keys) > 1:
            findings.append(("DUPLICATE", "/".join(keys),
                             "%s:%d-%d called under %d keys" % (coord[0], coord[1], coord[2], len(keys))))

    for r in rows:
        p = prot.get(r["fid"], "")
        spec = declared.get(r["key"], {})
        if r["end"] < r["start"]:
            findings.append(("BOUNDS", r["key"], "end %d < start %d" % (r["end"], r["start"])))
        cl = contig_len.get(r["contig"])
        if cl and r["end"] > cl:
            findings.append(("BOUNDS", r["key"], "end %d beyond contig length %d" % (r["end"], cl)))
        if not p:
            continue
        body = p[:-1] if p.endswith("*") else p
        # A cleaved N-terminus is not a start codon, so Met is not expected
        # there. upstream_ext == 0 on a mat_peptide is exactly that declaration.
        n_is_cut = (r["type"] == "mat_peptide" and spec.get("upstream_ext") == 0)
        if not body.startswith("M") and not n_is_cut:
            findings.append(("NO_MET", r["key"], "starts %s..., %d aa" % (body[:4], len(body))))
        if "*" in body:
            n = body.count("*")
            at = ",".join(str(i + 1) for i in range(len(body)) if body[i] == "*")[:40]
            allowed = spec.get("internal_stop")
            if not allowed:
                findings.append(("INTERNAL_STOP", r["key"],
                                 "%d internal stop(s) at %s" % (n, at)))
            elif n > MAX_INTERNAL_STOPS:
                findings.append(("INTERNAL_STOP_EXCESS", r["key"],
                                 "%d internal stop(s) at %s; internal_stop permits %d"
                                 % (n, at, MAX_INTERNAL_STOPS)))
            else:
                findings.append(("READTHROUGH", r["key"],
                                 "%d expected internal stop at %s (internal_stop declared)"
                                 % (n, at)))
        lo, hi = spec.get("min_len"), spec.get("max_len")
        if lo and hi and not (lo <= len(body) <= hi):
            findings.append(("LENGTH", r["key"],
                             "%d aa, JSON declares %d-%d" % (len(body), lo, hi)))

    called = {r["key"] for r in rows}
    missing = [k for k, e in declared.items()
               if e.get("copy_num") and k not in called and not e.get("special")]

    # A derived feature is named <PARENT>_MAT / <PARENT>_SP by convention. If the
    # parent CDS was called and the derived peptide was not, the derived profiles
    # cover fewer clusters than the parent does -- worth knowing, because the
    # genome silently gets a CDS with no mature peptide.
    derived_gap = []
    for k, e in declared.items():
        if e.get("feature_type") != "mat_peptide" or k in called:
            continue
        parent = re.sub(r"_(MAT|SP|MATURE)$", "", k)
        if parent != k and parent in called:
            derived_gap.append((k, parent))

    print("module: %s      %d features called (%d CDS, %d mat_peptide)"
          % (module, len(rows),
             sum(1 for r in rows if r["type"] == "CDS"),
             sum(1 for r in rows if r["type"] == "mat_peptide")))
    #  Say out loud that special features were not evaluated. This tool runs
    #  annotate_by_viral_pssm.pl only; transcript_edit and splice features are
    #  called downstream by get_transcript_edited_features.pl and
    #  get_splice_variant_features.pl, which it does not invoke. Without this
    #  line the count above reads as though the module declared features it
    #  failed to call -- on Togaviridae the full pipeline calls 14 and this
    #  tool shows 13, and the missing one is simply out of its scope.
    specials = sorted(k for k, e in declared.items() if e.get("special"))
    if specials:
        print("  %d special feature(s) NOT evaluated here: %s"
              % (len(specials), ", ".join(specials)))
        print("  they are called by get_transcript_edited_features.pl / "
              "get_splice_variant_features.pl,")
        print("  which this tool does not run. Use the full pipeline to see them.")
    print()
    print("  %-12s %-9s %-11s %-8s %-13s %s"
          % ("type", "key", "symbol", "length", "coords", "protein starts"))
    for r in rows:
        p = prot.get(r["fid"], "")
        body = p[:-1] if p.endswith("*") else p
        print("  %-12s %-9s %-11s %5d aa  %5d-%-7d %s"
              % (r["type"], r["key"], r["symbol"], len(body),
                 r["start"], r["end"], body[:24]))

    if findings:
        hard = [f for f in findings if f[0] in HARD]
        print("\n%d finding(s), %d hard:\n" % (len(findings), len(hard)))
        for kind, key, detail in sorted(findings, key=lambda x: x[0] not in HARD):
            print("  %-14s %-16s %s" % (kind, key, detail))
        if any(f[0] == "DUPLICATE" for f in findings):
            print("\n  A DUPLICATE means a cluster is binned under the wrong feature and its")
            print("  PSSM out-scores the correct one at that locus. Run qc_cross_feature.py.")
    else:
        print("\n  no findings")

    if missing:
        print("\nMISSING (%d) -- declared with copy_num but not called:" % len(missing))
        print("  %s" % ", ".join(sorted(missing)))
        print("  Either the genome genuinely lacks them, or a bit_cutoff is too high.")

    if derived_gap:
        print("\nDERIVED_GAP (%d) -- parent CDS called, derived peptide not:" % len(derived_gap))
        for k, parent in sorted(derived_gap):
            print("  %-16s derives from %s, which WAS called" % (k, parent))
        print("  The derived profiles cover fewer clusters than the parent. Check whether")
        print("  this genome's parent cluster has a matching derived profile at all.")

    if tmp:
        print("\n  (run output in %s)" % tmp)
    return 1 if any(f[0] in HARD for f in findings) else 0


if __name__ == "__main__":
    sys.exit(main())
