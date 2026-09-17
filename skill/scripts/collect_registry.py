#!/usr/bin/env python3
"""Measure how much of its own training data each feature's PSSMs actually call.

This is the "PSSM effectiveness" number. For every feature it searches that
feature's PSSMs against that feature's own collection and counts how many
sequences score at or above the JSON bit_cutoff. Self-recall is a floor, not a
success criterion -- a profile that cannot recover the sequences it was built
from will certainly not recover anything new -- but it is the one measurement
that needs no held-out data, and a low value points straight at a cluster that
should have been split.

Writes registry.json for gen_registry.py.

    python3 collect_registry.py --workdir . --out registry.json

Status assigned per feature:
    built           has PSSMs, has a collection
    derived         has PSSMs but no collection of its own (G_MAT, G_SP)
    special         called by an external program (transcript_edit / splice)
    too-few         collection exists but is smaller than --min-seqs
    no-collection   declared in the JSON with no sequences behind it at all
"""

import argparse
import glob
import json
import os
import re
import subprocess
import sys
import tempfile
from collections import OrderedDict


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


def collection_paths(workdir, module, feat):
    """Main collection plus any .outliers split, which is the same feature."""
    base = os.path.join(workdir, "collections", module)
    out = []
    for name in (feat + ".fasta", feat + ".outliers.fasta"):
        p = os.path.join(base, name)
        if os.path.exists(p):
            out.append(p)
    return out


def _has_departures(bp):
    """True only if BUILD_PARAMS records an actual departure from the defaults.

    The skill tells you to write a BUILD_PARAMS next to EVERY feature, and
    run_pipeline.sh does, with "departures   none".  Testing os.path.exists()
    therefore marked every feature in the module as non-default, which is the
    opposite of what the legend claims.
    """
    if not os.path.exists(bp):
        return False
    for line in open(bp, errors="replace"):
        if line.startswith("departures"):
            v = line.split(None, 1)[1].strip() if len(line.split(None, 1)) > 1 else ""
            return bool(v) and v.lower() not in ("none", "-", "n/a")
    return False


def self_recall(pssms, colls, bit_cutoff, tmp, threads):
    """(n_sequences, n_called) for one feature."""
    seqs = OrderedDict()
    for c in colls:
        for hdr, s in read_fasta(c):
            seqs[hdr.split()[0]] = s
    if not seqs:
        return 0, None
    db = os.path.join(tmp, "coll.faa")
    with open(db, "w") as fh:
        for i, (k, s) in enumerate(seqs.items()):
            fh.write(">s%d\n%s\n" % (i, s))
    keys = list(seqs)
    subprocess.run(["makeblastdb", "-in", db, "-dbtype", "prot",
                    "-out", os.path.join(tmp, "coll")], capture_output=True)
    hit = set()
    for p in pssms:
        r = subprocess.run(
            ["psiblast", "-in_pssm", p, "-db", os.path.join(tmp, "coll"),
             "-outfmt", "6 sseqid bitscore", "-evalue", "10",
             "-max_target_seqs", str(len(keys) + 10), "-num_threads", str(threads)],
            capture_output=True, text=True)
        for line in r.stdout.splitlines():
            sid, bits = line.split("\t")[:2]
            if float(bits) >= bit_cutoff:
                hit.add(sid)
    return len(keys), len(hit)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--json", help="module JSON (default: <workdir>/*_Viral_PSSM.json)")
    ap.add_argument("--out", default="registry.json")
    ap.add_argument("--min-seqs", type=int, default=3,
                    help="below this a collection cannot support a profile (default 3)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--no-recall", action="store_true",
                    help="skip the BLAST work; report structure only")
    args = ap.parse_args()

    workdir = os.path.abspath(args.workdir)
    jf = args.json
    if not jf:
        c = [p for p in glob.glob(os.path.join(workdir, "*_Viral_PSSM.json"))
             if os.path.basename(p) != "Viral_PSSM.json"]
        if not c:
            sys.exit("no module JSON found")
        jf = c[0]
    with open(jf) as fh:
        mod = json.load(fh, object_pairs_hook=OrderedDict)

    tmp = tempfile.mkdtemp(prefix="registry.")
    out = OrderedDict()
    for module, block in mod.items():
        feats = []
        for key, ent in block.get("features", {}).items():
            pdir = os.path.join(workdir, "Alignments", module, key, "pssms")
            pssms = sorted(glob.glob(os.path.join(pdir, "*.pssm")))
            colls = collection_paths(workdir, module, key)
            nseq = sum(1 for c in colls for _ in read_fasta(c))

            if ent.get("special"):
                status = "special"
            elif pssms and not colls:
                status = "derived"
            elif pssms:
                status = "built"
            elif not colls or nseq == 0:
                status = "no-collection"
            else:
                status = "too-few"

            called = None
            if pssms and colls and not args.no_recall:
                nseq, called = self_recall(pssms, colls,
                                           float(ent.get("bit_cutoff") or 0),
                                           tmp, args.threads)

            bp = os.path.join(workdir, "Alignments", module, key, "BUILD_PARAMS")
            feats.append(OrderedDict([
                ("key", key),
                ("anno", ent.get("anno", "")),
                ("gene", ent.get("gene_symbol", "")),
                ("ftype", ent.get("feature_type", "")),
                ("segment", ent.get("segment", "")),
                ("bit", ent.get("bit_cutoff")),
                ("minl", ent.get("min_len")),
                ("maxl", ent.get("max_len")),
                ("pmid", ent.get("PMID")),
                ("n", nseq),
                ("called", called),
                ("pssms", [os.path.basename(p) for p in pssms]),
                ("nondefault", _has_departures(bp)),
                ("status", status),
            ]))
            print("  %-22s %-14s %-14s n=%-5s called=%-5s %d pssm"
                  % (module, key, status, nseq,
                     "-" if called is None else called, len(pssms)))
        out[module] = OrderedDict([("features", feats)])

    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)
    tot = sum(len(v["features"]) for v in out.values())
    npss = sum(len(f["pssms"]) for v in out.values() for f in v["features"])
    print("\n  %d modules, %d features, %d PSSMs -> %s" % (len(out), tot, npss, args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
