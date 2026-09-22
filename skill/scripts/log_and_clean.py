#!/usr/bin/env python3
"""Log the genomes a module was tested on, then delete the bulk that produced them.

Annotating a taxon leaves hundreds of megabytes of GTOs behind. They are
regenerable from the contigs and the installed module, and every number the
audit quotes has already been extracted from them, so keeping them costs
storage and buys nothing. What must not be lost is *which genomes were
analysed* -- without that the measurements are unreproducible in principle,
not just expensive to redo.

So this writes a manifest first and deletes second, and refuses to delete
anything if the manifest cannot be written.

    python3 log_and_clean.py --workdir <dir> --module <M> --box <boxdir>
    python3 log_and_clean.py ... --write        # actually delete
"""
import argparse, glob, hashlib, json, os, shutil, subprocess, sys, time

#  Directories of annotator output. Regenerable, and large.
BULK = ["gto", "gto_out", "te_full", "ns1p_test", "fp_test", "fix_test",
        "coverage_eval", "gto_fix_out", "gto_rest_out", "_wd"]
#  Never touched: the module itself and anything measurements are read from.
KEEP = ["collections", "collections_preX", "Alignments", "Rep-Contigs",
        "Transcript-Editing", "artifacts", "measurements", "Contigs"]


def du(p):
    try:
        return int(subprocess.run(["du", "-sk", p], capture_output=True,
                                  text=True).stdout.split()[0]) * 1024
    except Exception:
        return 0


SUF = (".ann.gto", ".qual.gto", ".te.gto", ".sv.gto", ".feature.tbl", ".gto")


def genomes_in(d):
    """Genome ids represented anywhere under a directory.

    Walks rather than globbing one level: a test panel keeps its GTOs in
    subdirectories (ns1p_test/out, ns1p_test/te), so a flat glob reported
    0 genomes for 36 MB of output and would have deleted it unlogged.
    """
    ids = set()
    for root, _dirs, files in os.walk(d):
        for b in files:
            for suf in SUF:
                if b.endswith(suf):
                    ids.add(b[:-len(suf)]); break
    return ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--module", required=True)
    ap.add_argument("--box", help="where the manifest is written; defaults to workdir")
    ap.add_argument("--note", default="", help="one line on what this run measured")
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()
    W = os.path.abspath(a.workdir)
    box = os.path.abspath(a.box) if a.box else W

    present = [d for d in BULK if os.path.isdir(os.path.join(W, d))]
    if not present:
        print("nothing to clean in %s" % W)
        return 0

    manifest = {
        "module": a.module,
        "written": time.strftime("%Y-%m-%d %H:%M:%S"),
        "note": a.note,
        "workdir": W,
        "why": ("GTO output is regenerable from Contigs/ plus the installed module; "
                "every figure the audit quotes was extracted before deletion. This "
                "manifest records which genomes were analysed so the measurements "
                "remain attributable."),
        "regenerate": [
            "python3 $LOWVAN_KIT/skill/scripts/make_gto.py --fasta-dir Contigs "
            "--metadata <meta>.tsv --out gto --jobs 6",
            "python3 $LOWVAN_KIT/skill/scripts/run_gto_eval.py --gto-dir gto "
            "--repo $LOWVAN_DATA_DIR --out gto_out --report quality.tsv --jobs 6 "
            "--perl /Applications/BV-BRC.app/runtime/bin/perl",
        ],
        "panels": {},
        "kept": [],
        "deleted": [],
    }

    allg = set()
    for d in present:
        p = os.path.join(W, d)
        g = genomes_in(p)
        allg |= g
        manifest["panels"][d] = {"genomes": len(g), "bytes": du(p)}
        if g:
            manifest["panels"][d]["ids"] = sorted(g)

    #  what survives, and a checksum of each so a later reader can tell whether
    #  the thing the manifest describes is the thing still on disk
    for k in KEEP:
        p = os.path.join(W, k)
        if os.path.isdir(p):
            manifest["kept"].append({"path": k, "bytes": du(p)})
    for f in sorted(glob.glob(os.path.join(W, "*.tsv")) +
                    glob.glob(os.path.join(W, "*.json"))):
        if os.path.getsize(f) < 50_000_000:
            h = hashlib.sha256(open(f, "rb").read()).hexdigest()[:16]
            manifest["kept"].append({"path": os.path.basename(f),
                                     "bytes": os.path.getsize(f), "sha256_16": h})

    manifest["genomes_analysed"] = len(allg)
    total = sum(v["bytes"] for v in manifest["panels"].values())
    manifest["bytes_to_delete"] = total

    os.makedirs(box, exist_ok=True)
    mpath = os.path.join(box, "ANALYSED_GENOMES.%s.json" % a.module)
    with open(mpath, "w") as fh:
        json.dump(manifest, fh, indent=1)
    #  a flat list too, because that is what anyone actually greps
    lpath = os.path.join(box, "ANALYSED_GENOMES.%s.txt" % a.module)
    with open(lpath, "w") as fh:
        for g in sorted(allg):
            fh.write(g + "\n")
    print("manifest: %s" % mpath)
    print("          %s  (%d genome ids)" % (lpath, len(allg)))

    print("\n%-16s %10s %s" % ("directory", "MB", "genomes"))
    for d in present:
        v = manifest["panels"][d]
        print("%-16s %10.1f %d" % (d, v["bytes"] / 1e6, v["genomes"]))
    print("%-16s %10.1f  <- reclaimable" % ("TOTAL", total / 1e6))

    if not a.write:
        print("\nnothing deleted (rerun with --write)")
        return 0
    if not os.path.exists(mpath):
        sys.exit("manifest was not written; refusing to delete")
    for d in present:
        shutil.rmtree(os.path.join(W, d), ignore_errors=True)
        manifest["deleted"].append(d)
    json.dump(manifest, open(mpath, "w"), indent=1)
    print("\ndeleted %d directories, reclaimed %.1f MB" % (len(present), total / 1e6))
    return 0


sys.exit(main())
