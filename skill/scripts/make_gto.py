#!/usr/bin/env python3
"""Create Genome Type Objects for evaluation, via rast-create-genome.

The GTO path is the one that matters for quality scoring: a flat feature table
cannot represent special features, and it miscounts fragments, so a genome whose
L is split across two contigs scores as a short L rather than a fragmented one.
Scoring a module on feature tables understates it.

Use rast-create-genome. Do not hand-write the JSON: it issues the genome id
from the ID server, and `new_feature_id` numbers every feature the annotator
adds off that id. A GTO carrying a reused or invented id produces feature ids
that collide with a real genome's.

    python3 make_gto.py --fasta-dir gto_in --metadata exemplars.metadata \
            --out gto --jobs 4

`rast-create-genome` ships in the BV-BRC dev kit. On macOS the desktop bundle
has it at /Applications/BV-BRC.app/deployment/, though that release ships the
plbin script without its bin wrapper; --tool points at either.

NOTE: each call registers a new genome id with the ID server. That is a real,
outward-facing side effect, and it is the point -- but it means a run over a
whole taxon registers one id per exemplar.

Most of each call is the round trip to the ID server, not local work, so it
parallelises in principle. In practice the service rate-limits on concurrency
and answers **403 Forbidden**, not 429, so the failure reads like an expired
token and is not one: one Betaflexiviridae_MP run lost 89 of 120 genomes at
--jobs 4 with 4,200 hours left on the token. 403s are now retried with
exponential backoff and the default is 4; --jobs 1 avoids them entirely when
a run has to be reliable.
"""
import argparse, os, random, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

WRAPPER = """#!/bin/bash
dir=`dirname "$0"`
dir=`cd "$dir/../.."; pwd`
export KB_TOP="$dir/deployment"
export KB_RUNTIME="$dir/runtime"
export PATH="$KB_RUNTIME/bin:$KB_TOP/bin:$PATH"
export PERL5LIB="$KB_TOP/lib"
exec "$KB_RUNTIME/bin/perl" "$KB_TOP/plbin/rast-create-genome.pl" "$@"
"""


def resolve(tool):
    """Find rast-create-genome.

    Preference order: an explicit --tool, then the copy vendored in this
    repository, then one already on PATH, then a BV-BRC desktop bundle. The
    vendored copy is tried before PATH so a clone works with no BV-BRC install
    at all -- which is the point of shipping it.
    """
    if tool:
        return tool
    #  vendor/bv-brc/bin/, relative to this script: skill/scripts/ -> repo root
    here = os.path.dirname(os.path.abspath(__file__))
    for up in (3, 2, 4):
        root = os.path.abspath(os.path.join(here, *([".."] * up)))
        b = os.path.join(root, "vendor", "bv-brc", "bin", "rast-create-genome")
        if os.path.exists(b):
            return b
    from shutil import which
    w = which("rast-create-genome")
    if w:
        return w
    for base in ("/Applications/BV-BRC.app", os.path.expanduser("~/BV-BRC.app")):
        b = os.path.join(base, "deployment", "bin", "rast-create-genome")
        if os.path.exists(b):
            return b
        if os.path.exists(os.path.join(base, "deployment", "plbin",
                                       "rast-create-genome.pl")):
            #  the 1.046 bundle ships the plbin script but not its wrapper
            with open(b, "w") as fh:
                fh.write(WRAPPER)
            os.chmod(b, 0o755)
            print("  created missing wrapper %s" % b)
            return b
    return None


def preflight(tool):
    """Check the perl this tool will use can actually load what it needs.

    GenomeTypeObject loads its UUID module inside an eval, so a missing one is
    silent until create_uuid() is called and the run dies with "No UUID
    generator found" partway through. Fail here instead, with the fix.
    """
    import subprocess
    perl = os.environ.get("LOWVAN_PERL", "perl")
    root = os.path.abspath(os.path.join(os.path.dirname(tool), ".."))
    env = dict(os.environ, PERL5LIB=os.path.join(root, "lib")
               + (":" + os.environ["PERL5LIB"] if os.environ.get("PERL5LIB") else ""))
    missing = []
    for mod, why in (("File::Slurp", "required by GenomeTypeObject and CmdHelper"),
                     ("JSON::XS", "required by GenomeTypeObject")):
        if subprocess.run([perl, "-M" + mod, "-e", "1"],
                          capture_output=True, env=env).returncode:
            missing.append((mod, why))
    if all(subprocess.run([perl, "-M" + m, "-e", "1"],
                          capture_output=True, env=env).returncode
           for m in ("Data::UUID", "UUID")):
        missing.append(("Data::UUID", "GenomeTypeObject needs Data::UUID or UUID "
                                      "to mint feature ids"))
    if missing:
        print("ERROR: %s cannot load:" % perl, file=sys.stderr)
        for m, why in missing:
            print("   %-22s %s" % (m, why), file=sys.stderr)
        print("\n   cpanm %s" % " ".join(m for m, _ in missing), file=sys.stderr)
        print("   (or set LOWVAN_PERL to a perl that has them)", file=sys.stderr)
        return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta-dir", required=True)
    ap.add_argument("--metadata", help="headerless TSV: file, name, taxon_id")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tool", help="path to rast-create-genome")
    ap.add_argument("--taxon-default", default="11270")
    ap.add_argument("--genetic-code", type=int, default=1)
    ap.add_argument("--domain", default="Viruses")
    ap.add_argument("--jobs", type=int, default=4,
                    help="parallel rast-create-genome calls. The BV-BRC "
                         "service rate-limits above a handful and returns 403; "
                         "403s are retried with backoff, but --jobs 1 avoids "
                         "them entirely when a run has to be reliable.")
    args = ap.parse_args()

    tool = resolve(args.tool)
    if not tool:
        print("ERROR: rast-create-genome not found. Install the BV-BRC dev kit,\n"
              "       or pass --tool. Do not substitute hand-written JSON: the\n"
              "       genome id must be issued, not invented.", file=sys.stderr)
        return 1
    print("  tool: %s" % tool)
    if not preflight(tool):
        return 1

    meta = {}
    if args.metadata and os.path.exists(args.metadata):
        for line in open(args.metadata, errors="replace"):
            c = line.rstrip("\n").split("\t")
            if len(c) >= 3:
                meta[os.path.splitext(c[0])[0]] = (c[1].strip(), c[2].strip())

    os.makedirs(args.out, exist_ok=True)
    jobs = []
    for fn in sorted(os.listdir(args.fasta_dir)):
        if not fn.endswith((".fna", ".fasta", ".fa")):
            continue
        gid = os.path.splitext(fn)[0]
        out = os.path.join(args.out, gid + ".gto")
        if os.path.exists(out) and os.path.getsize(out) > 0:
            continue
        name, tax = meta.get(gid, (None, None))
        if not tax:
            head = gid.split(".")[0]
            tax = head if head.isdigit() else args.taxon_default
        if not name:
            name = gid.replace("MERGED_", "").replace("_", " ")
        jobs.append((os.path.join(args.fasta_dir, fn), out, name, tax, gid))
    print("  %d genome(s) to create" % len(jobs))

    def one(j):
        fa, out, name, tax, gid = j
        #  rast-create-genome calls the BV-BRC service, which rate-limits on
        #  concurrency and answers "403 Forbidden" rather than 429. At
        #  --jobs 4 that cost 89 of 120 genomes in one Betaflexiviridae_MP
        #  run while the token had 4,200 hours left on it, so the failure
        #  reads like an auth problem and is not one. Serial retries with
        #  backoff clear it.
        last = ""
        for attempt in range(5):
            r = subprocess.run([tool, "--scientific-name", name,
                                "--domain", args.domain,
                                "--genetic-code", str(args.genetic_code),
                                "--ncbi-taxonomy-id", str(tax),
                                "--source", "BV-BRC", "--source-id", gid,
                                "--contigs", fa, "-o", out],
                               capture_output=True, text=True)
            if r.returncode == 0 and os.path.exists(out):
                return None
            last = (r.stderr or r.stdout or "")[-200:]
            if "403" not in last and "Forbidden" not in last:
                break
            time.sleep(2 ** attempt + random.random())
        return (gid, last)

    bad = []
    done = 0
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for res in pool.map(one, jobs):
            done += 1
            if res:
                bad.append(res)
            if done % 100 == 0:
                print("    %d/%d" % (done, len(jobs)), flush=True)
    n = len([f for f in os.listdir(args.out) if f.endswith(".gto")])
    print("  %d GTO(s) in %s" % (n, args.out))
    if bad:
        print("  %d failed:" % len(bad))
        for g, e in bad[:5]:
            print("    %-24s %s" % (g, e.replace("\n", " ")[:110]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
