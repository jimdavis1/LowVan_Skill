#!/usr/bin/env python3
"""Backfill a genome map for every module that has a coverage audit.

Runs entirely off a LOCAL copy of the runtime repo. The first attempt read
PSSMs and rep contigs straight off Box under six-way concurrency and failed
two ways at once: shutil.copy raised ETIMEDOUT mid-run and killed the pool,
and tblastn reported "Unsupported format for PSSM" on modules whose profiles
are perfectly good -- a short read, not a bad file.

For each module: annotate rep contigs, routing against that module's own
references only, keep the genome that calls the most distinct declared
features, and draw the map from its feature table. Maps are written to a
local staging dir; copy them to reports/ afterwards in one pass.
"""
import json, os, re, subprocess, sys, glob, shutil, traceback
from concurrent.futures import ThreadPoolExecutor

VA = os.environ.get("LOWVAN_DATA_DIR",
                    os.path.expanduser("~/lowvan_runtime/Viral_Annotation"))
KIT = os.environ["LOWVAN_KIT"]
REPORTS = KIT + "/reports"
GEN = KIT + "/skill/scripts/gen_genome_map.py"
WORK = os.path.expanduser("~/lowvan_runtime/genome_maps")
OUTDIR = WORK + "/_maps"
REFROOT = WORK + "/_refs"
TRIES = int(os.environ.get("GM_TRIES", "25"))

VJSON = json.load(open(VA + "/Viral_PSSM.json"))
AUDITS = {f[:-len("-coverage-audit.html")]
          for f in os.listdir(REPORTS) if f.endswith("-coverage-audit.html")}
HAVE = {f[:-len("-genome-map.html")]
        for f in os.listdir(REPORTS) if f.endswith("-genome-map.html")}


def audit_name(mod):
    for c in (mod.lower().replace("_", "-"), mod.lower(),
              mod.lower().replace("-", "_")):
        if c in AUDITS:
            return c
    return None


def run_one(mod):
    try:
        return _run(mod)
    except Exception:
        return (mod, None, "EXCEPTION " + traceback.format_exc().strip().splitlines()[-1])


def _run(mod):
    stem = audit_name(mod)
    if not stem:
        return (mod, None, "no coverage audit page")
    if stem in HAVE:
        return (mod, stem, "already has a map")
    want = set(VJSON.get(mod, {}).get("features", {}))
    if not want:
        return (mod, stem, "no features in Viral_PSSM.json")

    refdir = "%s/%s" % (REFROOT, mod)
    os.makedirs(refdir, exist_ok=True)
    #  "One of the complete, good-quality reference genomes used in the
    #  blastn" (JD, 17 Nov 2026). The rep contigs ARE those references, but
    #  they are not all whole genomes -- many are single deposited records.
    #  Longest first is the cheap proxy for completeness, and the loop keeps
    #  whichever actually calls the most declared features, so a long partial
    #  loses to a shorter complete one on the measurement rather than on the
    #  guess. Taking rep contig #1 regardless is what drew Amdoparvovirus
    #  with 1 of its 5 features.
    refs = glob.glob("%s/Viral-Rep-Contigs/%s.*.dna" % (VA, mod))
    refs.sort(key=lambda p: -os.path.getsize(p))
    for r in refs:
        t = os.path.join(refdir, os.path.basename(r))
        if not os.path.exists(t):
            shutil.copy(r, t)
    if not refs:
        return (mod, stem, "no rep contigs installed")

    d = "%s/%s" % (WORK, mod)
    os.makedirs(d, exist_ok=True)
    best, errs = None, set()
    for ref in refs[:TRIES]:
        n = re.search(r"\.(\d+)\.dna$", ref).group(1)
        fa, pre = "%s/%s.fa" % (d, n), "%s/%s" % (d, n)
        shutil.copy(os.path.join(refdir, os.path.basename(ref)), fa)
        try:
            p = subprocess.run(
                ["perl", VA + "/annotate_by_viral_pssm.pl", "-i", fa,
                 "-p", pre, "-threads", "2", "-c", refdir,
                 "-max", "500000", "-min", "200"],
                cwd=d, capture_output=True, timeout=2400, text=True)
            for ln in (p.stderr or "").splitlines():
                if "error" in ln.lower() or "Unsupported" in ln:
                    errs.add(ln.strip()[:90])
        except subprocess.TimeoutExpired:
            errs.add("annotator timeout")
            continue
        tbl = pre + ".feature.tbl"
        if not os.path.exists(tbl) or not os.path.getsize(tbl):
            continue
        keys, rows = set(), 0
        for line in open(tbl):
            q = line.rstrip("\n").split("\t")
            if len(q) < 18:
                continue
            rows += 1
            t = q[16][:-5] if q[16].endswith(".pssm") else q[16]
            if t.startswith(mod + "."):
                t = t[len(mod) + 1:]
            keys.add(t.split(".")[0])
        hit = keys & want
        score = (len(hit), -(rows - len(hit)))
        if best is None or score > best[0]:
            best = (score, fa, tbl, len(hit))
        if len(hit) == len(want):
            break
    if not best:
        return (mod, stem, "NO CALLS" + (" | " + "; ".join(sorted(errs)[:2]) if errs else ""))

    os.makedirs(OUTDIR, exist_ok=True)
    out = "%s/%s-genome-map.html" % (OUTDIR, stem)
    r = subprocess.run(
        [sys.executable, GEN, "--tbl", best[2], "--fasta", best[1],
         "--module-json", VA + "/Viral_PSSM.json", "--module", mod,
         "--out", out], capture_output=True, text=True)
    if r.returncode:
        return (mod, stem, "gen failed: " + r.stderr.strip()[-160:])
    return (mod, stem, "map %d/%d features" % (best[3], len(want)))


def main():
    only = os.environ.get("GM_ONLY")
    mods = only.split(",") if only else sorted(VJSON)
    os.makedirs(WORK, exist_ok=True)
    with ThreadPoolExecutor(max_workers=int(os.environ.get("GM_JOBS", "6"))) as ex:
        for mod, stem, msg in ex.map(run_one, mods):
            print("%-34s %-32s %s" % (mod, stem or "-", msg), flush=True)


main()
