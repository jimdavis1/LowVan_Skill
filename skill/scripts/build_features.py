#!/usr/bin/env python3
"""Run fasta-cluster-pssm-2.pl over a module's collections.

The pipeline writes into a randomly-named directory and never tells you which
one, so this captures the set of directories before and after and renames the
new one to the feature key. It also records run.stdout, run.stderr and a
BUILD_PARAMS file next to every feature, which is what the registry artifact
reads to mark a feature as non-default.
"""
import os, sys, subprocess, shutil, argparse, time, json

DEFAULTS = dict(m="5", f="0.33", n="0.75", c="0", fd="0.20", e="3", efo="0.15",
                mi="0.8", mc="0.8", p_nterm="65", n_nterm="3",
                nterm_eval_length="10")
MI_FLOOR = 0.6   # references/pipeline.md: never redefine "same protein" below this

def run_one(workdir, module, key, anno, params, extra_flags):
    cdir = os.path.join(workdir, "collections", module)
    adir = os.path.join(workdir, "Alignments", module)
    os.makedirs(adir, exist_ok=True)
    src = os.path.join(cdir, key + ".fasta")
    if not os.path.exists(src) or os.path.getsize(src) == 0:
        print("  %-14s no collection -- skipped" % key); return None

    assert float(params["mi"]) >= MI_FLOOR, "-mi %s breaches the %.1f floor" % (params["mi"], MI_FLOOR)
    dest = os.path.join(adir, key)
    if os.path.exists(dest):
        shutil.rmtree(dest)
    before = set(os.listdir(adir))

    cmd = ["fasta-cluster-pssm-2.pl", "-a", anno, "-p", "%s.%s" % (module, key),
           "-m", params["m"], "-f", params["f"], "-n", params["n"], "-c", params["c"],
           "-fd", params["fd"], "-e", params["e"], "-efo", params["efo"],
           "-mi", params["mi"], "-mc", params["mc"],
           "-p_nterm", params["p_nterm"], "-n_nterm", params["n_nterm"],
           "-nterm_eval_length", params["nterm_eval_length"],
           "-nterm-mmseq-id", "0.7"] + extra_flags
    t = time.time()
    with open(src) as fin:
        p = subprocess.run(cmd, stdin=fin, cwd=adir, capture_output=True, text=True)
    after = set(os.listdir(adir)) - before
    new = [d for d in after if os.path.isdir(os.path.join(adir, d))]
    if len(new) != 1:
        print("  %-14s FAILED -- %d new dirs %s" % (key, len(new), sorted(new)[:4]))
        sys.stderr.write(p.stderr[-2000:]); return None
    os.rename(os.path.join(adir, new[0]), dest)
    # the pipeline names its scratch files after the tmp dir; drop them
    for fn in os.listdir(dest):
        if fn.startswith(new[0]):
            os.remove(os.path.join(dest, fn))
    open(os.path.join(dest, "run.stdout"), "w").write(p.stdout)
    open(os.path.join(dest, "run.stderr"), "w").write(p.stderr)

    departures = {k: v for k, v in params.items() if DEFAULTS.get(k) != v}
    with open(os.path.join(dest, "BUILD_PARAMS"), "w") as f:
        f.write("feature: %s\nannotation: %s\ncommand: %s\n" % (key, anno, " ".join(cmd)))
        f.write("departures: %s\n" % (json.dumps(departures) if departures else "none"))
        if extra_flags: f.write("extra_flags: %s\n" % " ".join(extra_flags))

    npssm = len(os.listdir(os.path.join(dest, "pssms"))) if os.path.isdir(os.path.join(dest, "pssms")) else 0
    nclu  = len([x for x in os.listdir(os.path.join(dest, "corrected_alis"))
                 if x.endswith(".fa")]) if os.path.isdir(os.path.join(dest, "corrected_alis")) else 0
    lo = os.path.join(dest, "Leftover_Seqs.aa")
    nleft = sum(1 for l in open(lo) if l.startswith(">")) if os.path.exists(lo) else 0
    nin = sum(1 for l in open(src) if l.startswith(">"))
    print("  %-14s %5d seqs -> %2d alignments, %2d PSSMs, %4d leftover   (%.0fs)%s"
          % (key, nin, nclu, npssm, nleft, time.time() - t,
             "  [" + ",".join("%s=%s" % kv for kv in sorted(departures.items())) + "]" if departures else ""))
    return dict(key=key, n=nin, alis=nclu, pssms=npssm, leftover=nleft)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--module", required=True)
    ap.add_argument("--features", required=True,
                    help='JSON file: {"KEY": {"anno": "...", "params": {...}, "flags": [...]}, ...}')
    ap.add_argument("--only", default="", help="comma-separated subset")
    a = ap.parse_args()
    spec = json.load(open(a.features))
    only = set(x for x in a.only.split(",") if x)
    print("building %s" % a.module)
    res = []
    for key, d in spec.items():
        if only and key not in only: continue
        params = dict(DEFAULTS); params.update(d.get("params", {}))
        r = run_one(a.workdir, a.module, key, d["anno"], params, d.get("flags", []))
        if r: res.append(r)
    print("\n%d features built, %d PSSMs total" % (len(res), sum(r["pssms"] for r in res)))

main()
