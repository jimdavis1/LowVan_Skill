#!/usr/bin/env python3
"""Run fasta-cluster-pssm-2.pl over a module's collections.

The pipeline writes into a randomly-named directory and never tells you which
one, so this captures the set of directories before and after and renames the
new one to the feature key. It also records run.stdout, run.stderr and a
BUILD_PARAMS file next to every feature, which is what the registry artifact
reads to mark a feature as non-default.
"""
import os
import glob, sys, subprocess, shutil, argparse, time, json

DEFAULTS = dict(m="5", f="0.33", n="0.75", c="0", fd="0.20", e="3", efo="0.15",
                mi="0.8", mc="0.8", p_nterm="65", n_nterm="3",
                nterm_eval_length="10")
MI_FLOOR = 0.6   # references/pipeline.md: never redefine "same protein" below this

def _normalise(given):
    """Accept {"-m": 2} or {"m": "5"} and return DEFAULTS-shaped string params.

    DEFAULTS is keyed WITHOUT the leading dash. A features.json written with
    dashes -- which is how the flags appear everywhere else, including in this
    script's own BUILD_PARAMS output -- therefore added unused keys and left
    the command on its defaults. Nothing failed: the run proceeded at -m 5 and
    -mi 0.8 while BUILD_PARAMS recorded `departures: {"-m": 2}`, so the file
    asserted a departure that had not happened. Avihepatovirus and Senecavirus
    were built at -m 5 for features holding 3 to 5 sequences, produced nothing,
    and the evidence said they had been built at -m 2.

    Unknown keys are fatal rather than ignored, for the same reason.
    """
    out = {}
    for k, v in (given or {}).items():
        kk = k[1:] if k.startswith("-") else k
        if kk not in DEFAULTS:
            raise SystemExit("build_features: unknown parameter %r "
                             "(known: %s)" % (k, ", ".join(sorted(DEFAULTS))))
        out[kk] = str(v)
    return out


def _write_params(dest, key, anno, params, cmd, extra_flags, note=None):
    """Record how a feature was built, next to the feature.

    build_leftover_pssms.py reads -mi from this file so the leftover pass
    cannot silently lower the identity floor, so it is load-bearing rather
    than documentation. It is written even when the output directory could
    not be identified, because the parameters are known before the run and a
    directory recovered by hand afterwards would otherwise have no record.
    """
    os.makedirs(dest, exist_ok=True)
    departures = {k: v for k, v in params.items() if DEFAULTS.get(k) != v}
    with open(os.path.join(dest, "BUILD_PARAMS"), "w") as f:
        f.write("feature: %s\nannotation: %s\ncommand: %s\n"
                % (key, anno, " ".join(cmd)))
        f.write("departures: %s\n"
                % (json.dumps(departures) if departures else "none"))
        if extra_flags:
            f.write("extra_flags: %s\n" % " ".join(extra_flags))
        if note:
            f.write("note: %s\n" % note)


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
        #  This destroys the feature's LEFTOVER profiles too, and nothing
        #  downstream notices: install_module.py checks that every declared
        #  feature has at least one PSSM, not that it still has as many as it
        #  had. Rebuilding three merged Carlavirus features took the module
        #  from 338 profiles to 234 and installed clean. Say how many are
        #  about to go, so a silent shrink is visible in the build log.
        gone = len(glob.glob(os.path.join(dest, "pssms", "*.pssm")))
        if gone:
            print("  %-14s rebuilding: %d existing profile(s) discarded; re-run "
                  "build_leftover_pssms.py for this key afterwards" % (key, gone))
        shutil.rmtree(dest)

    #  Clear any orphan output directory left by an earlier attempt at THIS
    #  feature. fasta-cluster-pssm-2.pl names its temp directory from the
    #  input, so re-running the same collection targets the same name, and it
    #  opens with `mkdir($tmp_dir) or die`. A leftover from a previous run
    #  therefore kills the next one before it starts, and the only symptom is
    #  this script reporting "FAILED -- no output directory carries <KEY>
    #  profiles (0 new dirs)" -- which reads as though nothing ran rather than
    #  as a collision. It cost several rebuild cycles across Coltivirus,
    #  Avian and Piscine Orthoreovirus before the pattern was clear: every
    #  time, deleting the stray hash directory by hand made the rebuild work.
    for d in list(os.listdir(adir)):
        p = os.path.join(adir, d)
        if d == key or not os.path.isdir(p):
            continue
        ps = glob.glob(os.path.join(p, "pssms", "*.pssm"))
        stem = os.path.basename(ps[0]).split(".") if ps else []
        owned = len(stem) > 2 and stem[1] == key
        #  An empty husk is a dead run of something; a husk is never useful.
        if owned or (not ps and not os.listdir(p)):
            shutil.rmtree(p, ignore_errors=True)
            print("  %-14s clearing stray output directory %s from an earlier "
                  "attempt" % (key, d))

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

    #  Identify the output by what is inside it, not by diffing the listing.
    #  fasta-cluster-pssm-2.pl can leave more than one directory behind for a
    #  single run, and the extra one sometimes appears late enough to be
    #  counted against the NEXT feature -- Alphaflexiviridae reported
    #  "TGB1 FAILED -- 3 new dirs ['93e9974de3', 'TGB2', 'fca9ad32d5']",
    #  naming an already-renamed feature as new. Four of seven features were
    #  declared failures while every one of them had built correctly.
    #
    #  The profiles carry the feature key in their filenames
    #  (<Module>.<KEY>.<n>.pssm), so that is the reliable signal. Duplicates
    #  are byte-identical; keep the most complete and drop the rest.
    def owner(d):
        ps = glob.glob(os.path.join(adir, d, "pssms", "*.pssm"))
        if not ps: return None, 0
        stem = os.path.basename(ps[0]).split(".")
        return (stem[1] if len(stem) > 2 else None), len(ps)

    mine = sorted(((d,) + owner(d) for d in new),
                  key=lambda x: -x[2])
    mine = [d for d, k, n in mine if k and key.startswith(k)]
    if not mine:
        #  Write the parameters anyway. They are known before the run and do
        #  not depend on finding the output, and a directory recovered by hand
        #  afterwards otherwise has no record of how it was built. That is not
        #  cosmetic: build_leftover_pssms.py reads -mi from here so the
        #  leftover pass cannot silently lower the identity floor, and with
        #  the file missing it refuses to run at all.
        _write_params(os.path.join(adir, key), key, anno, params, cmd,
                      extra_flags, note="output directory not identified; "
                      "parameters recorded from the invocation")
        print("  %-14s FAILED -- no output directory carries %s profiles "
              "(%d new dirs)" % (key, key, len(new)))
        sys.stderr.write(p.stderr[-2000:]); return None
    keep, dupes = mine[0], mine[1:]
    for d in dupes:
        shutil.rmtree(os.path.join(adir, d), ignore_errors=True)
    if dupes:
        print("  %-14s (dropped %d duplicate output dir(s))" % (key, len(dupes)))
    new = [keep]
    os.rename(os.path.join(adir, new[0]), dest)
    # the pipeline names its scratch files after the tmp dir; drop them
    for fn in os.listdir(dest):
        if fn.startswith(new[0]):
            os.remove(os.path.join(dest, fn))
    open(os.path.join(dest, "run.stdout"), "w").write(p.stdout)
    open(os.path.join(dest, "run.stderr"), "w").write(p.stderr)

    _write_params(dest, key, anno, params, cmd, extra_flags)

    npssm = len(os.listdir(os.path.join(dest, "pssms"))) if os.path.isdir(os.path.join(dest, "pssms")) else 0
    nclu  = len([x for x in os.listdir(os.path.join(dest, "corrected_alis"))
                 if x.endswith(".fa")]) if os.path.isdir(os.path.join(dest, "corrected_alis")) else 0
    lo = os.path.join(dest, "Leftover_Seqs.aa")
    nleft = sum(1 for l in open(lo) if l.startswith(">")) if os.path.exists(lo) else 0
    nin = sum(1 for l in open(src) if l.startswith(">"))
    #  recomputed here rather than carried from the writer, which now lives in
    #  _write_params -- the refactor left this line reading a variable that no
    #  longer existed in scope, and it only fires after a successful build so
    #  nothing caught it until a seven-feature rebuild reached the end.
    departures = {k: v for k, v in params.items() if DEFAULTS.get(k) != v}
    print("  %-14s %5d seqs -> %2d alignments, %2d PSSMs, %4d leftover   (%.0fs)%s"
          % (key, nin, nclu, npssm, nleft, time.time() - t,
             "  [" + ",".join("%s=%s" % kv for kv in sorted(departures.items())) + "]" if departures else ""))
    return dict(key=key, n=nin, alis=nclu, pssms=npssm, leftover=nleft)

def require_length_audit(workdir, module, force=False):
    """Refuse to build profiles over collections nobody has length-audited.

    THE PROFILE IS WHERE CONTAMINATION BECOMES PERMANENT. Once a mixed
    collection has been turned into a PSSM, every later check agrees with
    it: the profile recalls its own training data, the call rate is fine or
    better, and the annotation string is correct. The defect surfaces only as
    duplicate calls, which look like a cutoff problem -- so it gets chased
    with coverage runs, which is the most expensive place in the pipeline to
    discover it and the least likely place to understand it.

    Adenoviridae paid that twice. Mastadenovirus FIBER/FIBER1/FIBER2 came out
    95-99% identical to each other because "fiber-1" is positional and "long
    fiber" descriptive and the two schemes are INVERTED between the genera;
    UXP held 1,206 copies of a 55 aa exon beside 65 copies of the real 217 aa
    spliced protein (PMID 17881437). Neither was visible until panels had
    been run, and the curator's instruction, 5 Nov 2026, was to stop finding
    them there:

        "You need to decontaminate the pssms using the literature as your
         guide earlier in the pipeline before we start evaluating coverage
         etc. This is costing too much time."

    So the gate is here, at the build, not only at install. Run

        audit_collection_lengths.py --workdir <build> --module <M> \\
            --refs length_refs_<M>.json --accept accept_<M>.tsv \\
            --write-md <build>/LENGTH_AUDIT.<M>.md

    and resolve or sign off every finding first.
    """
    cands = [os.path.join(workdir, "LENGTH_AUDIT.%s.md" % module),
             os.path.join(workdir, module, "LENGTH_AUDIT.md"),
             os.path.join(workdir, "LENGTH_AUDIT.md")]
    path = next((p for p in cands if os.path.exists(p)), None)
    why = None
    if path is None:
        why = ("no length audit found. Looked for:\n     %s"
               % "\n     ".join(cands))
    else:
        import re as _re
        text = open(path, errors="replace").read()
        m = _re.search(r"\*\*Unresolved:\s*(\d+)", text)
        if not m:
            why = "%s does not state an unresolved count" % path
        elif int(m.group(1)) != 0:
            why = "%s reports %s unresolved finding(s)" % (path, m.group(1))
        else:
            #  A clean report that predates the collections describes
            #  something else. Staleness is the failure mode a gate like this
            #  invites, so it is checked rather than assumed.
            newest = 0
            for d in (os.path.join(workdir, "collections", module),
                      os.path.join(workdir, module, "collections")):
                for f in glob.glob(os.path.join(d, "*.fasta")):
                    newest = max(newest, os.path.getmtime(f))
            if newest and os.path.getmtime(path) < newest:
                why = ("%s is OLDER than the collections it describes"
                       % path)
    if why is None:
        print("length audit: clean (%s)" % path)
        return
    if force:
        print("WARNING: building over an uncurated collection set -- %s" % why)
        return
    sys.exit(
        "REFUSED: %s\n"
        "  A PSSM built from a mixed collection passes every later check --\n"
        "  it recalls its own training data, the call rate looks right, and\n"
        "  the annotation string is correct. The only symptom is duplicate\n"
        "  calls, and finding that in a coverage run costs a rebuild of\n"
        "  everything downstream.\n"
        "  Run audit_collection_lengths.py with --refs (literature ranges,\n"
        "  each with a PMID) and resolve or --accept every finding, then\n"
        "  build. --force-uncurated overrides, for an experiment only."
        % why)


def require_bin_audit(workdir, module, force=False):
    """Refuse to build profiles until somebody has asked whether each
    collection is ONE protein, and its OWN protein.

    The length audit asks whether a collection's members are the right size;
    audit_contamination.py moves individual mislabelled sequences between
    collections. Neither can see a collection that should not exist. The
    accessory-discovery pass declares every unassigned cluster that reaches
    the member floor as a feature, and it cannot tell a pile of partials, or a
    second copy of a named protein, from a new one.

    Mastadenovirus, 5 Nov 2026: of 115 collections, 3 were a protein that
    already had a bin (UNC16 was POL at 99.8% identity), 15 were partials
    (UNC23 was 100% identical to E1B55K over 17% of it; E3_CR1D, a NAMED
    feature, was a 46 aa fragment of E3_CR1B), and 7 matched a named feature
    at 29-37% across their whole length. 293 sequences, and they caused most
    of the module's 458 duplicate calls.

    This cannot be fixed downstream. A cutoff separates a profile's true loci
    from its false ones; these profiles score HIGHER on the foreign locus than
    on their own, because the foreign protein is what they were built from.
    UNC34 scored 80 bits on E1B19K's locus and 33 on its own. The symptom --
    duplicate calls -- looks exactly like a cutoff problem, so it gets chased
    through coverage runs instead.

        audit_bins.py --workdir <build> --module <M> \\
            --protect K22:K33 --protect UXP:U_EXON \\
            --accept accept_bins_<M>.tsv --write-md <build>/BIN_AUDIT.<M>.md
    """
    cands = [os.path.join(workdir, "BIN_AUDIT.%s.md" % module),
             os.path.join(workdir, module, "BIN_AUDIT.md"),
             os.path.join(workdir, "BIN_AUDIT.md")]
    path = next((p for p in cands if os.path.exists(p)), None)
    why = None
    if path is None:
        why = ("no bin audit found. Looked for:\n     %s"
               % "\n     ".join(cands))
    else:
        import re as _re
        text = open(path, errors="replace").read()
        m = _re.search(r"\*\*Unresolved:\s*(\d+)", text)
        if not m:
            why = "%s does not state an unresolved count" % path
        elif int(m.group(1)) != 0:
            why = "%s reports %s unresolved finding(s)" % (path, m.group(1))
        else:
            #  Same staleness trap as the length audit: a clean report that
            #  predates the collections describes a different set of bins.
            newest = 0
            for d in (os.path.join(workdir, "collections", module),
                      os.path.join(workdir, module, "collections")):
                for f in glob.glob(os.path.join(d, "*.fasta")):
                    newest = max(newest, os.path.getmtime(f))
            if newest and os.path.getmtime(path) < newest:
                why = "%s is OLDER than the collections it describes" % path
    if why is None:
        print("bin audit: clean (%s)" % path)
        return
    if force:
        print("WARNING: building over unaudited bins -- %s" % why)
        return
    sys.exit(
        "REFUSED: %s\n"
        "  A bin that is a pile of partials, or a second copy of a protein\n"
        "  that already has a name, becomes a profile that calls somebody\n"
        "  else's locus BETTER than its own. No cutoff separates that, and\n"
        "  the only symptom is duplicate calls -- which look like a cutoff\n"
        "  problem and get chased through coverage runs.\n"
        "  Run audit_bins.py, resolve or --accept every finding, then build.\n"
        "  --force-uncurated overrides, for an experiment only."
        % why)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--module", required=True)
    ap.add_argument("--features", required=True,
                    help='JSON file: {"KEY": {"anno": "...", "params": {...}, "flags": [...]}, ...}')
    ap.add_argument("--only", default="", help="comma-separated subset")
    ap.add_argument("--force-uncurated", action="store_true",
                    help="build profiles although the length or bin audit is "
                         "missing or unresolved. For a deliberate experiment "
                         "only.")
    a = ap.parse_args()
    require_length_audit(a.workdir, a.module, a.force_uncurated)
    require_bin_audit(a.workdir, a.module, a.force_uncurated)
    spec = json.load(open(a.features))
    only = set(x for x in a.only.split(",") if x)
    print("building %s" % a.module)
    res = []
    for key, d in spec.items():
        if only and key not in only: continue
        #  A `special` feature is called by its own pass -- the splice or
        #  transcript-edit path -- from references in Splice-Variants/, not
        #  by a profile. install_module.py treats a special that HAS profiles
        #  as fatal, because the profile and the special pass would both fire
        #  on the locus. So the declaration has to survive as far as the spec
        #  and then be skipped HERE, rather than being filtered out upstream:
        #  setup_build.py used to drop specials from features.json entirely to
        #  keep them away from this loop, and because gen_module_json.py
        #  builds the module JSON from features.json, one regeneration then
        #  deleted the feature from the JSON. Adenoviridae lost K33 (33K,
        #  spliced -- PMID 25636034) from all four modules that way, with its
        #  splice references still sitting on disk and nothing reporting a
        #  gap. The flag now travels in the spec and the skip lives here.
        if d.get("special"):
            print("  %-14s special:%s -- no profile, called by its own pass"
                  % (key, d["special"]))
            continue
        params = dict(DEFAULTS); params.update(_normalise(d.get("params")))
        flags = list(d.get("flags", []))
        #  A read-through feature carries the read-through codon as an X in
        #  EVERY member, by construction. fasta-cluster-pssm-2.pl drops any
        #  sequence matching /[BJXZ]/ unless -x is given, so the collection
        #  empties, mmseqs is handed a 0-byte FASTA and dies with "Cannot seek
        #  at the end of file" -- which names the temp file, not the feature,
        #  and reads like an mmseqs fault rather than a filtered input.
        #  Coltivirus VP9RDT lost all 17 of its sequences this way. Any
        #  feature declaring internal_stop needs -x, so set it here rather
        #  than relying on every module to remember.
        if d.get("internal_stop") and "-x" not in flags:
            flags.append("-x")
            print("  %-14s internal_stop: adding -x so the read-through codon "
                  "does not filter the collection to nothing" % key)
        r = run_one(a.workdir, a.module, key, d["anno"], params, flags)
        if r: res.append(r)
    print("\n%d features built, %d PSSMs total" % (len(res), sum(r["pssms"] for r in res)))

main()
