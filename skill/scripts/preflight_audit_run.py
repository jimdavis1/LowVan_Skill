#!/usr/bin/env python3
"""Pre-flight for a coverage/eval run: prove every input is final first.

A long annotation run reads the module JSON, the PSSMs, the rep contigs and
the annotator scripts fresh for EVERY genome. Editing any of them mid-run
does not restart anything -- it splits the run, and the output is part old
rules and part new, which is worse than either. The 26 Oct lentivirus run was
launched at 09:07:56 and EIAV S2 gained copy_num at 09:14:13; the first
twelve exemplars were scored against a different rule set than the rest.

So: run this, fix anything it reports, and launch only on a clean pass.
It prints the newest mtime of everything the run will read, which is the
timestamp a launch must come AFTER.

    python3 preflight_audit_run.py --kit LowVan_Skill --runtime Viral_Annotation
"""
import argparse, hashlib, json, os, sys, time


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def newest(paths):
    best = (0, None)
    for p in paths:
        if os.path.isfile(p):
            m = os.path.getmtime(p)
            if m > best[0]:
                best = (m, p)
        elif os.path.isdir(p):
            for dp, _dn, fns in os.walk(p):
                for fn in fns:
                    f = os.path.join(dp, fn)
                    try:
                        m = os.path.getmtime(f)
                    except OSError:
                        continue
                    if m > best[0]:
                        best = (m, f)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kit", required=True)
    ap.add_argument("--runtime", required=True)
    ap.add_argument("--modules", default="",
                    help="comma-separated modules whose copy_num to print")
    ap.add_argument("--scope", default="",
                    help="modules this run actually annotates. Parity is "
                         "ENFORCED for these and advisory for the rest -- a "
                         "run over lentivirus holdings is not made wrong by "
                         "Lagovirus drift, and blocking on it would train you "
                         "to ignore the check.")
    a = ap.parse_args()
    bad = []

    print("1. annotator scripts: kit vs runtime")
    kitd = os.path.join(a.kit, "annotate")
    for fn in sorted(os.listdir(kitd)):
        if not fn.endswith((".pl", ".pm")):
            continue
        k, r = os.path.join(kitd, fn), os.path.join(a.runtime, fn)
        if not os.path.exists(r):
            print("     MISSING in runtime  %s" % fn); bad.append(fn); continue
        if os.path.getsize(k) == 0 or os.path.getsize(r) == 0:
            print("     ZERO BYTES (Box)    %s" % fn); bad.append(fn); continue
        if md5(k) != md5(r):
            print("     DIFFERS             %s" % fn); bad.append(fn)
    print("     %d script(s) out of parity" % len(bad))

    scope = {x.strip() for x in a.scope.split(",") if x.strip()}
    advisory = []
    print("\n2. module JSON: kit vs runtime, every module the kit ships")
    if scope:
        print("     enforced for: %s" % ", ".join(sorted(scope)))
    rj = json.load(open(os.path.join(a.runtime, "Viral_PSSM.json")))
    modsd = os.path.join(a.kit, "modules")
    njson = 0
    for m in sorted(os.listdir(modsd)):
        kp = os.path.join(modsd, m, "%s_Viral_PSSM.json" % m)
        if not os.path.exists(kp):
            continue
        njson += 1
        try:
            kb = json.load(open(kp))[m]
        except Exception as e:
            print("     UNREADABLE %s: %s" % (m, e)); bad.append(m); continue
        if m not in rj:
            print("     NOT IN RUNTIME      %s" % m); bad.append(m); continue
        if json.dumps(kb, sort_keys=True) != json.dumps(rj[m], sort_keys=True):
            if scope and m not in scope:
                advisory.append(m)
            else:
                print("     DIFFERS (IN SCOPE)  %s" % m); bad.append(m)
    print("     compared %d module(s)" % njson)
    if advisory:
        print("     %d out of parity but OUT OF SCOPE for this run, not blocking:"
              % len(advisory))
        print("       %s" % ", ".join(advisory))

    print("\n3. rule conformance in the runtime JSON")
    spl = [(m, k) for m, b in rj.items() if isinstance(b, dict)
           if (not scope or m in scope)
           for k, v in (b.get("features") or {}).items()
           if isinstance(v, dict) and v.get("special") == "splice" and v.get("copy_num")]
    if spl:
        print("     SPLICED FEATURES MARKED ESSENTIAL: %s" % spl); bad.extend(spl)
    else:
        print("     no spliced feature carries copy_num")

    if a.modules:
        print("\n4. copy_num as it will be read")
        for m in a.modules.split(","):
            m = m.strip()
            if m in rj:
                print("     %-6s %s" % (m, {k: ("Y" if v.get("copy_num") else "-")
                                            for k, v in sorted(rj[m]["features"].items())}))

    print("\n5. newest mtime among everything the run reads")
    for label, paths in (
            ("runtime Viral_PSSM.json", [os.path.join(a.runtime, "Viral_PSSM.json")]),
            ("kit annotate/", [kitd]),
            ("runtime *.pl", [os.path.join(a.runtime, f)
                              for f in os.listdir(a.runtime) if f.endswith(".pl")]),
            ("runtime Viral-PSSMs/", [os.path.join(a.runtime, "Viral-PSSMs")]),
            ("runtime Viral-Rep-Contigs/", [os.path.join(a.runtime, "Viral-Rep-Contigs")]),
            ("runtime Splice-Variants/", [os.path.join(a.runtime, "Splice-Variants")]),
            ("runtime Transcript-Editing/", [os.path.join(a.runtime, "Transcript-Editing")]),
            ("kit skill/scripts/", [os.path.join(a.kit, "skill", "scripts")]),
    ):
        m, p = newest(paths)
        print("     %-28s %s  %s" % (label,
                                     time.strftime("%H:%M:%S", time.localtime(m)) if m else "--------",
                                     os.path.basename(p) if p else "(none)"))
    allm, allp = newest([os.path.join(a.runtime, "Viral_PSSM.json"), kitd,
                         os.path.join(a.runtime, "Viral-PSSMs"),
                         os.path.join(a.runtime, "Viral-Rep-Contigs"),
                         os.path.join(a.runtime, "Splice-Variants"),
                         os.path.join(a.runtime, "Transcript-Editing"),
                         os.path.join(a.kit, "skill", "scripts")])
    print("\n     LAUNCH MUST BE AFTER %s  (%s)"
          % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(allm)), allp))

    if bad:
        print("\n  NOT READY: %d problem(s) above. Fix, then re-run this." % len(bad))
        return 1
    print("\n  READY: inputs are consistent. Launch now and change nothing until it ends.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
