#!/usr/bin/env python3
"""Check the kit's annotator scripts against the runtime copies.

There are two copies of every annotator: the repository's, in
LowVan_Skill/annotate/, and the runtime target's, in Viral_Annotation/.
env.sh puts the KIT on PATH, and New-annotate-viral-taxon.pl invokes the
helpers by bare name, so the kit copy is the one that actually runs inside
the GTO pipeline -- while hand-driven runs that spell out a path usually
reach for the runtime one.

When the two drift, the same GTO gets two different answers and neither
copy reports an error. On 26 Oct 2026 the kit's get_splice_variant_features.pl
was five weeks old and had no splice_id/splice_cov support: it found no
splice sites, said "Splice variant features should exist for: EIAV, but none
were found", and exited 0. Every coverage audit run through that path
measured its module with the spliced features missing.

    python3 audit_script_parity.py --kit LowVan_Skill --runtime Viral_Annotation

Exits non-zero if any pair differs or is missing, so it can gate a build.
"""
import argparse, hashlib, os, sys


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kit", required=True, help="LowVan_Skill checkout")
    ap.add_argument("--runtime", required=True, help="Viral_Annotation runtime dir")
    a = ap.parse_args()

    kitd = os.path.join(a.kit, "annotate")
    if not os.path.isdir(kitd):
        sys.exit("no such directory: %s" % kitd)

    bad = 0
    for fn in sorted(os.listdir(kitd)):
        if not fn.endswith((".pl", ".pm")):
            continue
        k = os.path.join(kitd, fn)
        r = os.path.join(a.runtime, fn)
        if not os.path.exists(r):
            print("  MISSING in runtime   %s" % fn); bad += 1; continue
        ks, rs = os.path.getsize(k), os.path.getsize(r)
        #  Box serves short and empty reads while exiting 0, so a zero-byte
        #  side is a failed read, not an empty script.
        if ks == 0 or rs == 0:
            print("  UNREADABLE (0 bytes) %s  kit=%d runtime=%d" % (fn, ks, rs))
            bad += 1; continue
        if md5(k) != md5(r):
            newer = "runtime" if os.path.getmtime(r) > os.path.getmtime(k) else "kit"
            print("  DIFFERS  %-38s kit=%d runtime=%d  (%s is newer)"
                  % (fn, ks, rs, newer))
            bad += 1
        else:
            print("  ok       %s" % fn)

    if bad:
        print("\n  %d script(s) out of parity. The kit copy is what PATH resolves"
              "\n  inside New-annotate-viral-taxon.pl, so a stale kit silently"
              "\n  changes pipeline results. Establish which way the fix flowed"
              "\n  BEFORE copying either direction." % bad)
        return 1
    print("\n  all annotator scripts identical between kit and runtime")
    return 0


if __name__ == "__main__":
    sys.exit(main())
