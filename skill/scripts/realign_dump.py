#!/usr/bin/env python3
"""Re-align the three uniq.* files of a BV-BRC dump, and record what was lost.

The three files are LINE-INDEX ALIGNED -- that alignment is the join key -- and
BV-BRC has delivered them misaligned on every taxon this project has dumped:

    Hepaciviridae   100,855 / 103,166 / 100,290   565 md5 with no sequence
    Pestiviridae     11,318 /  11,318 /  11,311     7
    Pegivirus         2,767 /   2,767 /   2,764     3

The cause is always the same: the feature_sequence query returns nothing for a
few md5, so uniq.seq is short while uniq.md5 and uniq.id_ann are not. Left
alone, every downstream index is off by one from the first gap onward and the
wrong annotation is attached to the wrong sequence -- silently.

This keeps only rows whose md5 has a sequence, rewrites all three files in the
same order, and writes the dropped md5 to missing.md5 so the loss is on record
rather than inferred later from a line-count difference.

    python3 realign_dump.py --workdir . --prefix Pesti
    python3 realign_dump.py --workdir . --prefix Pesti --write
"""
import argparse, os, shutil, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--prefix", required=True, help="e.g. Pesti, Hepaci, Pegi")
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()
    H, T = os.path.abspath(a.workdir), a.prefix

    def p(x): return os.path.join(H, "%s.%s" % (T, x))
    for f in ("uniq.md5", "uniq.id_ann", "uniq.seq"):
        if not os.path.exists(p(f)):
            sys.exit("missing %s" % p(f))

    md5s = [l.rstrip("\n").split("\t")[0] for l in open(p("uniq.md5"))]
    ann = [l.rstrip("\n").split("\t") for l in open(p("uniq.id_ann"))]
    seq = {}
    for l in open(p("uniq.seq")):
        c = l.rstrip("\n").split("\t")
        if len(c) > 1 and c[0] and c[1]:
            seq[c[0]] = c[1]

    if len(ann) != len(md5s):
        print("WARNING: uniq.id_ann has %d rows against uniq.md5's %d. It is "
              "POSITIONAL against uniq.md5, so a length difference means the "
              "two were already out of step before any sequence went missing."
              % (len(ann), len(md5s)))

    keep = [i for i, m in enumerate(md5s) if m in seq]
    missing = [m for m in md5s if m not in seq]
    print("%s: %d rows, %d have a sequence, %d do not"
          % (T, len(md5s), len(keep), len(missing)))
    if not missing and len(ann) == len(md5s):
        print("already aligned; nothing to do"); return 0
    if not a.write:
        print("rerun with --write to rewrite the three files"); return 0

    bak = os.path.join(H, "pre_realign")
    os.makedirs(bak, exist_ok=True)
    for f in ("uniq.md5", "uniq.id_ann", "uniq.seq"):
        shutil.copy2(p(f), os.path.join(bak, "%s.%s" % (T, f)))

    with open(p("uniq.md5"), "w") as f:
        for i in keep: f.write(md5s[i] + "\n")
    with open(p("uniq.id_ann"), "w") as f:
        for i in keep:
            r = ann[i] if i < len(ann) else [""]
            f.write("%s\t%s\n" % (r[0] if r else "", r[1] if len(r) > 1 else ""))
    with open(p("uniq.seq"), "w") as f:
        for i in keep: f.write("%s\t%s\n" % (md5s[i], seq[md5s[i]]))
    open(os.path.join(H, "missing.md5"), "w").write(
        "\n".join(missing) + ("\n" if missing else ""))
    print("rewrote three files at %d rows; originals in pre_realign/, "
          "dropped md5 in missing.md5" % len(keep))
    return 0


sys.exit(main())
