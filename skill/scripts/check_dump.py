#!/usr/bin/env python3
"""Verify a BV-BRC dump before building anything from it.

Two properties the whole pipeline assumes and nothing else checks.

**Completeness.** Every distinct protein md5 named in `<T>.id_md5` must have a
sequence in `<T>.uniq.seq`. A missing one is invisible: the protein cannot be
binned, clustered, aligned or profiled, and no error is raised anywhere. This
is easy to get wrong, because the step that selects protein-coding features is
usually a positive filter on the feature id:

    grep "CDS\\|mat" ...        # silently drops every peg feature

On Rhabdoviridae that one omission cost 15,069 distinct sequences -- 45% of the
protein diversity -- and produced two modules that looked declared-but-empty
for reasons that had nothing to do with the data.

**Index alignment.** `uniq.md5`, `uniq.id_ann` and `uniq.seq` are joined by line
number, not by key. Row *n* of each must describe the same sequence. They are
built by three separate SOLR queries, so this holds only as long as the query
returns rows in input order -- a property worth verifying rather than assuming,
since pagination or parallelism would break it silently.

    python3 check_dump.py --workdir .
    python3 check_dump.py --workdir . --prefix Rhabdo

Exit is non-zero if either property fails.
"""

import argparse
import glob
import os
import re
import sys
from collections import Counter


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--prefix", help="dump prefix (default: found from *.uniq.md5)")
    args = ap.parse_args()

    workdir = os.path.abspath(args.workdir)
    pre = os.path.join(workdir, args.prefix) if args.prefix else None
    if not pre:
        c = glob.glob(os.path.join(workdir, "*.uniq.md5"))
        if not c:
            sys.exit("no <prefix>.uniq.md5 in %s" % workdir)
        pre = c[0][:-len(".uniq.md5")]
    name = os.path.basename(pre)
    print("dump: %s\n" % name)

    need = {"id_md5": pre + ".id_md5", "id_name_gs": pre + ".id_name_gs",
            "uniq.md5": pre + ".uniq.md5", "uniq.id_ann": pre + ".uniq.id_ann",
            "uniq.seq": pre + ".uniq.seq"}
    for k, p in need.items():
        if not os.path.exists(p):
            sys.exit("missing %s" % p)
        print("  %-14s %8d lines" % (k, sum(1 for _ in open(p, errors="replace"))))

    md5 = [l.strip() for l in open(need["uniq.md5"])]
    ann = [l.rstrip("\n").split("\t") for l in open(need["uniq.id_ann"], errors="replace")]
    seq = [l.rstrip("\n").split("\t") for l in open(need["uniq.seq"], errors="replace")]
    id2m, blank = {}, 0
    for l in open(need["id_md5"], errors="replace"):
        f = l.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        if not f[1].strip():
            blank += 1
            continue
        id2m[f[0]] = f[1].strip()

    fails = 0
    print("\n--- index alignment ---")
    if not (len(md5) == len(ann) == len(seq)):
        print("  FAIL  the three uniq.* files have different line counts")
        fails += 1
    else:
        print("  ok    all three uniq.* files are %d lines" % len(md5))
        bad = [i for i in range(len(md5)) if (seq[i][0] if seq[i] else "") != md5[i]]
        print("  %-5s uniq.seq row md5 matches uniq.md5 row  (%d mismatch)"
              % ("ok" if not bad else "FAIL", len(bad)))
        fails += bool(bad)
        bad2 = [i for i in range(len(md5)) if id2m.get(ann[i][0] if ann[i] else "") != md5[i]]
        print("  %-5s uniq.id_ann row resolves to uniq.md5 row  (%d mismatch)"
              % ("ok" if not bad2 else "FAIL", len(bad2)))
        fails += bool(bad2)
        for i in bad2[:3]:
            print("          row %d: %s -> %s, expected %s"
                  % (i, ann[i][0], id2m.get(ann[i][0]), md5[i]))

    print("\n--- completeness ---")
    have = set(md5)
    want = set(id2m.values())
    missing = want - have
    print("  distinct md5 named by id_md5   : %d" % len(want))
    print("  distinct md5 with a sequence   : %d" % len(have))
    if missing:
        fails += 1
        print("  MISSING                        : %d  (%.0f%% of the protein diversity)"
              % (len(missing), 100.0 * len(missing) / max(len(want), 1)))
        # which id shapes are affected -- this names the filter that dropped them
        shape = Counter()
        for fid, m in id2m.items():
            if m in missing:
                t = re.search(r"\.([A-Za-z_]+)\.\d+$", fid)
                shape[t.group(1) if t else "?"] += 1
        tot = Counter()
        for fid in id2m:
            t = re.search(r"\.([A-Za-z_]+)\.\d+$", fid)
            tot[t.group(1) if t else "?"] += 1
        print("\n  affected feature-id shapes:")
        for t, n in shape.most_common(8):
            print("      %-14s %6d of %-6d missing  (%.0f%%)" % (t, n, tot[t], 100.0 * n / tot[t]))
        worst = shape.most_common(1)[0][0]
        print("\n  Almost certainly a positive filter on the feature id that does not")
        print('  name "%s". Check the step that builds uniq.md5.' % worst)
    else:
        print("  ok    every md5 named by id_md5 has a sequence")

    if blank:
        print("\n  note: %d id_md5 rows have a blank md5 (RNA, misc_feature, gap ...)."
              % blank)
        print("        Those are not protein-coding and are correctly absent.")

    print("\n%s" % ("PASS" if not fails else "%d CHECK(S) FAILED" % fails))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
