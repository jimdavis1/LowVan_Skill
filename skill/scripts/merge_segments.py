#!/usr/bin/env python3
"""Reassemble multi-segment genomes split across one record per segment.

Genome collections store a segmented virus as one record per segment, each with
its own identifier and each a separate "genome" to any downstream count. Left
alone they inflate the denominator and make a fully annotated genome look like
several partial ones: 56 Dichorhavirus records are 28 genomes, and scoring them
as 56 halves the apparent completeness of a module that did nothing wrong.

Merging is gated on the module's own declared segment count, not on whatever
happens to share a name. A name with more parts than the module expects is a
name collision (different isolates, same strain string) and is left alone; a
name with fewer is a genuinely incomplete genome and is left alone too. Only an
exact match merges.

    python3 merge_segments.py --json module.json --fasta-dir contigs \
            --metadata genome_metadata.tsv --out merged/

Declare the count per module in the json as "segments": N. Modules without it
are treated as unsegmented and pass through untouched.
"""
import argparse, collections, json, os, re, sys


def norm(n):
    """Strip the segment designator so the parts of one genome share a key."""
    n = re.sub(r"\s+segment\s+\S+\s*$", "", n, flags=re.I)
    n = re.sub(r"\s+(RNA\s*)?[0-9]+\s*$", "", n)
    return re.sub(r"\s+", " ", n).strip().lower()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True)
    ap.add_argument("--fasta-dir", required=True)
    ap.add_argument("--metadata", required=True,
                    help="TSV with a genome id column and a genome name column")
    ap.add_argument("--id-col", default="genome.genome_id")
    ap.add_argument("--name-col", default="genome.genome_name")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    mod = json.load(open(args.json))
    counts = {m: b["segments"] for m, b in mod.items()
              if isinstance(b, dict) and b.get("segments", 1) > 1}
    if not counts:
        print("  no module declares segments > 1; nothing to merge")
        return 0
    for m, n in sorted(counts.items()):
        print("  %-24s expects %d segments" % (m, n))
    expected = set(counts.values())

    hdr, rows = None, []
    for line in open(args.metadata, errors="replace"):
        c = line.rstrip("\n").split("\t")
        if hdr is None:
            hdr = c; continue
        rows.append(dict(zip(hdr, c)))
    for col in (args.id_col, args.name_col):
        if col not in hdr:
            print("  ERROR: no column %r in %s" % (col, args.metadata), file=sys.stderr)
            return 1

    by = collections.defaultdict(list)
    for r in rows:
        gid = r.get(args.id_col, "").strip()
        nm = r.get(args.name_col, "").strip()
        if gid and nm:
            by[norm(nm)].append((gid, nm))

    os.makedirs(args.out, exist_ok=True)
    merged = kept = collided = partial = 0
    log = open(os.path.join(args.out, "merge_log.tsv"), "w")
    log.write("status\tn_records\tmerged_id\tgenome_name\tmembers\n")
    for key, members in sorted(by.items()):
        n = len(members)
        if n == 1:
            kept += 1; continue
        if n not in expected:
            st = "collision_not_merged" if n > max(expected) else "incomplete_not_merged"
            collided += n if n > max(expected) else 0
            partial += n if n <= max(expected) else 0
            log.write("%s\t%d\t\t%s\t%s\n"
                      % (st, n, members[0][1], ",".join(g for g, _ in members)))
            kept += n
            continue
        members.sort()
        out_id = members[0][0]
        path = os.path.join(args.out, out_id + ".fna")
        wrote = 0
        with open(path, "w") as fh:
            for gid, nm in members:
                src = os.path.join(args.fasta_dir, gid + ".fna")
                if not os.path.exists(src):
                    continue
                for line in open(src, errors="replace"):
                    fh.write(line)
                wrote += 1
        if wrote != n:
            os.remove(path); kept += n
            log.write("missing_fasta\t%d\t\t%s\t%s\n"
                      % (n, members[0][1], ",".join(g for g, _ in members)))
            continue
        merged += 1
        log.write("merged\t%d\t%s\t%s\t%s\n"
                  % (n, out_id, members[0][1], ",".join(g for g, _ in members)))
    log.close()

    print("\n  %d record(s) merged into %d genome(s)"
          % (sum(len(v) for v in by.values()) - kept, merged))
    print("  %d record(s) left as-is" % kept)
    if collided:
        print("  %d record(s) share a name but exceed the segment count "
              "(name collision, not merged)" % collided)
    if partial:
        print("  %d record(s) share a name but fall short of it "
              "(incomplete genome, not merged)" % partial)
    print("  log: %s" % os.path.join(args.out, "merge_log.tsv"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
