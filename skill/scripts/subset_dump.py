#!/usr/bin/env python3
"""Cut a family protein dump down to one genus, keeping the five files consistent.

A family dump is the efficient thing to pull from BV-BRC, but a module is built
per genus, and the engine processes whatever is in the workdir. Symlinking the
family dump into a genus workdir silently builds the module from the whole
family -- that happened once on Caliciviridae and was caught only because the
collections came back far too large.

Subsetting is not a grep, because of the three `uniq.*` files:

    <T>.uniq.md5      md5                                  one per unique sequence
    <T>.uniq.id_ann   representative_feature_id \t annotation
    <T>.uniq.seq      md5 \t sequence

`uniq.md5` and `uniq.id_ann` are **line-index aligned**, and the representative
feature id in `uniq.id_ann` names one arbitrary feature carrying that sequence.
Drop the genomes that feature belongs to and the representative is an orphan:
the md5 survives, the id it is paired with does not, and `check_dump.py` fails
on exactly those rows. So a representative outside the subset has to be
**re-pointed** at a feature inside it, not merely carried across.

    python3 subset_dump.py --workdir . --in-prefix Cali --out-prefix S \
            --genus Sapovirus --out-dir sapo

Writes the five files under --out-dir and reports what it kept. Run
check_dump.py on the result before building anything.
"""
import argparse, collections, os, sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".", help="where the family dump lives")
    ap.add_argument("--in-prefix", required=True, help="e.g. Cali")
    ap.add_argument("--out-prefix", required=True, help="e.g. S")
    ap.add_argument("--genus", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--genus-column", type=int, default=3,
                    help="0-based column of genus in <T>.id_name_gs (default 3)")
    a = ap.parse_args()

    W, I, O = a.workdir, a.in_prefix, a.out_prefix
    os.makedirs(a.out_dir, exist_ok=True)
    src = lambda ext: os.path.join(W, "%s.%s" % (I, ext))
    dst = lambda ext: os.path.join(a.out_dir, "%s.%s" % (O, ext))

    for ext in ("id_md5", "id_name_gs", "uniq.md5", "uniq.id_ann", "uniq.seq"):
        if not os.path.exists(src(ext)):
            sys.exit("missing %s" % src(ext))

    #  1. the genomes of this genus
    genomes, kept_meta = set(), []
    for line in open(src("id_name_gs")):
        p = line.rstrip("\n").split("\t")
        if len(p) > a.genus_column and p[a.genus_column] == a.genus:
            genomes.add(p[0]); kept_meta.append(line)
    if not genomes:
        sys.exit("no genomes with genus %r -- check --genus-column" % a.genus)

    #  2. the features on those genomes, and every md5 they use. A feature id
    #     is fig|<genome_id>.<type>.<n>, and genome_id itself contains a dot.
    def genome_of(fid):
        tail = fid.split("|")[-1]
        return ".".join(tail.split(".")[:2])

    md5_wanted = set()
    feat_by_md5 = collections.defaultdict(list)
    kept_feats = []
    for line in open(src("id_md5")):
        p = line.rstrip("\n").split("\t")
        if len(p) < 2:
            continue
        if genome_of(p[0]) in genomes:
            kept_feats.append(line)
            md5_wanted.add(p[1])
            feat_by_md5[p[1]].append(p[0])

    #  3. the uniq triple, line-aligned, with the representative re-pointed
    um = [l.strip() for l in open(src("uniq.md5"))]
    ua = [l.rstrip("\n") for l in open(src("uniq.id_ann"))]
    if len(um) != len(ua):
        sys.exit("uniq.md5 and uniq.id_ann are not line-aligned: %d vs %d"
                 % (len(um), len(ua)))

    out_md5, out_ann, repointed = [], [], 0
    for m, row in zip(um, ua):
        if m not in md5_wanted:
            continue
        parts = row.split("\t")
        rep = parts[0] if parts else ""
        anno = parts[1] if len(parts) > 1 else ""
        inside = feat_by_md5.get(m) or []
        if rep not in inside:
            if not inside:
                continue                       # nothing in-subset carries it
            rep = inside[0]; repointed += 1
        out_md5.append(m)
        out_ann.append("%s\t%s" % (rep, anno))

    keep_md5 = set(out_md5)
    seqs = 0
    with open(dst("uniq.seq"), "w") as fh:
        for line in open(src("uniq.seq")):
            p = line.rstrip("\n").split("\t")
            if p and p[0] in keep_md5:
                fh.write(line); seqs += 1

    open(dst("id_md5"), "w").writelines(kept_feats)
    open(dst("id_name_gs"), "w").writelines(kept_meta)
    open(dst("uniq.md5"), "w").write("".join(x + "\n" for x in out_md5))
    open(dst("uniq.id_ann"), "w").write("".join(x + "\n" for x in out_ann))

    print("  %s -> %s/%s.*" % (a.genus, a.out_dir, O))
    print("    genomes           %7d" % len(genomes))
    print("    protein features  %7d" % len(kept_feats))
    print("    unique sequences  %7d  (sequences written %d)" % (len(out_md5), seqs))
    print("    representatives re-pointed into the subset %d" % repointed)
    if seqs != len(out_md5):
        print("    WARNING: %d md5 kept with no sequence row" % (len(out_md5) - seqs))
    print("  now run: python3 check_dump.py --workdir %s" % a.out_dir)


main()
