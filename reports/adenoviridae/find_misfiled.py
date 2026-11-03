#!/usr/bin/env python3
"""Find genomes whose source genus puts them in the wrong module.

build_cols.py assigns a protein to a module by the BV-BRC genus of its genome.
A genome filed under the wrong genus therefore does two kinds of damage at
once, and nothing downstream can see either:

  1. its proteins join the wrong module's COLLECTIONS, so that module's PSSMs
     are trained on a foreign protein;
  2. it becomes eligible as one of that module's 25 REP CONTIGS, and a rep
     contig is a BLASTn subject -- so it wins routing for the whole clade it
     really belongs to.

Measured on Parvoviridae, 29 Oct 2026: `Protoparvovirus sp. HNXY21` is a waterfowl
dependoparvovirus (98.5% identical over 5,073 nt to goose parvovirus). It was
picked as a Protoparvovirus rep contig and took 37 of 300 sampled
dependoparvoviruses with it.

THE TEST. All-vs-all blastp of every module's anchor collection, then for each
genome count how many >=90%-identity partners live in its own module against
how many live elsewhere.

Do NOT use best-single-hit with a margin. A CLUSTER of misfiled genomes
vouches for itself inside the wrong module: the three "Protoparvovirus sp."
goose parvoviruses are each other's best own-module hit, so the margin test
passes them. It found 8 of the 18, and not HNXY21 -- the worst one.

HEXON is the anchor for Adenoviridae: every adenovirus has exactly one, it is
the major capsid protein and the most conserved gene in the family, and it is
the gene the whole typing literature is built on. POL would serve equally. A
genus-specific accessory would not, because its ABSENCE is not evidence of
anything -- and in this family the accessory sets are exactly what differ
between the four modules.
"""
import collections, json, os, subprocess, sys, tempfile

W = os.path.dirname(os.path.abspath(__file__)) or "."
ANCHOR = "HEXON"
MIN_ID = 90.0        # below this, two adenoviruses are not the same virus
MIN_LEN = 200        # a short local hit is not evidence of co-membership
FOREIGN_RATIO = 3    # foreign partners must outnumber own by this much


def modules():
    return sorted(d for d in os.listdir(W)
                  if os.path.isdir(os.path.join(W, d, "collections")))


def collect(path, root=None):
    """module|genome|seqname -> sequence, from every module's anchor.

    `root` reads a FLAT snapshot (<root>/<module>/NS1.fasta) instead of the
    live tree, so the detector can be re-run against the collections as they
    were BEFORE a fix was applied.
    """
    out = []
    for m in (sorted(os.listdir(root)) if root else modules()):
        p = (os.path.join(root, m, "%s.fasta" % ANCHOR) if root else
             os.path.join(W, m, "collections", "%s.fasta" % ANCHOR))
        if not os.path.exists(p):
            print("  %-22s no %s collection, skipped" % (m, ANCHOR))
            continue
        name = None
        for l in open(p):
            if l.startswith(">"):
                name = l[1:].strip()
            elif name:
                out.append((m, name, l.strip()))
                name = None
    with open(path, "w") as fh:
        for m, name, s in out:
            fh.write(">%s|%s\n%s\n" % (m, name, s))
    return len(out)


def UNNAMED(g):
    """True if the genome's organism name makes no taxonomic claim."""
    n = _NAMES.get(g, "").lower()
    return (not n) or n.startswith(("unidentified", "unclassified",
                                    "untyped")) or " sp." in n


_NAMES = {}


def main():
    names, genus = {}, {}
    gid = os.path.join(W, "A.id_name_gs.inscope")
    if os.path.exists(gid):
        for l in open(gid, errors="replace"):
            f = l.rstrip("\n").split("\t")
            if len(f) >= 2:
                names[f[0]] = f[1]
                genus[f[0]] = f[3] if len(f) >= 4 and f[3] else ""
    _NAMES.update(names)

    root = None
    if "--from" in sys.argv:
        root = sys.argv[sys.argv.index("--from") + 1]
        print("  reading the collection snapshot at %s" % root)
    t = tempfile.mkdtemp(prefix="adeno.misfiled.")
    try:
        faa = os.path.join(t, "anchor.faa")
        n = collect(faa, root)
        print("  %d %s sequences across %d modules" % (n, ANCHOR, len(modules())))
        subprocess.run(["makeblastdb", "-in", faa, "-dbtype", "prot",
                        "-out", os.path.join(t, "db")], check=True,
                       capture_output=True)
        #  stderr is NOT discarded. A `blastp: command not found` with stderr
        #  sent to /dev/null returns zero hits, which reads exactly like a
        #  clean negative result -- that happened once already in this build.
        r = subprocess.run(["blastp", "-query", faa, "-db", os.path.join(t, "db"),
                            "-outfmt", "6 qseqid sseqid pident length",
                            "-evalue", "1e-20", "-num_threads", "4"],
                           capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit("blastp failed: %s" % r.stderr.strip()[:400])
        lines = r.stdout.splitlines()
    finally:
        import shutil
        shutil.rmtree(t, ignore_errors=True)
    print("  %d alignments" % len(lines))

    part = collections.defaultdict(collections.Counter)
    pgen = collections.defaultdict(collections.Counter)
    top = {}
    for l in lines:
        f = l.split("\t")
        if len(f) < 4:
            continue
        if float(f[2]) < MIN_ID or int(f[3]) < MIN_LEN:
            continue
        q = f[0].split("|"); s = f[1].split("|")
        if q[1] == s[1]:            # same genome is not independent evidence
            continue
        #  AN UNNAMED GENOME IS NOT A VOUCHER. "unidentified adenovirus" is a
        #  single BV-BRC taxon (10535) covering 111 genomes of every genus,
        #  and the whole bucket is filed under Mastadenovirus by default. It
        #  carries no taxonomic claim at all, so it must never be used as
        #  EVIDENCE about where a named genome belongs -- only as a query.
        #
        #  MEASURED: this ran once without the guard and reported that
        #  `Penguin siadenovirus A` belonged in Mastadenovirus, on the
        #  strength of three >=98% partners. Those three partners were
        #  10535.510/511/512, all "unidentified adenovirus", and they are
        #  themselves siadenoviruses in the wrong module. The one named
        #  genome in the clade was accused by the three unnamed ones.
        if UNNAMED(s[1]):
            continue
        part[(q[0], q[1])][s[0]] += 1
        pgen[(q[0], q[1], s[0])][genus.get(s[1], "")] += 1
        k = (q[0], q[1], s[0])
        top[k] = max(top.get(k, 0.0), float(f[2]))

    mis = {}
    for (mod, g), c in part.items():
        own = c.get(mod, 0)
        foreign = sum(v for k, v in c.items() if k != mod)
        if not foreign:
            continue
        #  A SINGLE >=95% FULL-LENGTH HEXON MATCH IS CONCLUSIVE. Requiring
        #  two foreign vouchers assumes the right clade has more than one
        #  named member in the dump, and sometimes it has exactly one:
        #  10535.510/511/512 are 97.9-98.9% identical to Penguin siadenovirus
        #  A, the ONLY named siadenovirus of that clade here, so own=0 and
        #  foreign=1 and the two-voucher rule missed all three.
        best = max((v for (m, g2, t), v in top.items()
                    if m == mod and g2 == g and t != mod), default=0.0)
        if (foreign >= FOREIGN_RATIO * max(own, 1)
                or (own == 0 and foreign >= 2)
                or (own == 0 and best >= 95.0)):
            win = max(((v, k) for k, v in c.items() if k != mod))[1]
            #  Carry the GENUS its evidence points at, not just the module.
            #  build_rep_contigs.py selects references per genus, so a
            #  reassigned genome landing in a module with a made-up genus
            #  label would get its own private bucket and a reference slot
            #  it has not earned.
            gc = pgen.get((mod, g, win))
            mis[g] = {"from": mod, "to": win, "own": own, "foreign": foreign,
                      "genus": gc.most_common(1)[0][0] if gc else "",
                      "name": names.get(g, "")}

    by = collections.defaultdict(list)
    for g, d in mis.items():
        by[d["from"]].append((g, d))
    for mod in sorted(by):
        print("\n%-22s %d misfiled" % (mod, len(by[mod])))
        for g, d in sorted(by[mod], key=lambda kv: -kv[1]["foreign"]):
            print("   %-14s %-44s own=%-4d foreign=%-4d -> %s"
                  % (g, d["name"][:44], d["own"], d["foreign"], d["to"]))
    #  ---- MERGE, never replace ------------------------------------
    #  Once a reassignment has been applied, the genome sits in the right
    #  module and its evidence is GONE -- a re-run against the corrected
    #  collections sees nothing wrong with it. Overwriting therefore
    #  silently undoes every fix already made: this ran once and dropped 18
    #  genomes down to the 2 the first pass had missed. The file is the
    #  accumulated record of what homology said, so it accumulates.
    out = os.path.join(W, "misfiled.json")
    prev = json.load(open(out)) if os.path.exists(out) else {}
    fresh = {g: d for g, d in mis.items() if g not in prev}
    merged = dict(prev); merged.update(fresh)
    json.dump(merged, open(out, "w"), indent=1, sort_keys=True)
    print("\n  %d already recorded + %d new = %d genomes -> %s"
          % (len(prev), len(fresh), len(merged), out))


main()
