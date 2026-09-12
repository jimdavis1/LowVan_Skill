#!/usr/bin/env python3
"""Audit features declared in the JSON that have no collection behind them.

A feature declared with nothing binding to it is invisible at runtime: the
annotator loads no PSSM for it and never emits it. Some of those are honest
placeholders for a protein the dump does not carry. Others are keys invented
from a literature survey that no binding rule ever connected, sometimes naming
a "family" whose members are not homologous to one another at all.

For each unbuilt feature this finds the genus its key points at, gathers that
genus's proteins that reached no specific feature, and asks the only question
that matters: **are the candidates a homology group?** Three verdicts:

  BINDABLE     candidates exist and are homologous -- write a binding rule
  NOT-A-GROUP  candidates exist but are mutually unrelated -- a positional
               label, not a protein family. Drop the key; the proteins belong
               in UNCHAR.
  NO-DATA      the genus has no unaccounted proteins at all -- the key is a
               placeholder for something the dump does not contain

    python3 audit_declared.py --workdir .
    python3 audit_declared.py --workdir . --key ARURHA_GX     # one feature

Needs blastp on PATH.
"""

import argparse
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict


def load_dump(workdir):
    pre = None
    for p in glob.glob(os.path.join(workdir, "*.uniq.md5")):
        pre = p[:-len(".uniq.md5")]
    if not pre:
        sys.exit("no <T>.uniq.md5 in %s" % workdir)
    md5 = [l.strip() for l in open(pre + ".uniq.md5") if l.strip()]
    ann = [l.rstrip("\n").split("\t") for l in open(pre + ".uniq.id_ann")]
    seq = [l.rstrip("\n").split("\t") for l in open(pre + ".uniq.seq")]
    m2 = {}
    for i, m in enumerate(md5):
        m2[m] = {"id": ann[i][0],
                 "ann": ann[i][1] if len(ann[i]) > 1 else "",
                 "seq": seq[i][1] if len(seq[i]) > 1 else ""}
    gen = {}
    for l in open(pre + ".id_name_gs"):
        f = l.rstrip("\n").split("\t")
        if len(f) >= 4:
            gen[f[0]] = f[3].strip()
    id2m = {}
    for l in open(pre + ".id_md5"):
        f = l.rstrip("\n").split("\t")
        if len(f) > 1:
            id2m[f[0]] = f[1]
    return m2, gen, id2m


def collection_ids(workdir):
    """Every feature id that already reached a specific (non-UNCHAR) collection."""
    placed = set()
    for fa in glob.glob(os.path.join(workdir, "collections", "*", "*.fasta")):
        feat = os.path.basename(fa)[:-6].split(".")[0]
        if feat == "UNCHAR":
            continue
        for line in open(fa, errors="replace"):
            if line.startswith(">"):
                placed.add(line[1:].split()[0])
    return placed


GREEK = {"a": "alpha", "b": "beta", "g": "gamma", "d": "delta"}


def guess_genus(key, symbol, genera):
    """Map a key stem onto the genus it abbreviates.

    Two shapes in use. A plain prefix (ARURHA -> Arurhavirus), and a
    contraction that swaps a Greek qualifier for its initial (ACYTO ->
    Alphacytorhabdovirus, BRICIN -> Betaricinrhavirus, GNUCLEO ->
    Gammanucleorhabdovirus). Handle both, and refuse to guess when the stem is
    ambiguous rather than attributing a feature to the wrong genus.
    """
    stem = re.sub(r"[^a-z]", "", (symbol or key).split("_")[0].lower())
    if len(stem) < 3:
        return None
    low = {g.lower(): g for g in genera}
    hits = [g for l, g in low.items() if l.startswith(stem)]
    if len(hits) == 1:
        return hits[0]
    # contraction: first letter is a Greek initial, the rest follows it
    word = GREEK.get(stem[0])
    if word:
        rest = stem[1:]
        hits = [g for l, g in low.items()
                if l.startswith(word) and l[len(word):].startswith(rest)]
        if len(hits) == 1:
            return hits[0]
    hits = [g for l, g in low.items() if l.startswith(stem[:5])]
    return hits[0] if len(hits) == 1 else None


def groups(seqs, tmp, pident=25.0):
    """Single-linkage clusters of a candidate pool. Returns a list of index sets.

    The question an audit has to answer is not "does anything here match
    anything" but "how many distinct proteins are in this pool", because that
    is what bounds how many features the genus can support.
    """
    n = len(seqs)
    if n == 0:
        return []
    if n == 1:
        return [{0}]
    fa = os.path.join(tmp, "h.faa")
    with open(fa, "w") as fh:
        for i, s in enumerate(seqs):
            fh.write(">s%d\n%s\n" % (i, s))
    subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "prot",
                    "-out", os.path.join(tmp, "hdb")], capture_output=True)
    r = subprocess.run(
        ["blastp", "-query", fa, "-db", os.path.join(tmp, "hdb"), "-outfmt",
         "6 qseqid sseqid pident length", "-evalue", "1e-3", "-word_size", "2",
         "-max_target_seqs", str(n + 5)],
        capture_output=True, text=True)
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for line in r.stdout.splitlines():
        q, s, pid, alen = line.split("\t")
        if q == s:
            continue
        i, j = int(q[1:]), int(s[1:])
        if float(pid) >= pident and int(alen) >= 0.5 * min(len(seqs[i]), len(seqs[j])):
            parent[find(i)] = find(j)
    out = defaultdict(set)
    for i in range(n):
        out[find(i)].add(i)
    return sorted(out.values(), key=len, reverse=True)


def homology(seqs, tmp, pident=25.0):
    """Is this set a homology group? Returns (n_pairs_linked, max_pident)."""
    if len(seqs) < 2:
        return 0, 0.0
    fa = os.path.join(tmp, "h.faa")
    with open(fa, "w") as fh:
        for i, s in enumerate(seqs):
            fh.write(">s%d\n%s\n" % (i, s))
    subprocess.run(["makeblastdb", "-in", fa, "-dbtype", "prot",
                    "-out", os.path.join(tmp, "hdb")], capture_output=True)
    r = subprocess.run(
        ["blastp", "-query", fa, "-db", os.path.join(tmp, "hdb"), "-outfmt",
         "6 qseqid sseqid pident", "-evalue", "10", "-word_size", "2",
         "-max_target_seqs", str(len(seqs) + 5)],
        capture_output=True, text=True)
    best, linked = 0.0, set()
    for line in r.stdout.splitlines():
        q, s, p = line.split("\t")
        if q == s:
            continue
        p = float(p)
        best = max(best, p)
        if p >= pident:
            linked.add(tuple(sorted((q, s))))
    return len(linked), best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--json", help="module JSON")
    ap.add_argument("--key", help="audit one feature key")
    ap.add_argument("--pident", type=float, default=25.0,
                    help="identity above which two candidates count as homologous")
    args = ap.parse_args()

    workdir = os.path.abspath(args.workdir)
    jf = args.json
    if not jf:
        c = [p for p in glob.glob(os.path.join(workdir, "*_Viral_PSSM.json"))
             if os.path.basename(p) != "Viral_PSSM.json"]
        jf = c[0]
    mod = json.load(open(jf))
    m2, gen, id2m = load_dump(workdir)
    placed = collection_ids(workdir)
    genera = sorted(set(gen.values()) - {""})

    # unaccounted proteins per genus: everything not already in a named collection
    loose = defaultdict(list)
    for fid, m in id2m.items():
        rec = m2.get(m)
        if not rec or not rec["seq"]:
            continue
        gid = re.match(r"fig\|(\d+\.\d+)\.", fid)
        g = gen.get(gid.group(1)) if gid else None
        if not g or rec["id"] in placed:
            continue
        loose[g].append(rec)

    todo = []
    for module, block in mod.items():
        for key, e in block.get("features", {}).items():
            if args.key and key != args.key:
                continue
            if e.get("special"):
                continue
            # derived peptides belong to a parent CDS, not to a genus; they are
            # build_sp_pssms.py's business, not this audit's
            if re.search(r"_(MAT|SP|MATURE)$", key):
                continue
            has = glob.glob(os.path.join(workdir, "collections", module, key + ".fasta"))
            pss = glob.glob(os.path.join(workdir, "Alignments", module, key, "pssms", "*.pssm"))
            if has or pss:
                continue
            todo.append((module, key, e))

    bygenus = defaultdict(list)
    for module, key, e in todo:
        bygenus[guess_genus(key, e.get("gene_symbol", ""), genera)].append((module, key))
    print("%d declared feature(s) with no collection and no PSSM, "
          "across %d genus target(s)\n" % (len(todo), len(bygenus)))
    tmp = tempfile.mkdtemp(prefix="audit.")
    try:
        for g in sorted(bygenus, key=lambda x: (x is None, x or "")):
            keys = sorted(bygenus[g])
            if g is None:
                print("  UNMAPPED KEYS (%d) -- stem matches no genus in this dump:" % len(keys))
                print("     %s" % ", ".join("%s/%s" % (m, k) for m, k in keys))
                print()
                continue
            cands = loose.get(g, [])
            uniq, seen = [], set()
            for r in cands:
                if r["seq"] not in seen:
                    seen.add(r["seq"]); uniq.append(r)
            gs = groups([r["seq"] for r in uniq], tmp, args.pident)
            print("  %s  --  %d key(s) declared: %s"
                  % (g, len(keys), ", ".join(k for _m, k in keys)))
            if not uniq:
                print("     no unaccounted proteins in this genus  ->  NO-DATA, drop or keep as a placeholder\n")
                continue
            print("     %d unique unaccounted protein(s) forming %d homology group(s):"
                  % (len(uniq), len(gs)))
            for gi, idx in enumerate(gs, 1):
                lens = sorted(len(uniq[i]["seq"]) for i in idx)
                names = sorted({uniq[i]["ann"] or "(blank)" for i in idx})
                print("       group %d: %2d seq, %4d-%-4d aa, called %s"
                      % (gi, len(idx), lens[0], lens[-1],
                         "; ".join(n[:30] for n in names[:3])))
            if len(gs) < len(keys):
                print("     -> only %d real group(s) for %d keys: DROP %d, the labels are positional"
                      % (len(gs), len(keys), len(keys) - len(gs)))
            elif all(len(i) < 2 for i in gs):
                print("     -> every group is a singleton: NOT BUILDABLE, leave in UNCHAR")
            else:
                print("     -> %d buildable group(s): write binding rules" % sum(1 for i in gs if len(i) >= 2))
            print()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print("""  A genus with fewer homology groups than declared keys has been given
  positional labels for one protein family. Keep one key per group; the rest
  are names for a slot, not for a protein.""")
    return 0


if __name__ == "__main__":
    sys.exit(main())
