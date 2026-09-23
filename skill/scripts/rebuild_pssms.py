#!/usr/bin/env python3
"""Detect alignments edited since their PSSM was built, and rebuild those PSSMs.

Why this exists: Box sync flattens mtimes, so file timestamps cannot tell you
which alignments you have hand-curated since the last pipeline run. This instead
rebuilds every PSSM from its current alignment and compares bytes. A PSSM that
does not reproduce is stale, and the alignment behind it was edited.

Note that comparing only the PSSM's embedded query sequence is NOT sufficient --
deleting a non-master sequence leaves the master, and therefore the query,
untouched while still changing every score in the profile.

    python3 rebuild_pssms.py --workdir .            # report only
    python3 rebuild_pssms.py --workdir . --write    # rebuild the stale ones in place

Requires psiblast on PATH (conda env `lowvan`). Alignment files are never
modified -- only pssms/*.pssm are written.
"""
import argparse, glob, os, re, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from norm_pssm import norm as _norm   # noqa: E402
# order matters: the first directory holding <n>.fa is the PSSM's source
SUBDIRS = ("corrected_alis", "reclustered_alis", "truncated_alis", "alis")


def pssm_meta(path):
    """Return (local-id, title) from a PSSM's query block."""
    head = open(path, errors="replace").read(4000)
    lid = re.search(r'local str "([^"]*)"', head)
    title = re.search(r'title "([^"]*)"', head)
    return (lid.group(1) if lid else None), (title.group(1) if title else "")


def master_seq(fa):
    """First sequence of an alignment, still gapped."""
    hdr, cur = None, []
    for line in open(fa):
        if line.startswith(">"):
            if hdr is not None:
                return "".join(cur)
            hdr, cur = line, []
        else:
            cur.append(line.strip())
    return "".join(cur)


def jobs(root):
    for tax in sorted(os.listdir(root)):
        td = os.path.join(root, tax)
        if not os.path.isdir(td):
            continue
        for feat in sorted(os.listdir(td)):
            fd = os.path.join(td, feat)
            if not os.path.isdir(fd):
                continue
            for p in sorted(glob.glob(os.path.join(fd, "pssms", "*.pssm"))):
                num = re.sub(r"\.pssm$", "", os.path.basename(p)).split(".")[-1]
                src = next((c for sub in SUBDIRS
                            for c in [os.path.join(fd, sub, num + ".fa")]
                            if os.path.exists(c)), None)
                yield tax, feat, num, p, src


def build(src, lid, title, tmp):
    """Reproduce the pipeline's PSSM for one alignment; returns its text."""
    msa = os.path.join(tmp, "msa.fa")
    with open(msa, "w") as out:
        first = True
        for line in open(src):
            if line.startswith(">"):
                out.write(">%s %s\n" % (lid, title) if first else line)
                first = False
            else:
                out.write(line)
    subj = os.path.join(tmp, "subj.fa")
    ungapped = master_seq(src).replace("-", "").replace(".", "")
    open(subj, "w").write(">%s %s\n%s\n" % (lid, title, ungapped))
    raw = os.path.join(tmp, "raw.pssm")
    subprocess.run(["psiblast", "-subject", subj, "-in_msa", msa, "-out_pssm", raw],
                   capture_output=True)
    if not os.path.exists(raw):
        return None
    text = subprocess.run([sys.executable, os.path.join(HERE, "norm_pssm.py"), raw],
                          capture_output=True, text=True).stdout
    os.remove(raw)
    #  Force the query id back to the PSSM's own.
    #
    #  psiblast does not preserve it: rebuilding cluster 2 yields
    #  `local str "1"` however the subject is labelled. That is a name, not a
    #  score -- across a 30-PSSM sample, 16 reproduced exactly and the other
    #  14 differed in this one line and nothing else -- but leaving it made
    #  rebuild_pssms report 389 of 808 Hepaciviridae PSSMs stale when every
    #  score was identical.
    #
    #  It also matters for --write, not only for the comparison: writing
    #  psiblast's id would silently renumber a rebuilt profile's query, so the
    #  substitution has to happen here, in what gets written, rather than
    #  being neutralised at compare time.
    if lid is not None:
        text = re.sub(r'local str "[^"]*"', 'local str "%s"' % lid, text, count=1)
    return text


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".",
                    help="module working directory holding Alignments/ (default: cwd)")
    ap.add_argument("--write", action="store_true",
                    help="rebuild stale PSSMs in place (originals kept as .prev)")
    args = ap.parse_args()

    root = os.path.join(os.path.abspath(args.workdir), "Alignments")
    if not os.path.isdir(root):
        sys.exit("no Alignments/ under %s -- pass --workdir" % args.workdir)

    tmp = os.path.join(HERE, ".rebuild_tmp")
    os.makedirs(tmp, exist_ok=True)
    ok, stale, failed = 0, [], []
    for tax, feat, num, p, src in jobs(root):
        lid, title = pssm_meta(p)
        if src is None or lid is None:
            failed.append((tax, feat, num, "no source alignment" if src is None else "no query id"))
            continue
        text = build(src, lid, title, tmp)
        if text is None:
            failed.append((tax, feat, num, "psiblast produced nothing"))
        #  Normalise BOTH sides. Comparing a normalised rebuild against the raw
        #  stored file makes the comparison sensitive to anything norm_pssm
        #  strips, which is the opposite of the point. The `descr` block is
        #  emitted by some pipeline runs and not others and holds no scores;
        #  leaving it on one side of the comparison reported 399 of 808
        #  Hepaciviridae PSSMs stale when every score was identical.
        elif _norm(open(p, errors="replace").read()) == _norm(text):
            ok += 1
        else:
            stale.append((tax, feat, num, p, src, text))
    shutil.rmtree(tmp, ignore_errors=True)

    print("%d PSSMs reproduce from their current alignment" % ok)
    print("%d are STALE (alignment edited since the PSSM was built)" % len(stale))
    for tax, feat, num, p, src, _ in stale:
        print("    %-20s %-9s %-6s  src %s" % (tax, feat, num, os.path.basename(os.path.dirname(src))))
    for f in failed:
        print("    FAILED: %s" % (f,))

    if args.write and stale:
        for tax, feat, num, p, src, text in stale:
            shutil.copy2(p, p + ".prev")
            open(p, "w").write(text)
            print("  rebuilt %s" % p)
    elif stale:
        print("\nRe-run with --write to rebuild them.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
