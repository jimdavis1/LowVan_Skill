#!/usr/bin/env python3
"""Build a module's PSSMs from its alignments.

`rebuild_pssms.py` refreshes profiles that already exist -- it walks the pssms
directory and rebuilds each entry from its alignment, so it can detect drift but
cannot create a profile that is absent. That distinction matters when a module
is distributed as alignments only: the alignments are 2 MB of diffable text and
the profiles are 26 MB of derived matrices, so shipping the alignments is right,
but only if the profiles can actually be regenerated from them.

This is that step. It walks corrected_alis/ and reclustered_alis/ and writes one
profile per alignment, in the same format rebuild_pssms.py produces, so the two
agree and a freshly generated set reports 0 stale.

    python3 pssms_from_alignments.py --workdir . --module Tobamovirus
    python3 pssms_from_alignments.py --workdir . --module Tobamovirus --write
"""
import argparse, glob, os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_DIRS = ("corrected_alis", "reclustered_alis")


def read_fasta(p):
    h, s = None, []
    for line in open(p, errors="replace"):
        if line.startswith(">"):
            if h is not None: yield h, "".join(s)
            h, s = line[1:].strip(), []
        else: s.append(line.strip())
    if h is not None: yield h, "".join(s)


def build(ali_path, title, out_path, tmp):
    rows = list(read_fasta(ali_path))
    if len(rows) < 2:
        return None, "fewer than 2 rows"
    master = rows[0][1]
    msa = os.path.join(tmp, "m.msa")
    with open(msa, "w") as f:
        f.write(">1 %s\n%s\n" % (title, master))
        for i, (_h, s) in enumerate(rows[1:]):
            f.write(">%d\n%s\n" % (i + 2, s))
    subj = os.path.join(tmp, "m.subj")
    open(subj, "w").write(">1 %s\n%s\n" % (title, master.replace("-", "").replace(".", "")))
    raw = os.path.join(tmp, "m.raw")
    if os.path.exists(raw): os.remove(raw)
    subprocess.run(["psiblast", "-subject", subj, "-in_msa", msa, "-out_pssm", raw],
                   capture_output=True)
    if not os.path.exists(raw):
        return None, "psiblast produced nothing"
    txt = subprocess.run([sys.executable, os.path.join(HERE, "norm_pssm.py"), raw],
                         capture_output=True, text=True).stdout
    os.remove(raw)
    return txt, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--module", required=True)
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()

    base = os.path.join(a.workdir, "Alignments", a.module)
    if not os.path.isdir(base):
        sys.exit("no Alignments/%s in %s" % (a.module, os.path.abspath(a.workdir)))
    tmp = tempfile.mkdtemp()
    made = skipped = failed = 0
    for feat in sorted(os.listdir(base)):
        fd = os.path.join(base, feat)
        if not os.path.isdir(fd): continue
        alis = []
        for d in SRC_DIRS:
            alis += sorted(glob.glob(os.path.join(fd, d, "*.fa")))
        if not alis: continue
        pdir = os.path.join(fd, "pssms")
        n_ok = n_bad = 0
        for p in alis:
            cid = os.path.basename(p)[:-3]
            title = "%s.%s.%s" % (a.module, feat, cid)
            out = os.path.join(pdir, title + ".pssm")
            if os.path.exists(out) and not a.write:
                skipped += 1; continue
            txt, err = build(p, title, out, tmp)
            if txt is None:
                n_bad += 1; failed += 1
                print("    %-28s FAILED: %s" % (title, err))
                continue
            if a.write:
                os.makedirs(pdir, exist_ok=True)
                open(out, "w").write(txt)
            n_ok += 1; made += 1
        print("  %-14s %3d alignment(s) -> %d profile(s)%s"
              % (feat, len(alis), n_ok, "  %d failed" % n_bad if n_bad else ""))
    print("\n%d profile(s) %s, %d already present, %d failed"
          % (made, "written" if a.write else "would be written", skipped, failed))
    if not a.write:
        print("(dry run -- rerun with --write)")

main()
