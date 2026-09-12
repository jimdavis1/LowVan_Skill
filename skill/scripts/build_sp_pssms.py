#!/usr/bin/env python3
"""Build per-cluster signal-peptide PSSMs from a parent CDS and its mature peptide.

The signal peptide is whatever the mature peptide is missing from the front of
its precursor, so if both alignments exist the cleavage site is already known
per sequence and nothing needs predicting. For each cluster where <PARENT> and
<PARENT>_MAT share a name, this subtracts the mature sequence from the precursor,
aligns the resulting signal peptides, and builds a PSSM.

**One profile per cluster, never one for the whole feature.** A signal peptide
is ~19 residues of hydrophobic core with almost no conserved sequence, so a
profile only works inside the cluster it came from. Measured on the
Alpharhabdovirinae lyssavirus cluster: 47 bits on rabies, 17 bits on Irkut virus
and 16 on Australian bat lyssavirus -- both still lyssaviruses. A single
family-wide signal-peptide profile is genuinely hopeless; a per-cluster one is
ordinary.

Scope: this derives the **leading** product of a precursor -- the part the
mature peptide is missing from its front. A precursor cleaved more than once
(GPC into signal peptide + Gn + Gc, HA0 into HA1 + HA2) needs each internal
product bounded by *both* its neighbours, which this does not do; align those by
hand from the published cleavage sites and set both extensions to 0.

    python3 build_sp_pssms.py --workdir . --module Alpharhabdovirinae
    python3 build_sp_pssms.py --workdir . --module Alpharhabdovirinae --write

Setting bit_cutoff: pass --decoys with a FASTA of genomes and the tool measures
the highest score each profile reaches anywhere in them. The cutoff has to sit
above that noise ceiling and below the weakest true signal. Measured on
Alpharhabdovirinae over a seven-genome panel: the worst off-target hit across
all ten profiles was 19.9 bits, and true signal peptides score 45-47, so 30 is
comfortably inside the window. Without a decoy panel you are guessing.

Without --write it reports what it would build and how well each profile scores
against its own signal peptides. Requires mafft and psiblast.
"""

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SUBS = ("corrected_alis", "reclustered_alis", "truncated_alis")


def read_fasta(path):
    hdr, seq = None, []
    for line in open(path, errors="replace"):
        if line.startswith(">"):
            if hdr is not None:
                yield hdr, "".join(seq)
            hdr, seq = line[1:].rstrip("\n"), []
        else:
            seq.append(line.strip())
    if hdr is not None:
        yield hdr, "".join(seq)


def ungapped(workdir, module, feat):
    """{cluster: {seq_id: ungapped sequence}} across every alignment subdir."""
    out = {}
    for sub in SUBS:
        for p in glob.glob(os.path.join(workdir, "Alignments", module, feat, sub, "*.fa")):
            clus = os.path.basename(p)[:-3]
            d = {}
            for hdr, s in read_fasta(p):
                ident = re.sub(r"/\d+-\d+", "", hdr.split()[0])
                d[ident] = s.replace("-", "").replace(".", "").upper()
            if d:
                out[clus] = d
    return out


def build_pssm(seqs, title, tmp, outfile):
    """mafft + psiblast, normalised to the pipeline's PSSM format."""
    fa = os.path.join(tmp, "sp.fa")
    with open(fa, "w") as fh:
        for i, s in enumerate(seqs):
            fh.write(">%d\n%s\n" % (i + 1, s))
    ali = os.path.join(tmp, "sp.ali")
    with open(ali, "w") as fh:
        r = subprocess.run(["mafft", "--quiet", "--auto", fa],
                           stdout=fh, stderr=subprocess.DEVNULL)
    rows = list(read_fasta(ali))
    if not rows:
        return None
    msa = os.path.join(tmp, "sp.msa")
    with open(msa, "w") as fh:
        fh.write(">1 %s\n%s\n" % (title, rows[0][1]))
        for i, (_h, s) in enumerate(rows[1:]):
            fh.write(">%d\n%s\n" % (i + 2, s))
    subj = os.path.join(tmp, "subj.fa")
    with open(subj, "w") as fh:
        fh.write(">1 %s\n%s\n" % (title, rows[0][1].replace("-", "")))
    raw = os.path.join(tmp, "raw.pssm")
    subprocess.run(["psiblast", "-subject", subj, "-in_msa", msa, "-out_pssm", raw],
                   capture_output=True)
    if not os.path.exists(raw):
        return None
    text = subprocess.run([sys.executable, os.path.join(HERE, "norm_pssm.py"), raw],
                          capture_output=True, text=True).stdout
    os.remove(raw)
    with open(outfile, "w") as fh:
        fh.write(text)
    return len(rows[0][1])


def self_score(pssm, seqs, tmp):
    """Best bitscore of the profile against its own signal peptides, CBS off.

    Composition-based statistics halve the score of an alignment this short --
    22 bits vs 43 on the same data -- which is what makes a signal peptide look
    undetectable when it is not. tblastn against a genome is the search that
    matters, so measure without CBS to see the real signal.
    """
    db = os.path.join(tmp, "self.faa")
    with open(db, "w") as fh:
        for i, s in enumerate(seqs):
            fh.write(">s%d\n%s\n" % (i, s))
    subprocess.run(["makeblastdb", "-in", db, "-dbtype", "prot",
                    "-out", os.path.join(tmp, "selfdb")], capture_output=True)
    r = subprocess.run(
        ["psiblast", "-in_pssm", pssm, "-db", os.path.join(tmp, "selfdb"),
         "-comp_based_stats", "0", "-outfmt", "6 sseqid bitscore",
         "-evalue", "100", "-max_target_seqs", str(len(seqs) + 10)],
        capture_output=True, text=True)
    best = {}
    for line in r.stdout.splitlines():
        s, b = line.split("\t")[:2]
        best[s] = max(best.get(s, 0.0), float(b))
    if not best:
        return 0, 0.0
    v = sorted(best.values())
    return len(v), v[len(v) // 2]


def noise_ceiling(pssm, decoys_db, tmp):
    """Highest score this profile reaches anywhere in the decoy genomes.

    A cutoff must sit above this. tblastn is the search the annotator actually
    runs, so measure with tblastn against DNA rather than psiblast against
    protein -- the two score a 19-residue profile very differently.
    """
    r = subprocess.run(
        ["tblastn", "-in_pssm", pssm, "-db", decoys_db, "-outfmt", "6 bitscore",
         "-evalue", "1000", "-max_target_seqs", "500"],
        capture_output=True, text=True)
    v = sorted((float(x) for x in r.stdout.split()), reverse=True)
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--module", required=True)
    ap.add_argument("--parent", default="G", help="parent CDS feature key (default G)")
    ap.add_argument("--mature", help="mature feature key (default <parent>_MAT)")
    ap.add_argument("--out", help="signal-peptide feature key (default <parent>_SP)")
    ap.add_argument("--min-seqs", type=int, default=2)
    ap.add_argument("--decoys", help="FASTA of genomes; measures the off-target "
                                     "noise ceiling so a bit_cutoff can be chosen")
    ap.add_argument("--write", action="store_true", help="write alignments and PSSMs")
    args = ap.parse_args()

    workdir = os.path.abspath(args.workdir)
    mature = args.mature or (args.parent + "_MAT")
    spkey = args.out or (args.parent + "_SP")

    decoys_db = None
    if args.decoys:
        decoys_db = os.path.join(tempfile.mkdtemp(prefix="spdecoy."), "d")
        subprocess.run(["makeblastdb", "-in", os.path.abspath(args.decoys),
                        "-dbtype", "nucl", "-out", decoys_db], capture_output=True)

    par = ungapped(workdir, args.module, args.parent)
    mat = ungapped(workdir, args.module, mature)
    shared = sorted(set(par) & set(mat))
    if not shared:
        sys.exit("no cluster has both %s and %s in %s" % (args.parent, mature, args.module))

    print("%s: %d clusters have both %s and %s\n" % (args.module, len(shared), args.parent, mature))
    dest = os.path.join(workdir, "Alignments", args.module, spkey)
    tmp = tempfile.mkdtemp(prefix="sppssm.")
    made, skipped, noise_max = [], [], []
    try:
        for clus in shared:
            P, M = par[clus], mat[clus]
            sps, lens = [], set()
            for ident, pseq in P.items():
                mseq = M.get(ident)
                if not mseq or not pseq.endswith(mseq) or len(pseq) <= len(mseq):
                    continue
                sp = pseq[:len(pseq) - len(mseq)]
                if 8 <= len(sp) <= 60:
                    sps.append(sp)
                    lens.add(len(sp))
            uniq = sorted(set(sps))
            if len(uniq) < args.min_seqs:
                skipped.append((clus, len(sps), len(uniq)))
                print("  %-6s %3d pairs, %2d distinct  -- too few, skipped" % (clus, len(sps), len(uniq)))
                continue
            if args.write:
                os.makedirs(os.path.join(dest, "corrected_alis"), exist_ok=True)
                os.makedirs(os.path.join(dest, "pssms"), exist_ok=True)
            outp = os.path.join(dest, "pssms", "%s.%s.%s.pssm" % (args.module, spkey, clus))
            tgt = outp if args.write else os.path.join(tmp, "t.pssm")
            cols = build_pssm(uniq, "Signal peptide of %s" % args.parent, tmp, tgt)
            n, med = self_score(tgt, uniq, tmp)
            noise = ""
            top = None
            if decoys_db:
                v = noise_ceiling(tgt, decoys_db, tmp)
                # the profile's own genome is usually in the panel, so the top
                # hit is the true positive; the runner-up is the noise ceiling
                top = v[0] if v else 0.0
                second = v[1] if len(v) > 1 else 0.0
                noise = "  decoys: top %.0f, next %.0f" % (top, second)
                noise_max.append(second)
            print("  %-6s %3d pairs, %2d distinct, SP len %-9s -> %2d cols, self-recall %d/%d median %.0f bits%s"
                  % (clus, len(sps), len(uniq),
                     "-".join(str(x) for x in sorted(lens)[:3]) + ("+" if len(lens) > 3 else ""),
                     cols or 0, n, len(uniq), med, noise))
            if args.write:
                fa = os.path.join(dest, "corrected_alis", "%s.fa" % clus)
                shutil.copy(os.path.join(tmp, "sp.ali"), fa)
            made.append((clus, len(uniq), med))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print("\n  %d profile(s) %s, %d skipped"
          % (len(made), "written" if args.write else "would be built", len(skipped)))
    if made:
        lo = min(m[2] for m in made)
        print("  weakest true signal   : %.0f bits (lowest median self-score)" % lo)
        if noise_max:
            hi = max(noise_max)
            print("  noise ceiling         : %.0f bits (worst off-target hit in the decoys)" % hi)
            if hi < lo:
                print("  -> set bit_cutoff between %d and %d; %d is a reasonable middle"
                      % (int(hi) + 1, int(lo), int((hi + lo) / 2)))
            else:
                print("  -> NO SAFE CUTOFF: noise reaches the true signal. These profiles")
                print("     cannot be used as they stand -- split the clusters further.")
        else:
            print("  no decoy panel given, so the noise ceiling is unmeasured and any")
            print("  bit_cutoff is a guess. Rerun with --decoys.")
    if not args.write:
        print("  rerun with --write to build them")
    return 0


if __name__ == "__main__":
    sys.exit(main())
