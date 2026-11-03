#!/usr/bin/env python3
"""Build a profile for a feature whose cluster members are exact duplicates.

THE PROBLEM. fasta-cluster-pssm-2.pl deduplicates, so a cluster whose members
are identical collapses to one sequence, and `psiblast -in_msa` cannot build a
profile from a single row. Barthadenovirus RH4 has eight sequences in the
collection, but they are duplicate pairs across same-species genomes: every
cluster collapses to one and the feature builds nothing. RH6 is the same at
n=3. This is NOT divergence -- lowering the clustering identity does not help
and would be the wrong lever anyway.

CURATOR'S RULING, 3 Nov 2026: "when you know that a transcript is real for a
fact, but only have duplicates, you can keep 2 copies and still build the
pssm." So the duplicate row is retained deliberately, and only for features
whose product is established.

THE EVIDENCE REQUIRED, per feature, from PMID_claude_generated 12504581 (Both
2002, the right-end transcription map of OAdV-7):
  RH4  "Transcript 4 could produce RH4 protein"; one of the four related
       F-box genes named in the paper's title.
  RH6  "RH2 and RH6 proteins could also be produced from transcripts 3 and 2
       ... by reinitiation of translation", and at 24 h p.i. "modified
       splicing produces transcripts 6 and 7 that can only code for ORF RH6".

RH0 and RH3 are deliberately NOT in this list: RH0 "lacks an initiating
methionine" and no RH3-specific transcript was ever identified. They were
removed as features rather than rescued.

A duplicate-row profile is a lookup table for that sequence, not a
generalisation, and it is recorded as such in BUILD_PARAMS.
"""
import os, shutil, subprocess, sys, tempfile

W = os.path.dirname(os.path.abspath(__file__)) or "."
#  feature -> the transcript evidence that licenses a duplicate-row profile
LICENSED = {
 "Barthadenovirus": {
  "RH4": "PMID_claude_generated 12504581: 'Transcript 4 could produce RH4 "
         "protein'; one of the four F-box genes in the paper's title.",
  "RH6": "PMID_claude_generated 12504581: RH6 produced from transcript 2 by "
         "reinitiation, and transcripts 6 and 7 at 24 h p.i. 'can only code "
         "for ORF RH6'.",
 },
}
MI, MC = "0.8", "0.8"      # the settled identity. Not tuned here.


def read_fasta(p):
    out, name, seq = [], None, []
    for l in open(p):
        if l.startswith(">"):
            if name:
                out.append((name, "".join(seq)))
            name, seq = l[1:].strip(), []
        else:
            seq.append(l.strip())
    if name:
        out.append((name, "".join(seq)))
    return out


def main():
    mod = sys.argv[1]
    for key, why in sorted(LICENSED.get(mod, {}).items()):
        src = os.path.join(W, "build", mod, "collections", mod, "%s.fasta" % key)
        dest = os.path.join(W, "build", mod, "Alignments", mod, key)
        seqs = read_fasta(src)
        t = tempfile.mkdtemp(prefix="dup.")
        try:
            # cluster at the settled identity, exactly as the pipeline would
            fa = os.path.join(t, "in.faa")
            with open(fa, "w") as fh:
                for i, (n, s) in enumerate(seqs):
                    fh.write(">%d\n%s\n" % (i, s))
            subprocess.run(["mmseqs", "easy-cluster", fa, t + "/c", t + "/tmp",
                            "--min-seq-id", MI, "-c", MC, "--cov-mode", "0"],
                           capture_output=True)
            grp = {}
            for l in open(t + "/c_cluster.tsv"):
                r, m = l.split()
                grp.setdefault(r, []).append(int(m))
            os.makedirs(os.path.join(dest, "corrected_alis"), exist_ok=True)
            os.makedirs(os.path.join(dest, "pssms"), exist_ok=True)
            for f in os.listdir(os.path.join(dest, "corrected_alis")):
                os.remove(os.path.join(dest, "corrected_alis", f))
            n_ali = 0
            for members in sorted(grp.values(), key=len, reverse=True):
                rows = [seqs[i] for i in members]
                distinct = {s for _n, s in rows}
                if len(distinct) == 1:
                    #  THE DUPLICATE ROW, kept on purpose. Without it psiblast
                    #  has one sequence and refuses.
                    rows = rows[:1] * 2
                n_ali += 1
                af = os.path.join(t, "a%d.faa" % n_ali)
                with open(af, "w") as fh:
                    for j, (nm, s) in enumerate(rows):
                        fh.write(">%s_%d\n%s\n" % (nm.split()[0], j, s))
                out = os.path.join(dest, "corrected_alis", "%d.fa" % n_ali)
                r = subprocess.run(["mafft", "--auto", "--quiet", af],
                                   capture_output=True, text=True)
                open(out, "w").write(r.stdout)
            print("  %-6s %d sequences -> %d alignment(s) written"
                  % (key, len(seqs), n_ali))
            with open(os.path.join(dest, "BUILD_PARAMS"), "w") as fh:
                fh.write("# %s: duplicate-row profile, built by "
                         "build_duplicate_pssms.py\n" % key)
                fh.write("# mi=%s mc=%s (settled, not tuned)\n" % (MI, MC))
                fh.write("# A cluster whose members are exact duplicates keeps "
                         "two copies so psiblast has two rows.\n")
                fh.write("# This is a lookup table for that sequence, not a "
                         "generalisation.\n")
                fh.write("# Licensed by: %s\n" % why)
        finally:
            shutil.rmtree(t, ignore_errors=True)


main()
