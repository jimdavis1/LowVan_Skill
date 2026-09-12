#!/usr/bin/env python3
"""Turn SignalP predictions into verified G_SP / G_MAT alignments and PSSMs.

A precursor's mature chain is not a prediction problem once the cleavage site is
known -- it is a coordinate operation on an alignment you already curated. The
hard part is deciding whether a site predicted from one consensus sequence
actually applies to every member of the cluster. This does that check and
refuses the clusters that fail it.

    # 1. run SignalP on the staged consensus N-termini (see --help-install)
    signalp6 --fastafile Alignments/G_SIGNALP_QUERIES.faa --organism eukarya \
             --output_dir signalp_out --format none --mode fast

    # 2. apply, verify, report
    python3 apply_signalp.py --workdir . --signalp signalp_out/prediction_results.txt

    # 3. once the report looks right
    python3 apply_signalp.py --workdir . --signalp ... --write

The verification, which is the whole point:

  * SignalP must call it a signal peptide at or above --min-prob.
  * The consensus this tool rebuilds from the alignment must equal the sequence
    SignalP was given. If the alignments have been re-clustered since the query
    file was written, cluster 3 is no longer the same cluster 3, and a cut taken
    on the old consensus would land in the wrong place. Mismatch = refuse.
  * The residue immediately after the cut must be conserved across the cluster
    at or above --min-cons. A signal peptide has almost no conserved sequence,
    but the mature N-terminus does, because it is a protease recognition site.
    This is the same standard the Alpharhabdovirinae G_MAT PROVENANCE used --
    the RABV mature motif FP[ILM]YTIP present at the same column in 100% of
    sequences in six of seven clusters.
  * Every member must actually have residues on both sides of the cut, so a
    fragment cannot contribute a mature chain it does not have.

Clusters that fail are listed with the reason and left unbuilt. That is the
correct outcome: an unbuilt feature is a statement about the evidence, and it
is better than a profile cut at a guessed position.
"""

import argparse
import collections
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
INSTALL = """
SignalP 6.0 is academic-licensed: the download is behind a form you have to
submit yourself (name, institution, and acceptance of the licence terms).

  1. https://services.healthtech.dtu.dk/services/SignalP-6.0/  ->  Downloads
  2. choose 6.0 "fast", submit the form, download signalp-6.0h.fast.tar.gz
  3. tar -xzf signalp-6.0h.fast.tar.gz
     pip install signalp-6-package/
     python -c "import signalp; print(signalp.__file__)"

The web server at the same address accepts pasted FASTA and reports the same
cleavage sites, which is enough for 20 sequences if you would rather not install.
Save its table and pass it with --signalp.
"""


def read_fasta(path):
    h, s = None, []
    for line in open(path, errors="replace"):
        if line.startswith(">"):
            if h is not None:
                yield h, "".join(s)
            h, s = line[1:].strip(), []
        else:
            s.append(line.strip())
    if h is not None:
        yield h, "".join(s)


def consensus_with_columns(rows, length=60):
    """Column-wise consensus plus the alignment column each residue came from.

    Returns (consensus_string, [column index per consensus residue]). Gap-majority
    columns are skipped, which is what makes the consensus ungapped and gives the
    mapping back from a SignalP position to a column in the real alignment.
    """
    if not rows:
        return "", []
    width = len(rows[0][1])
    cons, cols = [], []
    for c in range(width):
        col = [r[1][c] for r in rows if c < len(r[1])]
        counts = collections.Counter(x for x in col if x not in "-.")
        if not counts or len(counts) == 0:
            continue
        if sum(counts.values()) * 2 < len(col):      # gap-majority column
            continue
        cons.append(counts.most_common(1)[0][0])
        cols.append(c)
        if len(cons) >= length:
            break
    return "".join(cons), cols


def conservation(rows, col):
    """(modal residue, fraction carrying it, n occupied) at one alignment column."""
    col_chars = [r[1][col] for r in rows if col < len(r[1])]
    occ = [x for x in col_chars if x not in "-."]
    if not occ:
        return "-", 0.0, 0
    res, n = collections.Counter(occ).most_common(1)[0]
    return res, n / len(occ), len(occ)


def parse_signalp(path):
    """SignalP 6 prediction_results.txt, or the web server's table."""
    out = {}
    for line in open(path, errors="replace"):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        sid = f[0].split()[0]
        m = re.search(r"CS pos:\s*(\d+)-(\d+)", line)
        if not m:
            continue
        prob = 0.0
        pm = re.search(r"Pr:\s*([\d.]+)", line)
        if pm:
            prob = float(pm.group(1))
        pred = f[1].strip()
        out[sid] = {"cut": int(m.group(1)), "prob": prob, "pred": pred}
    return out



def motif_predictions(W, args):
    """Locate a mature N-terminal motif and turn it into cleavage predictions.

    A measured cleavage site is usually recorded as a column number, and a
    column number is only true of one clustering. Re-cluster the parent and it
    silently points somewhere else -- the Alpharhabdovirinae G_MAT derivation
    named ten clusters (1_1, 1_2, 1_3, 3_1..3_4, 9, 11, 14) and after one
    rebuild five no longer existed and the other five had different membership.

    The motif is the same fact in a form that survives. FP[ILM]YTIP re-anchored
    at column 20 across every cluster that carries it, covering more sequences
    than the derivation it replaces, with nothing to update by hand.
    """
    import collections as _c
    rx = re.compile(args.motif)
    out = {}
    for base in sorted(glob.glob(os.path.join(W, "Alignments", "*", args.parent))):
        module = base.split(os.sep)[-2]
        for sub in ("corrected_alis", "reclustered_alis"):
            for fa in sorted(glob.glob(os.path.join(base, sub, "*.fa"))):
                rows = list(read_fasta(fa))
                if len(rows) < args.min_seqs:
                    continue
                cols = []
                for _h, s in rows:
                    m = rx.search(s.replace("-", ""))
                    if not m:
                        continue
                    u = 0
                    for i, ch in enumerate(s):
                        if ch not in "-.":
                            if u == m.start():
                                cols.append(i)
                                break
                            u += 1
                if not cols:
                    continue
                col, n = _c.Counter(cols).most_common(1)[0]
                if n < args.motif_frac * len(rows):
                    continue
                #  express the hit as a consensus position so the shared
                #  verification path below applies unchanged
                cons, cmap = consensus_with_columns(rows, length=max(col + 2, 60))
                if col not in cmap:
                    continue
                out["%s|%s" % (module, os.path.basename(fa)[:-3])] = {
                    "cut": cmap.index(col), "prob": n / len(rows),
                    "pred": "SP(motif)"}
    return out



def self_recall(pssm, rows, tmp):
    """Lowest bitscore this profile gives any of its own training sequences.

    Composition-based statistics halve the score of a profile this short, so
    they are off: tblastn against a genome is the search the annotator runs,
    and it does not apply them either. Measuring with CBS on is what made
    signal peptides look undetectable when they are not.
    """
    seqs = [s.replace("-", "") for _h, s in rows]
    db = os.path.join(tmp, "self.faa")
    with open(db, "w") as fh:
        for i, s in enumerate(seqs):
            fh.write(">s%d\n%s\n" % (i, s))
    subprocess.run(["makeblastdb", "-in", db, "-dbtype", "prot",
                    "-out", os.path.join(tmp, "selfdb")], capture_output=True)
    r = subprocess.run(
        ["psiblast", "-in_pssm", pssm, "-db", os.path.join(tmp, "selfdb"),
         "-comp_based_stats", "0", "-outfmt", "6 sseqid bitscore",
         "-evalue", "1000", "-max_target_seqs", str(len(seqs) + 10)],
        capture_output=True, text=True)
    best = {}
    for line in r.stdout.splitlines():
        a, b = line.split("\t")[:2]
        best[a] = max(best.get(a, 0.0), float(b))
    if len(best) < len(seqs):
        return 0.0, len(seqs) - len(best)
    v = sorted(best.values())
    return v[0], 0


def build_pssm(rows, title, outfile, tmp):
    """Identical code path to build_sp_pssms.py: psiblast over the MSA as given."""
    msa = os.path.join(tmp, "m.msa")
    with open(msa, "w") as fh:
        fh.write(">1 %s\n%s\n" % (title, rows[0][1]))
        for i, (_h, s) in enumerate(rows[1:]):
            fh.write(">%d\n%s\n" % (i + 2, s))
    subj = os.path.join(tmp, "m.subj")
    with open(subj, "w") as fh:
        fh.write(">1 %s\n%s\n" % (title, rows[0][1].replace("-", "")))
    raw = os.path.join(tmp, "m.raw")
    subprocess.run(["psiblast", "-subject", subj, "-in_msa", msa, "-out_pssm", raw],
                   capture_output=True)
    if not os.path.exists(raw):
        return False
    norm = os.path.join(HERE, "norm_pssm.py")
    text = subprocess.run([sys.executable, norm, raw],
                          capture_output=True, text=True).stdout
    os.remove(raw)
    with open(outfile, "w") as fh:
        fh.write(text)
    return True


def main():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=INSTALL)
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--signalp", help="SignalP prediction_results.txt")
    ap.add_argument("--queries", default=None,
                    help="the FASTA SignalP was run on "
                         "(default <workdir>/Alignments/G_SIGNALP_QUERIES.faa)")
    ap.add_argument("--parent", default="G")
    ap.add_argument("--min-prob", type=float, default=0.70)
    ap.add_argument("--min-cons", type=float, default=0.90,
                    help="fraction of the cluster that must carry the modal "
                         "residue at the first mature column (default 0.90)")
    ap.add_argument("--min-seqs", type=int, default=3)
    ap.add_argument("--motif", default=None,
                    help="regex for the MATURE N-terminus. Cuts every cluster of "
                         "--parent where it lands at a consistent column, with no "
                         "SignalP needed. This is the durable form of a measured "
                         "cleavage site: a column number dies at the next "
                         "re-clustering, a motif does not.")
    ap.add_argument("--motif-frac", type=float, default=0.90,
                    help="fraction of a cluster that must carry the motif at the "
                         "modal column (default 0.90)")
    ap.add_argument("--motif-evidence", default="",
                    help="one line naming the measurement the motif stands for, "
                         "written into PROVENANCE")
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--help-install", action="store_true")
    ap.add_argument("--stage-queries", metavar="FASTA", default=None,
                    help="write the consensus N-terminus of every cluster of "
                         "--parent to FASTA, ready for SignalP, and exit. Always "
                         "re-stage after a rebuild: the ids carry cluster numbers, "
                         "and a prediction keyed to a cluster that has since been "
                         "renumbered would be applied to the wrong sequences.")
    ap.add_argument("--stage-len", type=int, default=60)
    ap.add_argument("--self-cutoff", type=float, default=None, metavar="BITS",
                    help="after building, drop any signal-peptide profile that "
                         "cannot score every one of its own training sequences "
                         "at or above BITS. Pass the feature's bit_cutoff. A "
                         "profile that misses its own members would ship a "
                         "feature the annotator can never call -- the reason "
                         "Dichorhavirus G_SP was withdrawn in Sep 2026.")
    ap.add_argument("--corroborate", default=None, metavar="REGEX",
                    help="regex for a mature N-terminus known from direct "
                         "measurement. Clusters matching it are marked in "
                         "PROVENANCE as corroborated, so a reader can tell a "
                         "prediction that agrees with an experiment from a "
                         "prediction standing on its own.")
    args = ap.parse_args()

    if args.stage_queries:
        W = os.path.abspath(args.workdir)
        n = 0
        with open(args.stage_queries, "w") as fh:
            for base in sorted(glob.glob(os.path.join(W, "Alignments", "*", args.parent))):
                module = base.split(os.sep)[-2]
                for sub in ("corrected_alis", "reclustered_alis"):
                    for fa in sorted(glob.glob(os.path.join(base, sub, "*.fa"))):
                        rows = list(read_fasta(fa))
                        if len(rows) < args.min_seqs:
                            continue
                        cons, _c = consensus_with_columns(rows, length=args.stage_len)
                        if not cons:
                            continue
                        fh.write(">%s|%s n=%d\n%s\n"
                                 % (module, os.path.basename(fa)[:-3], len(rows), cons))
                        n += 1
        print("  staged %d consensus N-termini -> %s" % (n, args.stage_queries))
        return 0

    if args.help_install or (not args.signalp and not args.motif):
        print(INSTALL)
        return 0

    W = os.path.abspath(args.workdir)
    qpath = args.queries or os.path.join(W, "Alignments", "G_SIGNALP_QUERIES.faa")
    queries = {h.split()[0]: s for h, s in read_fasta(qpath)} if os.path.exists(qpath) else {}
    preds = parse_signalp(args.signalp) if args.signalp else {}
    if args.motif:
        preds.update(motif_predictions(W, args))
    print("  %d prediction(s), %d staged query sequence(s)\n" % (len(preds), len(queries)))

    built, refused = [], []
    for sid, p in sorted(preds.items()):
        if "|" not in sid:
            refused.append((sid, "id is not <Module>|<cluster>")); continue
        module, clus = sid.split("|", 1)
        base = os.path.join(W, "Alignments", module, args.parent)
        cand = [os.path.join(base, d, clus + ".fa")
                for d in ("corrected_alis", "reclustered_alis")]
        ali = next((c for c in cand if os.path.exists(c)), None)
        if not ali:
            refused.append((sid, "no alignment %s/{corrected,reclustered}_alis/%s.fa"
                            % (args.parent, clus))); continue
        rows = list(read_fasta(ali))
        if len(rows) < args.min_seqs:
            refused.append((sid, "only %d sequence(s)" % len(rows))); continue
        if "SP" not in p["pred"].upper():
            refused.append((sid, "SignalP calls it %s" % p["pred"])); continue
        if p["prob"] < args.min_prob:
            refused.append((sid, "probability %.2f below %.2f" % (p["prob"], args.min_prob)))
            continue

        cons, cols = consensus_with_columns(rows, length=60)
        q = None if p["pred"] == "SP(motif)" else queries.get(sid)
        if q and cons[:len(q)] != q:
            refused.append((sid, "consensus no longer matches the query SignalP saw "
                                 "-- clusters were renumbered, re-stage and re-run"))
            continue
        cut = p["cut"]                      # last residue of the signal peptide
        if cut >= len(cols):
            refused.append((sid, "cut %d beyond the %d consensus residues" % (cut, len(cols))))
            continue
        cut_col = cols[cut - 1]             # alignment column of that last residue
        first_mat_col = cols[cut]           # alignment column of the first mature residue
        res, frac, occ = conservation(rows, first_mat_col)
        if frac < args.min_cons:
            refused.append((sid, "first mature residue %s conserved in only %.0f%% "
                                 "(%d/%d)" % (res, 100 * frac, int(frac * occ), occ)))
            continue
        both = sum(1 for _h, s in rows
                   if s[:cut_col + 1].strip("-.") and s[first_mat_col:].strip("-."))
        if both < len(rows):
            refused.append((sid, "%d of %d sequences lack residues on both sides of the cut"
                            % (len(rows) - both, len(rows))))
            continue

        sp_rows = [(h, s[:cut_col + 1]) for h, s in rows]
        mat_rows = [(h, s[first_mat_col:]) for h, s in rows]
        sp_len = sorted(len(s.replace("-", "")) for _h, s in sp_rows)
        mat_len = sorted(len(s.replace("-", "")) for _h, s in mat_rows)
        built.append(dict(sid=sid, module=module, clus=clus, n=len(rows),
                          cut=cut, prob=p["prob"], res=res, frac=frac,
                          sp_rows=sp_rows, mat_rows=mat_rows,
                          sp_len=(sp_len[0], sp_len[-1]),
                          mat_len=(mat_len[0], mat_len[-1]),
                          motif="".join(cons[cut:cut + 8]),
                          corrob=bool(args.corroborate and
                                      re.match(args.corroborate, "".join(cons[cut:cut + 8])))))

    print("  BUILDABLE (%d)" % len(built))
    print("  %-26s %5s %5s %5s  %-9s %-10s %s"
          % ("cluster", "n", "cut", "Pr", "SP aa", "mature aa", "mature motif / conservation"))
    for b in built:
        print("  %-26s %5d %5d %5.2f  %-9s %-10s %-9s %s at %.0f%%"
              % (b["sid"], b["n"], b["cut"], b["prob"],
                 "%d-%d" % b["sp_len"], "%d-%d" % b["mat_len"],
                 b["motif"], b["res"], 100 * b["frac"]))
    if refused:
        print("\n  REFUSED (%d) -- left unbuilt on purpose" % len(refused))
        for sid, why in refused:
            print("    %-26s %s" % (sid, why))

    if not args.write:
        print("\n  nothing written (rerun with --write)")
        return 0

    dropped = []
    by_mod = collections.defaultdict(list)
    for b in built:
        by_mod[b["module"]].append(b)
    for module, bs in sorted(by_mod.items()):
        for feat, key in ((args.parent + "_SP", "sp_rows"), (args.parent + "_MAT", "mat_rows")):
            d = os.path.join(W, "Alignments", module, feat)
            #  Clear before writing. These directories are wholly derived from
            #  the current parent clustering, so anything already in them came
            #  from a clustering that no longer exists -- exactly how the old
            #  G_MAT ended up citing five alignments that had been renumbered
            #  away. A derived directory must describe exactly one derivation.
            for sub in ("corrected_alis", "reclustered_alis", "pssms"):
                shutil.rmtree(os.path.join(d, sub), ignore_errors=True)
            for sub in ("corrected_alis", "pssms"):
                os.makedirs(os.path.join(d, sub), exist_ok=True)
            tmp = tempfile.mkdtemp(prefix="sp.")
            try:
                for b in bs:
                    rows = b[key]
                    fa = os.path.join(d, "corrected_alis", b["clus"] + ".fa")
                    with open(fa, "w") as fh:
                        for h, s in rows:
                            fh.write(">%s\n%s\n" % (h, s))
                    title = "%s.%s.%s" % (module, feat, b["clus"])
                    out = os.path.join(d, "pssms", title + ".pssm")
                    build_pssm(rows, title, out, tmp)
                    #  A short profile that cannot find its own training
                    #  sequences is a feature the annotator can never call.
                    if (args.self_cutoff and feat.endswith("_SP")
                            and os.path.exists(out)):
                        lo, missing = self_recall(out, rows, tmp)
                        if missing or lo < args.self_cutoff:
                            os.remove(out)
                            os.remove(fa)
                            dropped.append((module, b["clus"], lo, missing,
                                            len(rows)))
            finally:
                shutil.rmtree(tmp, ignore_errors=True)
            #  Report only what survived. A PROVENANCE table listing profiles
            #  that were dropped for failing self-recall claims coverage the
            #  module does not have.
            gone = {x[1] for x in dropped if x[0] == module} if feat.endswith("_SP") else set()
            kept = [b for b in bs if b["clus"] not in gone]
            with open(os.path.join(d, "PROVENANCE"), "w") as fh:
                fh.write("%s for %s -- %d of the %s alignments.\n\n"
                         % (feat, module, len(kept), args.parent))
                if args.motif:
                    fh.write("Cut where the mature N-terminal motif /%s/ falls, which is a\n"
                             "measured site expressed so that it re-locates itself after any\n"
                             "re-clustering of %s.\n%s\n\n"
                             % (args.motif, args.parent,
                                args.motif_evidence or
                                "  (record the measurement this motif stands for "
                                "with --motif-evidence)"))
                fh.write("Cut at cleavage sites, each verified against the\n"
                         "alignment before cutting rather than trusted from the consensus:\n"
                         "the first mature residue had to be the modal residue in at least\n"
                         "%.0f%% of the cluster, and every member had to carry residues on\n"
                         "both sides of the cut.\n\n" % (100 * args.min_cons))
                fh.write("  %-10s %5s %5s %5s  %-9s %-10s %-9s %s\n"
                         % ("cluster", "n", "cut", "Pr", "SP aa", "mature aa",
                            "motif", "evidence"))
                for b in kept:
                    fh.write("  %-10s %5d %5d %5.2f  %-9s %-10s %-9s %s\n"
                             % (b["clus"], b["n"], b["cut"], b["prob"],
                                "%d-%d" % b["sp_len"], "%d-%d" % b["mat_len"],
                                b["motif"],
                                "prediction + measured site"
                                if b.get("corrob") else "prediction"))
                nc = sum(1 for b in kept if b.get("corrob"))
                if nc:
                    fh.write("\n%d of %d clusters carry a mature N-terminus matching a\n"
                             "directly measured cleavage site, so for those the prediction\n"
                             "is corroborated by experiment rather than standing alone.\n"
                             % (nc, len(kept)))
                if gone:
                    fh.write("\nWITHHELD -- these clusters were cut successfully but their\n"
                             "signal-peptide profile could not score its own training\n"
                             "sequences at the feature's %.0f-bit cutoff, so the profile\n"
                             "would never fire. The mature chain is still built.\n"
                             % (args.self_cutoff or 0))
                    for _m, _c, _lo, _mi, _n in sorted(x for x in dropped if x[0] == module):
                        fh.write("  cluster %-8s %d seqs, lowest self-score %.1f bits\n"
                                 % (_c, _n, _lo))
                fh.write("\nNOT a pipeline product: there is no %s collection, because the\n"
                         "split is a coordinate operation on curated %s rather than a\n"
                         "binning rule.\n" % (feat, args.parent))
            _d = [x for x in dropped if x[0] == module]
            if _d and feat.endswith("_SP"):
                for _m, _c, _lo, _mi, _n in _d:
                    print("  dropped %s/%s %s: self-recall %.1f bits over %d seqs, "
                          "below the %.0f-bit cutoff"
                          % (_m, feat, _c, _lo, _n, args.self_cutoff))
            print("  wrote %s/%s: %d alignment(s), %d PSSM(s)"
                  % (module, feat,
                     len(glob.glob(os.path.join(d, "corrected_alis", "*.fa"))),
                     len(glob.glob(os.path.join(d, "pssms", "*.pssm")))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
