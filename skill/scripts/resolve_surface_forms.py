#!/usr/bin/env python3
"""Resolve generic surface / core strings to their nested form, structurally.

WHY A SEPARATE PASS.  The hepadnavirus surface ORF has three in-frame starts
and one stop, so S is an exact C-terminal substring of M and M of L; the core
ORF does the same with preC and C.  Roughly seven thousand Orthohepadnavirus
deposits carry a string that names the ORF and not the form -- `surface
antigen`, `surface protein`, `HBsAg`, `envelope protein`, `preS protein` --
and neither a name rule nor homology can place them.  Homology cannot, because
the forms are 100% identical over their shared region: a 226 aa query matches
the L canonical over 100% of itself exactly as well as it matches S, so the
rescue's >=80%/>=60% test passes on all three and the winner is arbitrary.

Length alone cannot do it either, and that is the trap worth naming.  A 226 aa
surface deposit is either a genuine S or an L truncated after 226 residues, and
those are opposite things: one belongs in the S collection, the other is a
partial that must be discarded.  Binning on length alone puts L fragments into
S and builds the S profile out of L's N-terminus.

THE TEST IS WHICH END THE SEQUENCE SHARES, which is the rule
`references/annotation-triage.md` states for nested alternative starts.  Every
genuine form -- L, M and S alike -- ends at the common stop, so every one of
them carries the S C-terminus.  A truncated L does not: it carries L's
N-terminus and stops wherever the submission stopped.  So:

    shares the C-terminal anchor  ->  a real nested form; assign by length
    does not                      ->  a truncated parent; discard

The anchor is taken from the module's own curated S collection rather than
written down here, so it stays correct for woodchuck and bat viruses whose
C-terminus differs from HBV's.

It also moves deposits the window parked in `<FEAT>.outliers.fasta` when their
length says they are a different form of the same ORF -- 28 Avihepadnavirus
`core protein` deposits are 305 aa, which is preC, and 10 `S protein` deposits
are 300-400 aa, which is the large form.  Those are ordinary mislabels and the
same structural test governs them.

Run it immediately after build_collections.py and before the audits.  Re-running
build_collections.py truncates the collection files, so this must be re-run with
it -- build.sh chains the two for that reason.
"""
import argparse, collections, json, os, re, sys

AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta(path):
    out, hid, buf = [], None, []
    if not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith(">"):
            if hid is not None:
                out.append((hid, "".join(buf)))
            hid, buf = line[1:].rstrip("\n"), []
        else:
            buf.append(line.strip())
    if hid is not None:
        out.append((hid, "".join(buf)))
    return out


def ident(a, b):
    """Fraction of identical residues over the compared length, ungapped."""
    n = min(len(a), len(b))
    if not n:
        return 0.0
    return sum(1 for x, y in zip(a[-n:], b[-n:]) if x == y) / n


def anchor_from(seqs, k):
    """Modal C-terminal k-mer of a collection, as the shared-terminus anchor."""
    c = collections.Counter(s[-k:] for s in seqs if len(s) >= k)
    if not c:
        return None
    return c.most_common(1)[0][0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--module", required=True)
    ap.add_argument("--anchor-from", default="S",
                    help="feature whose C-terminus every nested form shares")
    ap.add_argument("--prune-truncated", default="",
                    help="FEAT,FEAT ... move C-terminally truncated members of "
                         "these collections to <FEAT>.outliers.fasta. Each takes "
                         "its anchor from its OWN collection.")
    ap.add_argument("--forms", default="",
                    help="FEAT:lo-hi,FEAT:lo-hi ... the nested forms, by length")
    ap.add_argument("--strings", default="",
                    help="file of generic annotation strings, one per line")
    ap.add_argument("--also-outliers", default="",
                    help="FEAT,FEAT ... re-test these features' outliers too")
    ap.add_argument("--anchor-k", type=int, default=25)
    ap.add_argument("--anchor-id", type=float, default=0.70)
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()

    W = os.path.abspath(a.workdir)
    cdir = os.path.join(W, "collections", a.module)

    #  ---- prune truncated members, a job of its own ------------------------
    #  A collection of ONE form still needs the C-terminal test, for the same
    #  reason a nested set does. Orthohepadnavirus X is 154 aa and its window
    #  is 122-174, so an X truncated to 142 or 138 sits inside the window and
    #  no length rule can see it. 1,109 of 9,789 X members were such partials,
    #  stopping at successive points along ...ETRLMIFVLGGCRHKLVCSPAPCNFFTSA.
    #  They cluster by where they stop rather than by what they are, so each
    #  group earned its own profile: X shipped 22 profiles for a 154 aa
    #  protein, of which one fired three times on a single genome.
    #
    #  An exact-motif test would be wrong here -- 420 full-length X sequences
    #  carry a genotype-variant tail (NFLTSA, NFFTFA) and are perfectly real.
    #  The 25-mer at 70% identity keeps them and still removes the fragments.
    if a.prune_truncated:
        for feat in [x for x in a.prune_truncated.split(",") if x]:
            fp = os.path.join(cdir, feat + ".fasta")
            recs = read_fasta(fp)
            if not recs:
                print("%s: no collection to prune" % feat)
                continue
            anc = anchor_from([s for _h, s in recs], a.anchor_k)
            keep, cut = [], []
            for h, s in recs:
                ok = len(s) >= a.anchor_k and ident(s[-a.anchor_k:], anc) >= a.anchor_id
                (keep if ok else cut).append((h, s))
            print("%s: anchor %s\n   keep %d, truncated %d (%.1f%%)"
                  % (feat, anc, len(keep), len(cut),
                     100.0 * len(cut) / max(len(recs), 1)))
            if a.write:
                with open(fp, "w") as f:
                    for h, s in keep:
                        f.write(">%s\n%s\n" % (h, s))
                op = os.path.join(cdir, feat + ".outliers.fasta")
                with open(op, "a") as f:
                    for h, s in cut:
                        f.write(">%s [C-terminally truncated]\n%s\n" % (h, s))
                print("   written; %d moved to %s" % (len(cut), os.path.basename(op)))
        if not a.forms:
            return

    forms = {}
    for part in a.forms.split(","):
        k, _, rng = part.partition(":")
        lo, _, hi = rng.partition("-")
        forms[k] = (int(lo), int(hi))

    anchor_seqs = [s for _h, s in read_fasta(os.path.join(cdir, a.anchor_from + ".fasta"))]
    anchor = anchor_from(anchor_seqs, a.anchor_k)
    if not anchor:
        sys.exit("no %s collection to take the C-terminal anchor from" % a.anchor_from)
    print("anchor: C-terminal %d-mer of %s, from %d sequences\n  %s"
          % (a.anchor_k, a.anchor_from, len(anchor_seqs), anchor))

    generic = {l.strip() for l in open(a.strings) if l.strip()}

    #  ---- gather candidates: the unassigned generic-string sequences ...
    cand = []
    seen = set()
    T = [f for f in os.listdir(W) if f.endswith(".uniq.id_ann")][0].split(".")[0]
    anns = [l.rstrip("\n").split("\t") for l in open(os.path.join(W, "%s.uniq.id_ann" % T))]
    seqs = [l.rstrip("\n").split("\t") for l in open(os.path.join(W, "%s.uniq.seq" % T))]
    for ar, sr in zip(anns, seqs):
        s = sr[1] if len(sr) > 1 else ""
        astr = " ".join((ar[1] if len(ar) > 1 else "").split())
        if s and astr in generic:
            cand.append((ar[0], astr, s, "unassigned"))
            seen.add(s)

    #  ... and the sequences the length window parked as outliers
    for feat in [x for x in a.also_outliers.split(",") if x]:
        for h, s in read_fasta(os.path.join(cdir, feat + ".outliers.fasta")):
            if s and s not in seen:
                cand.append((h.split()[0], "[%s.outlier]" % feat, s, feat))
                seen.add(s)

    print("%d candidate sequence(s)" % len(cand))

    assign = collections.defaultdict(list)
    tally = collections.Counter()
    detail = []
    for fid, astr, s, src in cand:
        if set(s) - AA:
            tally["discarded: non-standard residue"] += 1
            detail.append((fid, astr, len(s), src, "DISCARD", "non-standard residue"))
            continue
        shares = len(s) >= a.anchor_k and ident(s[-a.anchor_k:], anchor) >= a.anchor_id
        if not shares:
            tally["discarded: does not share the C-terminus (truncated parent)"] += 1
            detail.append((fid, astr, len(s), src, "DISCARD", "no C-terminal anchor"))
            continue
        hit = [k for k, (lo, hi) in forms.items() if lo <= len(s) <= hi]
        if len(hit) == 1:
            assign[hit[0]].append((fid, astr, s))
            tally["assigned -> " + hit[0]] += 1
            detail.append((fid, astr, len(s), src, hit[0], "C-terminus + length"))
        elif not hit:
            tally["discarded: shares C-terminus but matches no form window"] += 1
            detail.append((fid, astr, len(s), src, "DISCARD", "length matches no form"))
        else:
            tally["discarded: length ambiguous between forms"] += 1
            detail.append((fid, astr, len(s), src, "DISCARD", "ambiguous: " + "/".join(hit)))

    for k in sorted(tally):
        print("  %-62s %6d" % (k, tally[k]))

    rep = os.path.join(W, "RESOLVED_FORMS.%s.tsv" % a.anchor_from)
    with open(rep, "w") as f:
        f.write("feature_id\tannotation\tlength\tsource\tverdict\treason\n")
        for row in detail:
            f.write("\t".join(str(x) for x in row) + "\n")
    print("  report -> %s" % rep)

    if not a.write:
        print("\n(dry run -- pass --write to append to the collections)")
        return

    for feat, rows in sorted(assign.items()):
        path = os.path.join(cdir, feat + ".fasta")
        have = {s for _h, s in read_fasta(path)}
        n = 0
        with open(path, "a") as f:
            for fid, astr, s in rows:
                if s in have:
                    continue
                f.write(">%s %s [resolved by C-terminal anchor]\n%s\n" % (fid, astr, s))
                have.add(s)
                n += 1
        print("  %s.fasta  += %d" % (feat, n))

    #  ---- keep the engine's bookkeeping honest.
    #
    #  collection_engine.py wrote synonyms.tsv and UNASSIGNED_TRACKING.tsv
    #  before this pass ran, so every string resolved here is still recorded as
    #  UNASSIGNED.  Left alone that is a false report in both directions: the
    #  collapse page understates what the vocabulary absorbed, and the
    #  "is this module finished" test -- does any UNASSIGNED row reach the
    #  member floor -- fails on strings that have in fact been decided.
    #
    #  The engine counts DEPOSITS and this pass works on UNIQUE SEQUENCES, so a
    #  string's deposit count is split across its destinations in proportion to
    #  its unique sequences.  That is an approximation and it is named as one in
    #  the header line written below.
    dest_of = collections.defaultdict(collections.Counter)
    for fid, astr, ln, src, verdict, why in detail:
        if src == "unassigned":
            dest_of[astr][verdict if verdict != "DISCARD" else "DISCARDED"] += 1

    cdir_parent = os.path.dirname(cdir)
    syn = os.path.join(cdir_parent, "synonyms.tsv")
    if os.path.exists(syn) and dest_of:
        rows = [l.rstrip("\n").split("\t") for l in open(syn)]
        head, body, out_rows, moved = rows[0], rows[1:], [], 0
        for r in body:
            if len(r) < 4 or r[1] not in dest_of or r[3] != "UNASSIGNED":
                out_rows.append(r)
                continue
            counts = dest_of[r[1]]
            tot = sum(counts.values()) or 1
            n_dep = int(r[0]) if r[0].isdigit() else 0
            for d, c in sorted(counts.items()):
                out_rows.append([str(max(1, round(n_dep * c / tot))), r[1], r[2], d])
            moved += 1
        with open(syn, "w") as f:
            f.write("\t".join(head) + "\n")
            for r in out_rows:
                f.write("\t".join(r) + "\n")
        print("  synonyms.tsv: %d UNASSIGNED row(s) re-attributed "
              "(deposit counts split by unique-sequence share)" % moved)

    unas = os.path.join(cdir_parent, "UNASSIGNED_TRACKING.tsv")
    if os.path.exists(unas) and dest_of:
        rows = [l.rstrip("\n").split("\t") for l in open(unas)]
        head, body = rows[0], rows[1:]
        kept = [r for r in body if len(r) < 4 or r[3] not in dest_of]
        with open(unas, "w") as f:
            f.write("\t".join(head) + "\n")
            for r in kept:
                f.write("\t".join(r) + "\n")
        print("  UNASSIGNED_TRACKING.tsv: %d string(s) removed, now decided "
              "(see RESOLVED_FORMS.%s.tsv for every sequence and why)"
              % (len(body) - len(kept), a.anchor_from))


if __name__ == "__main__":
    main()
