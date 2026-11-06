#!/usr/bin/env python3
"""Ask of every collection: is this ONE protein, and is it ITS OWN protein?

    python3 audit_bins.py --workdir . --module Mastadenovirus \
        --write-md BIN_AUDIT.Mastadenovirus.md
    python3 audit_bins.py ... --accept accept_bins_Mastadenovirus.tsv --write

`audit_collection_lengths.py` asks whether the MEMBERS of a collection are the
right size, and `audit_contamination.py` moves individual mislabelled
sequences between collections. Neither asks whether the collection should
exist at all, and that is a different defect with the same invisible
signature.

WHAT GOES WRONG. The accessory-discovery pass (`filter_unchar.py`) declares
every unassigned cluster that reaches the member floor as UNC1..UNCn. It has
no way to know that a cluster is a pile of partials, or a second copy of a
protein that already has a name, or three unrelated things that happened to
cluster together. Each of those becomes a feature, gets a PSSM, and then calls
real loci that belong to somebody else.

Measured on Adenoviridae, 5 Nov 2026 -- Mastadenovirus, 115 collections:

  3 bins were a protein that already had a bin.   UNC16 was POL at 99.8%
    identity, 1187 aa against POL's 1190.
  15 bins were partials.                          UNC23 was 100% identical to
    E1B55K over 95% of itself and 17% of E1B55K. E3_CR1D -- a NAMED feature --
    was 29 members of a 46 aa fragment of E3_CR1B.
  7 bins matched a named feature at 29-37% across their full length: too
    divergent to merge, too close to profile apart.

Those 25 bins held 293 sequences and caused most of the module's 458
co-located duplicate calls. UNC34, seven sequences that mmseqs splits into
three groups sharing 23-32%, scored 80 bits on E1B19K's locus against 33 on
its own and produced 90 false E1B19K calls on its own.

WHY IT IS NOT A CUTOFF PROBLEM, which is the trap. A cutoff separates a
profile's true loci from its false ones. These profiles score HIGHER on the
foreign locus than on their own, because the foreign protein is what they were
built from. There is no threshold between the two, and the only symptom is
duplicate calls -- which look exactly like a cutoff problem and get chased
through coverage runs, the most expensive place in the pipeline to find this
and the least likely place to understand it.

THE VERDICTS, in the order they are tried. `tgt` is the feature that most of
the bin's members hit best, blastp, excluding the bin itself.

  DUPLICATE   >=95% identity AND (same median length as tgt, or full coverage
              both ways). The protein already has a bin. Retire this one.
  FRAGMENT    >=90% identity, >=85% of the QUERY aligned, <=75% of the
              SUBJECT. The members are complete as deposited and partial as
              biology -- a slice of tgt. Also caught at >=84% when the hit is
              short on both sides (<=35% of tgt).  The curator's standing rule
              is that partials are discarded, and a partial that became its
              own feature is still a partial.
  INCOHERENT  the members resemble the foreign feature about as much as they
              resemble EACH OTHER -- within-bin identity is less than
              --coherence points above the identity to tgt. The bin is not one
              protein, so its profile is a fold detector. Discard.

Everything else keeps its bin, and a bin that matches a named feature at
25-55% across its full length WHILE being internally coherent is reported as
a PARALOGUE: information, not a defect.

THE TRAP IN THE INCOHERENT TEST, which cost a rule rewrite on 5 Nov 2026.
The first version asked only "does this bin match a named feature at 25-55%
across its full length". That test is SYMMETRIC, so it flags both halves of
every divergent paralogue pair -- Aviadenovirus ORF12/UNC8 at 44-46%,
ORF20/UNC47, ORF9/UNC4, Barthadenovirus RH4 and RH6 against RH2 -- and it
proposed retiring ORF20, a named feature with 99 members. Divergent
paralogues are SUPPOSED to sit at 30-45%; that is what makes them paralogues.
What separates UNC34 from ORF20 is not its distance to E1B19K but its
distance to itself: UNC34's seven members fall into three groups sharing
23-32%, which is no better than their 34% to E1B19K, while ORF20's 99 members
agree with each other far better than with anything else. So the bin is
measured against itself first, and the foreign hit only decides what it is a
fold detector FOR.

PAIRS THAT MUST NOT BE JUDGED. Nested and spliced forms are SUPPOSED to be
near-identical -- K22/K33, UXP/U_EXON, E1A/E1A_243R, E4ORF67/E4ORF6, and the
rotavirus/parvovirus nested capsid forms all trip FRAGMENT by construction.
Declare them with --protect (or in protected_pairs.json beside the
collections) and they are skipped. Getting this list wrong deletes real
biology, so the default is to skip nothing and make the curator name them.

A verdict is a FINDING, not an action. Resolve each one -- re-bin it, merge
it, or sign it off in --accept with a reason -- and only then does --write
retire what is left. build_features.py refuses to build profiles until this
report exists and reads **Unresolved: 0**.
"""
import argparse, collections, glob, json, os, shutil, statistics
import subprocess, sys, tempfile

BV = "/Applications/BV-BRC.app/runtime/bin"
if os.path.isdir(BV):
    os.environ["PATH"] = BV + ":" + os.environ.get("PATH", "")


def read_fasta(path):
    seqs, name = collections.OrderedDict(), None
    for line in open(path, errors="replace"):
        if line.startswith(">"):
            name = line[1:].rstrip("\n")
            seqs.setdefault(name, [])
        elif name is not None:
            seqs[name].append(line.strip())
    return collections.OrderedDict((k, "".join(v)) for k, v in seqs.items())


def collections_dir(workdir, module):
    for d in (os.path.join(workdir, module, "collections"),
              os.path.join(workdir, "collections", module)):
        if os.path.isdir(d):
            return d
    sys.exit("no collections directory for %s under %s" % (module, workdir))


def load_protected(cdir, extra):
    """Pairs that are near-identical BY DESIGN and must never be judged."""
    pairs = set()
    p = os.path.join(cdir, "protected_pairs.json")
    if os.path.exists(p):
        for a, b in json.load(open(p)):
            pairs.add((a, b)); pairs.add((b, a))
    for item in extra:
        a, _, b = item.partition(":")
        if b:
            pairs.add((a, b)); pairs.add((b, a))
    return pairs


def nesting(cdir, feats):
    """Which bins are an exact substring of another bin, IN THE SAME GENOME.

    This is the test the capsid audit already applies to VP1/VP2, and it
    answers in one pass what a literature search answers slowly and often not
    at all. A nested form is an alternative START: it carries the parent's
    C-terminus and is missing the parent's N-terminus. A deposit that begins
    at the parent's residue 1 and stops short is the opposite -- a truncated
    parent under another name -- and no identity or length rule separates the
    two, because both are merely "shorter".

    Aviadenovirus UNC81, 5 Nov 2026: three deposits named `Lip-3`, 187 aa,
    each an exact C-terminal substring of ITS OWN genome's 560 aa `ORF19`.
    Fowl aviadenovirus C really does carry a second right-end lipase gene,
    ORF19A, so the question was fair -- but ORF19A is a distinct homolog with
    its own sequence, not a verbatim copy of ORF19's tail, and no paper uses
    `Lip-3`. Same genome plus exact substring means ONE LOCUS annotated
    twice, whatever the labels say.

    Returns {bin: (parent, n_nested, n_tested, end)} where end is "C" for a
    shared C-terminus (alternative start), "N" for a shared N-terminus
    (truncated parent), or "mid".
    """
    by_feat = {}
    for f in feats:
        d = {}
        for h, sq in read_fasta(os.path.join(cdir, f + ".fasta")).items():
            d.setdefault(h.split()[0], []).append(sq)
        by_feat[f] = d
    out = {}
    for f in feats:
        per = {}
        for g in feats:
            if g == f:
                continue
            nested = tested = 0
            ends = collections.Counter()
            example = None
            for gid, seqs in by_feat[f].items():
                others = by_feat[g].get(gid)
                if not others:
                    continue
                tested += 1
                for sq in seqs:
                    hit = next((o for o in others
                                if len(o) > len(sq) and sq in o), None)
                    if hit:
                        nested += 1
                        ends["C" if hit.endswith(sq) else
                             "N" if hit.startswith(sq) else "mid"] += 1
                        if example is None:
                            example = (gid, len(sq), len(hit))
                        break
            if tested:
                per[g] = (nested, tested,
                          ends.most_common(1)[0][0] if ends else "mid",
                          example)
        if per:
            out[f] = per
    return out


def measure(cdir, feats, link_id=40.0, link_cov=0.5):
    """Best FOREIGN feature for every member of every collection."""
    tmp = tempfile.mkdtemp(prefix="binaudit_")
    try:
        allf = os.path.join(tmp, "all.faa")
        with open(allf, "w") as o:
            for f in feats:
                for i, (_, s) in enumerate(read_fasta(
                        os.path.join(cdir, f + ".fasta")).items()):
                    o.write(">%s|%d\n%s\n" % (f, i, s))
        subprocess.run(["makeblastdb", "-in", allf, "-dbtype", "prot",
                        "-out", os.path.join(tmp, "db")],
                       check=True, capture_output=True)
        out = {}
        for f in feats:
            fa = os.path.join(cdir, f + ".fasta")
            seqs = read_fasta(fa)
            if not seqs:
                continue
            lens = sorted(len(s) for s in seqs.values())
            r = subprocess.run(
                ["blastp", "-query", fa, "-db", os.path.join(tmp, "db"),
                 "-outfmt", "6 qseqid sseqid pident length qlen slen bitscore",
                 "-max_target_seqs", "5000", "-evalue", "1e-3",
                 "-num_threads", "4"], capture_output=True, text=True)
            best = {}
            for line in r.stdout.splitlines():
                q, s, pid, aln, qlen, slen, bits = line.split("\t")
                feat = s.split("|")[0]
                if feat == f:
                    continue
                bits = float(bits)
                if q not in best or bits > best[q][0]:
                    best[q] = (bits, feat, float(pid), int(aln),
                               int(qlen), int(slen))
            #  WITHIN-BIN IDENTITY, the thing the first version of this script
            #  did not measure. For each member, the best identity to any
            #  OTHER member of the same bin; the bin's score is the median of
            #  those. A member with no within-bin hit at all scores 0, which
            #  is the correct answer -- it belongs to nothing here.
            selfid, biggest, ncomp = 0.0, 1.0, 1
            if len(seqs) > 1:
                sr = subprocess.run(
                    ["blastp", "-query", fa, "-subject", fa,
                     "-outfmt", "6 qseqid sseqid pident length qlen",
                     "-evalue", "10", "-num_threads", "4"],
                    capture_output=True, text=True)
                bestself = {}
                order = list(seqs)
                for line in sr.stdout.splitlines():
                    q, s, pid, aln, qlen = line.split("\t")
                    #  ids repeat inside a collection (one genome depositing
                    #  two copies), so a name match is not proof of self --
                    #  but a full-length 100% hit is, and dropping it is the
                    #  conservative error here.
                    if q == s and float(pid) > 99.5 and int(aln) >= int(qlen):
                        continue
                    pid = float(pid)
                    if pid > bestself.get(q, 0):
                        bestself[q] = pid
                vals = [bestself.get(k.split()[0], 0.0) for k in order]
                selfid = statistics.median(vals) if vals else 0.0
                #  CONNECTIVITY, because the median cannot see a SPLIT bin.
                #  UNC34 is three groups -- one trio at 95-99%, a second at
                #  39-48%, one sequence alone -- sharing 23-32% ACROSS the
                #  groups. Its median within-bin identity is 48.5%, which
                #  clears any margin over its 34% to E1B19K, so the median
                #  calls it coherent and it is not. Single-linkage at a
                #  modest threshold answers the real question: is this bin
                #  ONE thing?
                parent = {k: k for k in order}
                def find(x):
                    while parent[x] != x:
                        parent[x] = parent[parent[x]]; x = parent[x]
                    return x
                idx = {}
                for k in order:
                    idx.setdefault(k.split()[0], k)
                for line in sr.stdout.splitlines():
                    q, s, pid, aln, qlen = line.split("\t")
                    if q == s or float(pid) < link_id:
                        continue
                    if int(aln) < link_cov * int(qlen):
                        continue
                    qq, ss = idx.get(q), idx.get(s)
                    if qq and ss:
                        a, b = find(qq), find(ss)
                        if a != b:
                            parent[a] = b
                comp = collections.Counter(find(k) for k in order)
                biggest = comp.most_common(1)[0][1] / float(len(order))
                ncomp = len(comp)
            rec = dict(n=len(seqs), medlen=int(statistics.median(lens)),
                       tgt="-", pid=0.0, qcov=0.0, scov=0.0, frac=0.0,
                       selfid=selfid, biggest=biggest, ncomp=ncomp)
            if best:
                votes = collections.Counter(v[1] for v in best.values())
                tgt, nv = votes.most_common(1)[0]
                hits = [v for v in best.values() if v[1] == tgt]
                rec.update(tgt=tgt, frac=nv / float(len(seqs)),
                           pid=statistics.median(h[2] for h in hits),
                           qcov=statistics.median(h[3] / float(h[4]) for h in hits),
                           scov=statistics.median(h[3] / float(h[5]) for h in hits))
            out[f] = rec
        return out
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def verdict(f, rec, medlens, protected, coherence, min_component):
    #  SAME-GENOME EXACT NESTING FIRST. It is structural rather than
    #  statistical: one genome deposited the same locus twice, so no identity
    #  threshold has to be argued about. It outranks the blast verdicts
    #  because it is the only one that cannot be a coincidence of similarity.
    #  BUT IT MUST BE THE RULE IN THAT BIN, NOT ONE DEPOSIT. A form that is
    #  really an alternative start of its parent nests in essentially every
    #  genome depositing both -- UNC81 was 2 of 2, UNC17 is 1 of 1. A single
    #  nesting out of many co-deposits is the opposite finding: one bad
    #  deposit, and usually in the PARENT.
    #
    #  Mastadenovirus PROTEASE, 6 Nov 2026, caught in pre-flight. Flagged as
    #  an exact substring of HEXON in 1 of 138 shared genomes, which would
    #  have retired the 252-member adenain bin (median 205 aa, and the
    #  literature's 23K adenain is 204-206). The nesting is real and the
    #  conclusion is backwards: genome 3115880.8 deposited a 1122 aa protein
    #  called "Hexain" with the 204 aa adenain verbatim inside it, and its
    #  sibling 3115874.8 deposited another at 1124 aa, against a HEXON median
    #  of 942 over 570 members. Two fused gene models in HEXON accused
    #  PROTEASE. Barthadenovirus RH1 was the same shape at 1 of 4.
    per = rec.get("nest") or {}
    best = max(((g, v) for g, v in per.items()
                if v[0] and (v[0] >= 2 or v[1] == 1)
                and v[0] / float(v[1]) >= 0.5),
               key=lambda kv: kv[1][0], default=None)
    if best and (f, best[0]) not in protected:
        parent = best[0]
        n, tested, end = best[1][:3]
        shape = {"C": "shares %s's C-terminus and is missing its N-terminus "
                      "-- shaped like an alternative start, so it needs a "
                      "literature name before it can be its own feature",
                 "N": "begins at %s residue 1 and stops short -- a truncated "
                      "%s under another name",
                 "mid": "sits inside %s"}[end]
        return ("NESTED",
                "exact substring of %s in %d of %d shared genome(s); %s"
                % (parent, n, tested,
                   shape % ((parent, parent) if end == "N" else parent)))

    #  A minority nesting is still a defect -- it is just not THIS bin's.
    #  Name the deposit so it reaches the length and contamination audits
    #  instead of disappearing with the retracted verdict.
    odd = max(((g, v) for g, v in per.items()
               if v[0] and v[1] > 1 and v[0] / float(v[1]) < 0.5),
              key=lambda kv: kv[1][0], default=None)
    if odd and (f, odd[0]) not in protected:
        n, tested, _end, ex = odd[1]
        return ("ODD-DEPOSIT",
                "nests inside %s in only %d of %d shared genomes -- not an "
                "alternative start. Genome %s deposited %s at %d aa inside a "
                "%s of %d aa; check that deposit, not this bin"
                % (odd[0], n, tested, ex[0] if ex else "?", f,
                   ex[1] if ex else 0, odd[0], ex[2] if ex else 0))
    tgt = rec["tgt"]
    if tgt == "-" or rec["frac"] < 0.5 or (f, tgt) in protected:
        #  A bin can be several proteins without resembling any NAMED one.
        #  Discarding it would throw away novel proteins, so it is reported
        #  as a split to be made, not as something to retire.
        if rec["n"] >= 4 and rec["biggest"] < min_component:
            return ("SPLIT",
                    "%d groups, largest holds %.0f%% of the bin -- this is "
                    "more than one protein and resembles no named feature"
                    % (rec["ncomp"], rec["biggest"] * 100))
        return None
    pid, qc, sc = rec["pid"], rec["qcov"], rec["scov"]
    tl = medlens.get(tgt, 0)
    same_len = bool(tl) and 0.80 <= rec["medlen"] / float(tl) <= 1.25
    if pid >= 95 and (same_len or (qc >= 0.85 and sc >= 0.85)):
        return ("DUPLICATE", "same protein as %s (%.1f%% id)" % (tgt, pid))
    #  A FRAGMENT VERDICT NEEDS STRUCTURAL SUPPORT, OR AT LEAST NO STRUCTURAL
    #  CONTRADICTION. Cross-genome identity alone cannot tell "a piece of that
    #  protein" from "a relative that shares a conserved region", and for a
    #  conserved-domain family it reliably says the wrong thing.
    #
    #  Mastadenovirus E3_CR1D, 5 Nov 2026: 29 deposits, median 46 aa, 96.6%
    #  identical to E3_CR1B over 7% of it -- a textbook FRAGMENT by identity.
    #  But NOT ONE of the 29 is a substring of a CR1B from its own genome, and
    #  nine of them had a same-genome CR1B available to be a substring of. CR1
    #  stands for conserved region 1: the CR1-alpha/beta/gamma/delta proteins
    #  share a conserved region BY DEFINITION, so the rule was measuring the
    #  very thing that makes them a family and calling it contamination.
    #  Retiring it would have deleted a real protein.
    tnest = (rec.get("nest") or {}).get(tgt)
    contradicted = bool(tnest) and tnest[1] >= 3 and tnest[0] == 0
    if not contradicted:
        if pid >= 90 and qc >= 0.85 and sc <= 0.75:
            return ("FRAGMENT", "partial of %s (%.1f%% id over %.0f%% of %s)%s"
                    % (tgt, pid, sc * 100, tgt,
                       "" if tnest and tnest[0] else
                       " -- identity only, no co-deposit to confirm it"))
        if pid >= 84 and qc < 0.85 and sc <= 0.35:
            return ("FRAGMENT", "partial of %s (%.1f%% id, %.0f%% of %s)%s"
                    % (tgt, pid, sc * 100, tgt,
                       "" if tnest and tnest[0] else
                       " -- identity only, no co-deposit to confirm it"))
    elif pid >= 84:
        return ("CONSERVED-REGION",
                "%.1f%% to %s over %.0f%% of it, but 0 of %d same-genome "
                "co-deposits nest inside %s -- a shared conserved region, "
                "not a piece of it"
                % (pid, tgt, sc * 100, tnest[1], tgt))
    if 25 <= pid <= 55 and qc >= 0.85 and sc >= 0.80 and rec["frac"] >= 0.85:
        #  Divergent from tgt. Only a defect if it is equally divergent from
        #  ITSELF -- otherwise it is tgt's paralogue, which is a fact about
        #  the taxon rather than a fault in the bin.
        if rec["biggest"] < min_component:
            return ("INCOHERENT",
                    "not one protein -- %d groups, largest holds %.0f%% of "
                    "the bin; best foreign match %s at %.0f%%"
                    % (rec["ncomp"], rec["biggest"] * 100, tgt, pid))
        if rec["selfid"] < pid + coherence:
            return ("INCOHERENT",
                    "members agree with each other no better (%.0f%%) than "
                    "with %s (%.0f%%)" % (rec["selfid"], tgt, pid))
        return ("PARALOGUE",
                "%.0f%% to %s full length, but one group, internally "
                "coherent at %.0f%%" % (pid, tgt, rec["selfid"]))
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--module", required=True)
    ap.add_argument("--protect", action="append", default=[],
                    metavar="A:B",
                    help="pair that is near-identical by design and must not "
                         "be judged (K22:K33). Repeatable; also read from "
                         "protected_pairs.json beside the collections.")
    ap.add_argument("--link-id", type=float, default=40.0,
                    help="percent identity at which two members of a bin are "
                         "linked when testing whether the bin is ONE thing "
                         "(default 40)")
    ap.add_argument("--link-cov", type=float, default=0.5,
                    help="fraction of the query a link must span (default 0.5)")
    ap.add_argument("--min-component", type=float, default=0.70,
                    help="fraction of a bin that must fall in its largest "
                         "connected component (default 0.70)")
    ap.add_argument("--coherence", type=float, default=10.0,
                    help="points by which within-bin identity must exceed the "
                         "identity to the foreign feature for the bin to "
                         "count as one protein (default 10)")
    ap.add_argument("--accept", help="TSV of signed-off bins: FEATURE<TAB>why")
    ap.add_argument("--write-md", help="write the report here")
    ap.add_argument("--write", action="store_true",
                    help="retire every unaccepted bin to collections/retired/ "
                         "-- the collection is MOVED, never deleted, so the "
                         "decision can be read back and reversed.")
    a = ap.parse_args()

    cdir = collections_dir(a.workdir, a.module)
    feats = sorted(os.path.basename(p)[:-6]
                   for p in glob.glob(os.path.join(cdir, "*.fasta"))
                   if not p.endswith(".outliers.fasta"))
    protected = load_protected(cdir, a.protect)
    accepted = {}
    if a.accept and os.path.exists(a.accept):
        for line in open(a.accept):
            if line.strip() and not line.startswith("#"):
                k, _, why = line.rstrip("\n").partition("\t")
                accepted[k.strip()] = why.strip() or "signed off"

    recs = measure(cdir, feats, a.link_id, a.link_cov)
    for f, nst in nesting(cdir, feats).items():
        if f in recs:
            recs[f]["nest"] = nst
    medlens = {f: r["medlen"] for f, r in recs.items()}

    findings = []
    for f in feats:
        v = (verdict(f, recs[f], medlens, protected, a.coherence,
                     a.min_component) if f in recs else None)
        if v:
            findings.append((f, v[0], v[1], recs[f]))

    #  A MUTUAL DUPLICATE MUST NOT RETIRE BOTH HALVES, which is how this
    #  tool would delete a protein outright. Aviadenovirus UNC11 (15
    #  members) and UNC15 (13) are the same thing at 98.9%, so each is the
    #  other's duplicate and each earns a DUPLICATE verdict. Retiring the
    #  pair removes the protein from the module. The survivor is the bin
    #  with more members -- more evidence, better profile -- and ties go to
    #  the name that sorts first, so the choice is at least reproducible.
    dup = {f: rec["tgt"] for f, cls, _, rec in findings if cls == "DUPLICATE"}
    keep = set()
    for f, tgt in dup.items():
        if dup.get(tgt) == f:
            winner = max((f, tgt), key=lambda x: (recs[x]["n"], x == min(f, tgt)))
            keep.add(winner)
    if keep:
        findings = [(f, cls, why, rec) for f, cls, why, rec in findings
                    if f not in keep]
        for f in sorted(keep):
            findings.append((f, "KEPT", "survivor of a mutual-duplicate pair "
                             "with %s -- retiring both would delete the "
                             "protein" % dup[f], recs[f]))

    #  PARALOGUE and SPLIT are reported, never blocking. A coherent
    #  divergent bin is a fact about the taxon. A bin that splits without
    #  resembling any named feature may be several NOVEL proteins -- and at
    #  any linkage threshold loose enough to hold a real family together
    #  across host species, the test stops finding anything, while at 40% it
    #  flags Barthadenovirus FIBER, K22, LH1 and LH2, all of them real. So it
    #  goes to the curator as a question, not to the build as a refusal.
    defects = [x for x in findings
               if x[1] not in ("PARALOGUE", "SPLIT", "KEPT",
                               "CONSERVED-REGION", "ODD-DEPOSIT")]
    unresolved = [x for x in defects if x[0] not in accepted]
    lines = ["# Bin audit -- %s" % a.module, "",
             "Does every collection hold ONE protein, and its OWN protein?",
             "", "**Unresolved: %d**  (defects %d, signed off %d, "
             "paralogues %d, collections %d)"
             % (len(unresolved), len(defects), len(defects) - len(unresolved),
                len(findings) - len(defects), len(feats)), ""]
    if protected:
        lines += ["Protected pairs (not judged): %s" % ", ".join(
            sorted("%s/%s" % p for p in protected if p[0] < p[1])), ""]
    for cls in ("NESTED", "DUPLICATE", "FRAGMENT", "INCOHERENT", "KEPT",
                "ODD-DEPOSIT", "CONSERVED-REGION", "SPLIT", "PARALOGUE"):
        sel = [x for x in findings if x[1] == cls]
        if not sel:
            continue
        lines += ["## %s -- %d bin(s), %d sequence(s)%s"
                  % (cls, len(sel), sum(x[3]["n"] for x in sel),
                     "  (reported, not retired)"
                     if cls in ("PARALOGUE", "SPLIT", "KEPT",
                                "CONSERVED-REGION", "ODD-DEPOSIT") else ""),
                  "",
                  "| bin | n | median aa | finding | signed off |",
                  "|---|---|---|---|---|"]
        for f, _, why, rec in sel:
            lines.append("| `%s` | %d | %d | %s | %s |"
                         % (f, rec["n"], rec["medlen"], why,
                            accepted.get(f, "")))
        lines.append("")
    if not defects:
        lines += ["Every collection is one protein and its own protein.", ""]

    report = "\n".join(lines)
    print(report)
    if a.write_md:
        open(a.write_md, "w").write(report + "\n")
        print("wrote %s" % a.write_md)

    if a.write:
        rdir = os.path.join(cdir, "retired")
        os.makedirs(rdir, exist_ok=True)
        moved = 0
        for f, cls, why, _ in unresolved:
            src = os.path.join(cdir, f + ".fasta")
            shutil.move(src, os.path.join(rdir, f + ".fasta"))
            with open(os.path.join(rdir, "RETIRED.tsv"), "a") as o:
                o.write("%s\t%s\t%s\n" % (f, cls, why))
            moved += 1
        print("retired %d bin(s) to %s -- REBUILD profiles and the JSON" %
              (moved, rdir))
    sys.exit(1 if unresolved and not a.write else 0)


if __name__ == "__main__":
    main()
