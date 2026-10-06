#!/usr/bin/env python3
"""Project mature-peptide boundaries onto polyprotein-only genomes.

Most picornavirus genomes in BV-BRC carry one "polyprotein" CDS and no
mat_peptides, so a mature-peptide collection is built from the minority of
genomes somebody annotated -- 2 to 5 sequences per feature in Avihepatovirus,
against 196 polyproteins. This recovers the rest.

References are genomes whose polyprotein contains each of its own collection
mat_peptides as an exact substring: their cut coordinates are a fact of the
source annotation, not an inference. Every other polyprotein is globally
aligned to its nearest reference (blastp bitscore), and each reference cut is
carried across the alignment. Nothing is predicted; a cut is only ever moved
along an alignment.

Accepted only when
  * the target is >= --min-id identical to its reference over the alignment
  * the extracted segment is inside the feature's JSON min_len..max_len
  * the segment is >= --min-seg-id identical to the reference segment
so a cut drifting across a divergent region is refused, not shipped.

--split KEY:NEWKEY:REGEX cuts KEY in two at the first REGEX match (the NEWKEY
product starts at the match). Used where the source annotation lumps two
products that a documented site separates, e.g. DHAV 2A2|2A3.

Writes <workdir>/Extracted/<G>.<KEY>.faa, which the collection merge picks up.
"""
import argparse, collections, json, os, re, subprocess, sys, tempfile
from Bio import Align
from Bio.Align import substitution_matrices

ap = argparse.ArgumentParser()
ap.add_argument("--workdir", required=True)
ap.add_argument("--genus", required=True)
ap.add_argument("--dump", default="G")
ap.add_argument("--min-id", type=float, default=0.80,
                help="minimum polyprotein identity to transfer a cut across. "
                     "NEVER below 0.80: at 0.50 the alignment is no longer "
                     "trustworthy residue-by-residue and the cut lands wherever "
                     "the gaps fell.")
ap.add_argument("--min-seg-id", type=float, default=0.80,
                help="minimum identity between the projected segment and the "
                     "reference segment. NEVER below 0.80.")
ap.add_argument("--min-site-cons", type=float, default=0.90,
                help="a feature is only projectable if its P1 cleavage residue is "
                     "this conserved across the reference genomes.")
ap.add_argument("--p1", default="QE",
                help="residues allowed at P1 (the residue BEFORE the cut). "
                     "Picornavirus 3C cleaves Q|G almost universally, sometimes "
                     "E|G or Q|S; the S1 subsite that binds P1 is the most "
                     "conserved part of the protease across the family "
                     "(Meng 2022, J Virol 96:e0073622, PMID 35727031).")
ap.add_argument("--split", action="append", default=[])
ap.add_argument("--between", action="append", default=[],
                help="KEY:PREV:NEXT -- derive KEY on a reference as the span from PREV's end "
                     "to NEXT's start, for an internal product the source annotation omits "
                     "while annotating both its neighbours. Only fires where KEY is absent "
                     "and the gap is non-empty; the cleavage-site test still applies.")
ap.add_argument("--first-product", default="",
                help="KEY:NEXTKEY -- KEY is the polyprotein N-terminus up to the start of "
                     "NEXTKEY. For a leader, which few genomes annotate but every polyprotein "
                     "carries. Derived, not projected: the coordinates are the reference's own.")
ap.add_argument("--genus-map", default="",
                help="TSV <genome_id>\\t<genus> (e.g. column 1 and 4 of Picorna.full.tsv). "
                     "For a GROUPED module, which holds several genera: references, the "
                     "cleavage-site table and the projection all run inside one genus at a "
                     "time. Pooling them would average a genus-specific cut convention away "
                     "-- two genera that cleave the same junction one residue apart look "
                     "like one unconserved site, and a junction every genus agrees on looks "
                     "trustworthy even where the only reference is in a different genus.")
ap.add_argument("--threads", type=int, default=4)
a = ap.parse_args()
W, G = a.workdir, a.genus
C = os.path.join(W, "collections", G)

def fasta(p):
    out, h = collections.OrderedDict(), None
    for l in open(p):
        l = l.rstrip("\n")
        if l.startswith(">"): h = l[1:].split()[0]; out[h] = []
        elif h: out[h].append(l.strip())
    return {k: "".join(v) for k, v in out.items()}

mod = json.load(open(os.path.join(W, "%s_Viral_PSSM.json" % G)))[G]["features"]
skip = {"POLY", "POOL", "PRECURSORS"}
coll = {}
for f in os.listdir(C):
    k = f[:-6]
    #  only keys the module declares -- a stale collection (e.g. L before the Ldr rename) is ignored
    if f.endswith(".fasta") and k not in skip and k in mod: coll[k] = set(fasta(os.path.join(C, f)).values())
poly = fasta(os.path.join(C, "POLY.fasta"))
splits = []
for s in a.split:
    k, nk, rx = s.split(":", 2); splits.append((k, nk, re.compile(rx)))

# genome -> sequences, from the dump
md5seq = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(os.path.join(W, a.dump + ".uniq.seq"))}
bygen = collections.defaultdict(set)
for l in open(os.path.join(W, a.dump + ".id_md5")):
    fid, m = l.rstrip("\n").split("\t")
    gid = fid.split("|")[1].rsplit(".", 2)[0]
    if m in md5seq: bygen[gid].add(md5seq[m])

polyset = set(poly.values())

#  The unit that projection runs INSIDE. For a genus module that is the whole
#  module. For a grouped module it must be one genus at a time: pooling the
#  cleavage-site table across genera averages a genus-specific convention away.
#  Two genera that annotate the same junction one residue apart look like a
#  single unconserved site and the feature is dropped for everyone; and a
#  junction every genus happens to agree on looks trustworthy even where the
#  only reference sits in a different genus the target is nowhere near 80%
#  identical to. Kobuvirus VP3 is exactly the first pattern inside one genus
#  (P1 P 53%, P1' split Q61/H55), which is what makes the risk concrete.
if a.genus_map:
    g2gen = {}
    for l in open(a.genus_map):
        f = l.rstrip("\n").split("\t")
        if len(f) >= 2 and f[1]: g2gen[f[0]] = f[1]
    units, unplaced = collections.defaultdict(dict), 0
    for gid, seqs in bygen.items():
        gen = g2gen.get(gid)
        if gen is None: unplaced += 1; continue
        units[gen][gid] = seqs
    units = dict(sorted(units.items()))
    print("  genus map: %d genera, %d genomes with no genus (skipped)" % (len(units), unplaced))
else:
    units = {G: bygen}

al = Align.PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                           open_gap_score=-10, extend_gap_score=-0.5)
ALPHA = set(al.substitution_matrix.alphabet)
def clean(x): return "".join(ch if ch in ALPHA else "X" for ch in x.upper())

def mapping(r, t):
    aln = al.align(clean(r), clean(t))[0]
    m = [None] * (len(r) + 1)
    same = 0
    for (rs, re_), (ts, te) in zip(*aln.aligned):
        for d in range(re_ - rs):
            m[rs + d] = ts + d
            same += r[rs + d] == t[ts + d]
    nxt = len(t)                      # unaligned ref residues take the next aligned target position
    for i in range(len(r), -1, -1):
        if m[i] is None: m[i] = nxt
        else: nxt = m[i]
    return m, same / float(max(len(r), len(t)))

def pid(x, y):
    a2 = al.align(clean(x), clean(y))[0]
    same = sum(x[i] == y[j] for (rs, re_), (ts, te) in zip(*a2.aligned)
               for i, j in zip(range(rs, re_), range(ts, te)))
    return same / float(max(len(x), len(y)))

OUT = collections.defaultdict(dict)       # key -> {seq: header}
WHY = collections.Counter()
SPLITKEYS = {sk for sk, _, _ in splits}

for UNIT, unit_bygen in units.items():
    #  ---- references: genomes that place their own mat_peptides exactly ------
    refs = {}                             # poly seq -> {key: (start, end)}
    for gid, seqs in unit_bygen.items():
        for p in seqs & polyset:
            cuts = {}
            for k, cs in coll.items():
                for s in seqs & cs:
                    i = p.find(s)
                    if i >= 0 and p.find(s, i + 1) < 0: cuts[k] = (i, i + len(s))
            if a.first_product:
                k0, k1 = a.first_product.split(":")
                if k1 in cuts and cuts[k1][0] > 0: cuts[k0] = (0, cuts[k1][0])
            #  An INTERNAL product the source never annotates, recovered as the
            #  span its two annotated neighbours leave between them. Like
            #  --first-product this is a derivation from the reference's own
            #  coordinates, not a prediction -- but it is only sound where the
            #  neighbours are genuinely adjacent to it, so it fires only when the
            #  key is absent and the gap is non-empty, and the cleavage-site
            #  table below still has to accept the resulting cut.
            #  Hunnivirus is the case it was written for: its two reference
            #  genomes annotate every product except VP1 and tile the whole
            #  polyprotein, leaving VP3.end -> 2A.start = 234 aa of exactly one
            #  unannotated product, which the genome organisation
            #  (L-VP4-VP2-VP3-VP1-2A) identifies as VP1. 122 Hunnivirus
            #  polyproteins had VP1 called on none of them.
            for spec in a.between:
                k0, kp, kn = spec.split(":")
                if k0 in cuts or kp not in cuts or kn not in cuts: continue
                lo_, hi_ = cuts[kp][1], cuts[kn][0]
                if lo_ < hi_: cuts[k0] = (lo_, hi_)
            if len(cuts) >= 5 and (p not in refs or len(cuts) > len(refs[p])): refs[p] = cuts
    unit_poly = {p for seqs in unit_bygen.values() for p in seqs & polyset}
    print("\n  == %s: %d polyproteins, %d references with >=5 placed mat_peptides"
          % (UNIT, len(unit_poly), len(refs)))
    if not refs:
        if unit_poly: WHY["%s: no reference genome places >=5 mat_peptides" % UNIT] += len(unit_poly)
        continue

    #  ---- cleavage-site conservation, measured on THIS unit's references -----
    #  A cut is only trustworthy if the protease site it sits on is conserved.
    #  For each feature take P1 (the residue immediately BEFORE the cut) and P1'
    #  (the first residue of the product) across every reference that places it.
    site = collections.defaultdict(lambda: (collections.Counter(), collections.Counter()))
    for p, cuts in refs.items():
        for k, (st, en) in cuts.items():
            site[k][0][p[st - 1] if st > 0 else "-"] += 1
            site[k][1][p[st]] += 1
    ok_site, site_req, site_p1p = {}, {}, {}
    print("  cleavage-site conservation on the references (P1 | P1'):")
    for k in sorted(site):
        p1, p1p = site[k]
        n = sum(p1.values()); top, cnt = p1.most_common(1)[0]
        frac = cnt / float(n)
        #  CONSERVATION is the criterion, not the identity of the residue. Q|G is
        #  the canonical 3C site and covers most junctions, but a picornavirus
        #  polyprotein also carries cuts 3C does not make -- the VP0 -> VP4+VP2
        #  maturation cleavage is autocatalytic and sits on a conserved K|A, and
        #  the StopGo site is G|P. A site conserved at 99% over 158 references is
        #  trustworthy whatever residue it uses; one at 54% is not, whatever
        #  residue it uses.
        keep = frac >= a.min_site_cons or top == "-"
        ok_site[k] = keep
        if keep and top != "-": site_req[k] = top
        kind = ("polyprotein N-terminus" if top == "-"
                else "3C-type" if top in a.p1 else "conserved, non-3C")
        shown = "%s %.0f%% (n=%d)" % (top, 100 * frac, n)
        #  3C accepts EITHER Q or E at P1 -- Q|G is canonical and E|G is the
        #  documented minor form (Meng 2022, J Virol 96:e0073622, PMID 35727031),
        #  and the S1 pocket that reads P1 does not distinguish them. So a
        #  junction that is Q in some references and E in others is one conserved
        #  3C site, not an unconserved one, and demanding a single exact residue
        #  is stricter than the biology: Hunnivirus VP1 is PYLLLE|GEDSQV and
        #  PYLALQ|GEDTQI -- the same site, scored 50% by exact residue.
        #  This only ever applies as a FALLBACK, and only with corroboration:
        #  P1' must itself be conserved, so the second half of the motif has to
        #  agree before the first half is allowed to vary. Anything already
        #  passing on an exact residue is untouched.
        if not keep:
            qe = sum(p1[r] for r in a.p1) / float(n)
            ptop, pcnt = p1p.most_common(1)[0]
            if qe >= a.min_site_cons and pcnt / float(n) >= a.min_site_cons:
                keep = True
                ok_site[k] = True
                site_req[k] = set(a.p1)
                site_p1p[k] = ptop          # the corroborating half of the motif
                kind = "3C-type (%s, P1' %s)" % ("/".join(sorted(a.p1)), ptop)
                shown = "%s %.0f%% (n=%d)" % ("/".join(sorted(a.p1)), 100 * qe, n)
        print("    %-8s P1 %-20s P1' %-20s %-26s %s" % (
            k, shown,
            " ".join("%s%d" % kv for kv in p1p.most_common(3)),
            kind if keep else "",
            "ok" if keep else "NOT CONSERVED (%.0f%%) -- feature not projected" % (100 * frac)))

    #  ---- nearest reference per target, by blastp bitscore -------------------
    td = tempfile.mkdtemp()
    ref_ids = {"r%d" % i: p for i, p in enumerate(refs)}
    tgt_ids = {"t%d" % i: p for i, p in enumerate(sorted(unit_poly))}
    open(td + "/r.faa", "w").write("".join(">%s\n%s\n" % kv for kv in ref_ids.items()))
    open(td + "/t.faa", "w").write("".join(">%s\n%s\n" % kv for kv in tgt_ids.items()))
    subprocess.run(["makeblastdb", "-in", td + "/r.faa", "-dbtype", "prot", "-out", td + "/r"],
                   check=True, capture_output=True)
    bl = subprocess.run(["blastp", "-query", td + "/t.faa", "-db", td + "/r",
                         "-outfmt", "6 qseqid sseqid bitscore", "-max_target_seqs", "2000",
                         "-evalue", "1e-20", "-num_threads", str(a.threads)],
                        check=True, capture_output=True, text=True).stdout
    hits = collections.defaultdict(dict)
    for l in bl.splitlines():
        q, s, b = l.split("\t"); b = float(b)
        hits[q][s] = max(b, hits[q].get(s, 0))

    #  ---- project ------------------------------------------------------------
    allkeys = sorted({k for c in refs.values() for k in c})
    for tid, t in tgt_ids.items():
        if tid not in hits: WHY["no reference hit"] += 1; continue
        ranked = sorted(hits[tid], key=lambda x: -hits[tid][x])
        cache = {}
        for k in allkeys:
            rid = next((x for x in ranked if k in refs[ref_ids[x]]), None)
            if rid is None: WHY[k + " no reference places it"] += 1; continue
            r = ref_ids[rid]
            if rid not in cache:
                cache[rid] = (list(range(len(r) + 1)), 1.0) if t == r else mapping(r, t)
            m, ident = cache[rid]
            if ident < a.min_id: WHY[k + " polyprotein below %.0f%% id" % (100 * a.min_id)] += 1; continue
            s, e = refs[r][k]
            seg, rseg = t[m[s]:m[e]], r[s:e]
            lo, hi = mod[k]["min_len"] if k in mod else 0, mod[k]["max_len"] if k in mod else 10 ** 6
            #  A key about to be SPLIT still carries the LUMPED length here; its
            #  JSON window describes the post-split product. Checking a pre-split
            #  segment against a post-split window rejects every one of them
            #  (Avihepatovirus 2A2: 183 refused on length alone).
            if k in SPLITKEYS: lo, hi = 0, 10 ** 6
            #  Inside a GROUPED module the JSON window describes every genus
            #  lumped together, and a genus whose product is genuinely a
            #  different size fails it on length alone: PICO_NOLDR_VP4 declares
            #  VP4 54-82 while Tremovirus VP4 is ~20. Measure against the
            #  REFERENCE segment instead -- it is in the same genus as the
            #  target. The >= --min-seg-id identity test already does the real
            #  work here; this is only a sanity bound. (Same failure mode as the
            #  Betaflexiviridae CP window, which was set from the two biggest
            #  genera and silently deleted 90 real Citrivirus coat proteins.)
            if a.genus_map: lo, hi = int(0.75 * len(rseg)), int(1.34 * len(rseg)) + 1
            if not seg or not (lo <= len(seg) <= hi) or "X" in seg or "*" in seg:
                WHY[k + " length/ambiguity"] += 1; continue
            if seg != rseg and pid(rseg, seg) < a.min_seg_id: WHY[k + " segment id"] += 1; continue
            if not ok_site.get(k, False): WHY[k + " cleavage site not conserved on refs"] += 1; continue
            #  the projected cut must land on the P1 the references use
            need = site_req.get(k)
            if need and m[s] > 0:
                got = t[m[s] - 1]
                if (got not in need) if isinstance(need, set) else (got != need):
                    WHY[k + " P1 is not %s at the projected cut"
                        % ("/".join(sorted(need)) if isinstance(need, set) else need)] += 1
                    continue
            #  where P1 was allowed to vary across the Q/E class, P1' must match
            needp = site_p1p.get(k)
            if needp and t[m[s]] != needp:
                WHY[k + " P1' is not %s at the projected cut" % needp] += 1; continue
            OUT[k][seg] = "%s|%s|%s|ref=%s|id=%.2f" % (G, UNIT, k, rid, ident)

out, why = OUT, WHY
#  splits: applied to the projected segments and to the existing collection alike
for k, nk, rx in splits:
    pool = dict(out.pop(k, {}))
    for s in coll.get(k, ()): pool.setdefault(s, "%s|%s|collection" % (G, k))
    for s, h in pool.items():
        mt = rx.search(s)
        if not mt or mt.start() == 0: why[k + " split motif absent"] += 1; continue
        out[k][s[:mt.start()]] = h + "|split"
        out[nk][s[mt.start():]] = h.replace("|%s|" % k, "|%s|" % nk) + "|split"

os.makedirs(os.path.join(W, "Extracted"), exist_ok=True)
print("\n  %-8s %8s %8s %8s" % ("key", "extract", "new", "in coll"))
for k in sorted(out):
    p = os.path.join(W, "Extracted", "%s.%s.faa" % (G, k))
    with open(p, "w") as f:
        for i, (s, h) in enumerate(sorted(out[k].items())):
            f.write(">ext%d|%s\n%s\n" % (i + 1, h, s))
    new = len(set(out[k]) - coll.get(k, set()))
    print("  %-8s %8d %8d %8d" % (k, len(out[k]), new, len(coll.get(k, ()))))
print("\n  refused: %s" % {k: v for k, v in sorted(why.items(), key=lambda kv: -kv[1]) if v})
