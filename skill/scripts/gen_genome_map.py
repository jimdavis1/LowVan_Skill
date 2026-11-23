#!/usr/bin/env python3
"""To-scale genome map of one annotated reference genome, as an inline SVG
<figure> for the coverage audit.

The fifth page-figure, and the one the kit had no generator for. Adenoviridae
shipped `reports/generators/adenoviridae/genome_map.py`, which was written
against that build's own directory layout and palette, so no other module ever
got a map: ten families were built and asked for one after the fact. This is
that script generalised -- paths come from the command line, segmented genomes
get one axis per segment, and the colours are the coverage-audit page's own
CSS variables rather than the adenovirus page's.

    python3 gen_genome_map.py --gto panel/145856.193.qual.gto \
        --module-json Picobirnaviridae_Viral_PSSM.json \
        --module Picobirnaviridae --out picobirnaviridae-genome-map.html

Output is a bare `<figure>` fragment, which is what `gen_coverage_audit.py`
lifts verbatim for a `{"type": "figure", "from": "..."}` block. `--standalone`
wraps it in a minimal page so it can be opened on its own.

Drawing conventions:
  + strand above the axis, - strand below; genes as arrows labelled by
  gene_symbol (falling back to the feature key)
  essential (copy_num) = --built, other called CDS = --muted,
  mature peptides = --warn
  spliced features: exons joined by a thin line
  segmented genomes: one axis per contig, all drawn on ONE shared nt->px
  scale, so a short segment looks short

Three things this gets right that the adenovirus original did not, each
because it broke on the first family that was not an adenovirus:

1. **A feature key without the module prefix.** A `non_pssm_partner` feature
   is assigned as a bare `CP_NP`, not `Picobirnaviridae.CP_NP.3`, so
   `fa[0][1].split('.')[1]` raises IndexError. Both shapes are accepted.
2. **Multi-contig genomes.** The original summed every contig's length into
   one axis, which for a bisegmented family draws the two segments as one
   false 4.3 kb molecule with a seam in the middle.
3. **The palette.** The coverage-audit page defines --built/--missed/--warn/
   --muted; the adenovirus page defined --mod. A figure referencing a variable
   the page does not define renders with the browser's initial value and the
   chart is invisible on screen while every QC step passes -- exactly the
   failure the CSS comment at gen_coverage_audit.py:136 was written for.
"""
import json
import argparse
import html


def feature_key(fa, module):
    """LOWVAN family assignment -> feature key in the module JSON."""
    tag = fa[0][1]
    if tag.startswith(module + "."):
        tag = tag[len(module) + 1:]
    #  strip the profile index: CP_MAT.6_1 -> CP_MAT, RDRP.1_1 -> RDRP
    return tag.split(".")[0]


def load_gto(path):
    return json.load(open(path))


def read_fasta(path):
    seqs, cid, buf = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if cid:
                seqs.append((cid, "".join(buf)))
            cid, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if cid:
        seqs.append((cid, "".join(buf)))
    return seqs


def gto_from_tbl(tbl, fasta):
    """A GTO-shaped dict from annotate_by_viral_pssm.pl's .feature.tbl.

    The feature table is the cheap path: it comes straight out of the
    annotator with no network round trip, where make_gto.py registers a new
    genome id with the ID server for every genome it builds. For a figure
    that is the wrong trade -- 70-odd registrations to draw 70 pictures.

    What the table cannot carry is exon structure: a spliced feature appears
    as the single span the PSSM matched, because the splice pass is a
    separate script that writes its own output. A module with spliced
    features should be mapped from a GTO produced by that full pass; for
    everything else the table and the GTO agree exactly.

    Columns (0-based): 0 genome id, 1 domain, 2 contig, 3 src, 4 type,
    5 feature id, 6 gene symbol, 7 begin, 8 end, 9 strand, 10 length,
    11 module, 12 rep contig, 13 bits, 14 genome ids, 15 genome name,
    16 pssm, 17 function.
    """
    contigs = [{"id": cid, "dna": seq} for cid, seq in read_fasta(fasta)]
    feats, gname, gid = [], "", ""
    for line in open(tbl):
        p = line.rstrip("\n").split("\t")
        if len(p) < 18:
            continue
        gid, gname = p[0], p[15]
        b, e, strand = int(p[7]), int(p[8]), p[9]
        lo, hi = min(b, e), max(b, e)
        #  pssm column -> feature key: Picobirnaviridae.CP_MAT.2_1.pssm, or a
        #  bare CP_NP for a non_pssm_partner that has no profile at all.
        tag = p[16][:-5] if p[16].endswith(".pssm") else p[16]
        feats.append({
            "type": p[4],
            "function": p[17],
            #  GTO convention: the coordinate is the feature's 5' end and the
            #  length runs 3'-ward, so a minus-strand feature starts at hi.
            "location": [[p[2], str(lo if strand == "+" else hi), strand,
                          hi - lo + 1]],
            "family_assignments": [["LOWVAN", tag, p[17], "tbl"]],
        })
    return {"id": gid, "scientific_name": gname,
            "contigs": contigs, "features": feats}


def tick_step(n):
    for s in (100, 200, 500, 1000, 2000, 5000, 10000, 20000):
        if n / s <= 12:
            return s
    return 50000


#  ---------------------------------------------------------------------------
#  Circular replicons
#
#  A linear axis cannot draw a circular genome honestly. A gene crossing the
#  origin folds to a two-piece join(), and on a straight line those pieces land
#  at opposite ends with nothing connecting them: the Orthohepadnavirus map drew
#  P, L and M twice each, and M's 9 nt first piece came out 1.3 px wide against
#  the right margin. The reader sees two genes where there is one, and the very
#  thing these modules were built to express is the thing the figure hides.
#
#  On a ring the same feature is one unbroken arc that happens to pass twelve
#  o'clock, which is what it is. The origin is drawn as a marked spoke rather
#  than left implicit, because it carries no biology -- it is wherever the
#  submitter cut the circle, and for hepadnaviruses that is a restriction site
#  that most strains no longer have.
#  ---------------------------------------------------------------------------
import math


def circ_pieces(loc, L):
    """[(lo, hi), ...] on the circle, with an origin-crossing join rejoined.

    fold_circular_features.pl emits a wrapped CDS as two pieces, one ending at
    L and one starting at 1. Those are not two exons; they are one run of
    sequence reported either side of the cut. Rejoin them into a single arc and
    leave everything else -- including a genuine splice -- as separate pieces.
    """
    ps = []
    for c, s, st, n in loc:
        s, n = int(s), int(n)
        ps.append((s, s + n - 1) if st == "+" else (s - n + 1, s))
    if len(ps) == 2:
        a, b = sorted(ps)
        if a[0] == 1 and b[1] == L:          # b runs into the origin, a out of it
            return [(b[0], b[1] + (a[1] - a[0] + 1))], True
    return sorted(ps), False


def circ_overlap(s1, n1, s2, n2, L):
    return ((s2 - s1) % L) < n1 or ((s1 - s2) % L) < n2


def orf_groups(ring, L):
    """Group features by the stop codon they share -> one group per reading frame.

    Genes in one frame sharing one stop have the same 3' coordinate, and that
    is the robust test on a circle: absolute position mod 3 is NOT the frame,
    because crossing the origin shifts it whenever the replicon length is not
    a multiple of three. HBV's surface ORF is the case in point -- L at 2850,
    M at 3174 and S at 157 are one frame, but 3182 % 3 == 2, so their start
    coordinates fall in two different residue classes. Their ends do not:
    all three stop at 837.
    """
    g = {}
    for f in ring:
        end = ((f["s"] + f["n"] - 2) % L) + 1 if f["strand"] == "+" else f["s"]
        g.setdefault((f["strand"], end), []).append(f)
    #  Biggest ORF first, so the frames carrying most of the genome take the
    #  strong hues and a lone accessory gene takes the neutral. Gene count
    #  ranks first, then nucleotides covered -- without that second term HBV's
    #  polymerase, one gene over 78% of the circle, ties with X at 465 nt and
    #  loses on coordinate. Coordinate breaks what is left, so the assignment
    #  is stable across runs.
    return sorted(g.values(),
                  key=lambda m: (-len(m), -sum(f["n"] for f in m), m[0]["s"]))


#  Four slots, one per reading frame, validated as a set rather than chosen.
#  Light: worst all-pairs CVD dE 9.9, worst normal-vision 19.6.
#  Dark:  worst all-pairs CVD dE 8.7, worst normal-vision 18.7.
#  (OKLab x100; >=8 CVD target, >=15 normal-vision floor, all 6 pairs in play
#  because arcs on a ring have no fixed neighbour order.)
#
#  Seven hues, one per gene, is what this figure wants and cannot have: at
#  seven categories no ordering clears the floors -- a light-tint set measured
#  worst-CVD dE 0.3, and a saturated set 1.1. Identity is carried by the
#  direct label on every arc, which is required anyway; colour groups the
#  genes by the frame they are read in, which is the thing a reader cannot
#  get from the labels.
ORF_LIGHT = ["#2a78d6", "#eb6834", "#1baf7a", "#c3c2b7"]
ORF_DARK  = ["#3987e5", "#d95926", "#25c48d", "#cfcec5"]


def draw_circular(G, F, a, cid, L, feats, COL):
    """One contig, declared circular, as a ring. Returns (body, width, height)."""
    CX, R = 360, 186
    BAND, LANE_STEP = 17, 22
    TOPPAD = 86                      # legend + origin label
    CY = TOPPAD + R + 26

    def pt(nt, r):
        t = 2 * math.pi * ((nt - 1) % L) / L - math.pi / 2
        return CX + r * math.cos(t), CY + r * math.sin(t)

    def band(s_nt, n_nt, r0, r1):
        frac = min(n_nt / L, 0.9999)
        big = 1 if frac > 0.5 else 0
        x0o, y0o = pt(s_nt, r1)
        x1o, y1o = pt(s_nt + n_nt, r1)
        x1i, y1i = pt(s_nt + n_nt, r0)
        x0i, y0i = pt(s_nt, r0)
        return ("M%.1f,%.1f A%.1f,%.1f 0 %d 1 %.1f,%.1f "
                "L%.1f,%.1f A%.1f,%.1f 0 %d 0 %.1f,%.1f Z"
                % (x0o, y0o, r1, r1, big, x1o, y1o,
                   x1i, y1i, r0, r0, big, x0i, y0i))

    ring = []
    for f in feats:
        ps, wrapped = circ_pieces(f["loc"], L)
        ring.append(dict(f, pieces=ps, wrapped=wrapped, s=ps[0][0],
                         n=sum(hi - lo + 1 for lo, hi in ps)))

    #  ---- colour by reading frame ------------------------------------------
    groups = orf_groups(ring, L)
    for gi, members in enumerate(groups):
        slot = gi if gi < 3 else 3          # never cycle; the rest go neutral
        for f in members:
            f["slot"] = slot
            f["group"] = members

    #  ---- lane packing ------------------------------------------------------
    lanes = {"+": [], "-": []}
    for f in sorted(ring, key=lambda f: -f["n"]):
        ls = lanes[f["strand"]]
        for i, used in enumerate(ls):
            if all(not circ_overlap(f["s"], f["n"], u[0], u[1], L) for u in used):
                used.append((f["s"], f["n"])); f["lane"] = i; break
        else:
            ls.append([(f["s"], f["n"])]); f["lane"] = len(ls) - 1
    nout, nin = max(len(lanes["+"]), 1), len(lanes["-"])
    R_LAB = R + nout * LANE_STEP + 30

    body = []
    body.append('<circle cx="%d" cy="%d" r="%d" fill="none" '
                'stroke="currentColor" stroke-width="1.5" opacity="0.4"/>'
                % (CX, CY, R))

    step = tick_step(L)
    k = 0
    while k < L:
        x0, y0 = pt(max(k, 1), R - 5); x1, y1 = pt(max(k, 1), R + 5)
        body.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" '
                    'stroke="currentColor" opacity="0.4"/>' % (x0, y0, x1, y1))
        if k:
            tx, ty = pt(k, R - 17)
            body.append('<text x="%.1f" y="%.1f" text-anchor="middle" dy="3.5" '
                        'font-size="11" fill="currentColor" opacity="0.7">%s</text>'
                        % (tx, ty, "%.1f" % (k / 1000.0) if step < 1000
                           else "%d" % (k // 1000)))
        k += step

    body.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="currentColor" '
                'stroke-width="1.5" stroke-dasharray="4 3"/>'
                % (CX, CY - R - 30, CX, CY - R + 12))
    body.append('<text x="%d" y="%.1f" text-anchor="middle" font-size="12" '
                'fill="currentColor">origin &middot; nt 1</text>'
                % (CX, CY - R - 38))

    #  ---- arcs ---------------------------------------------------------------
    placed = []
    for f in sorted(ring, key=lambda f: (f["strand"], f["lane"])):
        col = "var(--orf%d)" % f["slot"]
        out = f["strand"] == "+"
        off = 6 + f["lane"] * LANE_STEP
        r0 = (R + off) if out else (R - off - BAND)
        r1 = r0 + BAND
        if f["wrapped"]:
            body.append('<path d="%s" style="fill:%s"/>'
                        % (band(f["s"], f["n"], r0, r1), col))
        else:
            for lo, hi in f["pieces"]:
                body.append('<path d="%s" style="fill:%s"/>'
                            % (band(lo, hi - lo + 1, r0, r1), col))
            if len(f["pieces"]) > 1:
                rm = (r0 + r1) / 2
                for (l0, h0), (l1, _) in zip(f["pieces"], f["pieces"][1:]):
                    p0 = pt(h0, rm); p1 = pt(l1, rm)
                    body.append('<path d="M%.1f,%.1f A%.1f,%.1f 0 0 1 %.1f,%.1f" '
                                'fill="none" stroke="%s" stroke-width="1.2"/>'
                                % (p0[0], p0[1], rm, rm, p1[0], p1[1], col))
        end = (f["s"] + f["n"]) if out else f["s"]
        ex0, ey0 = pt(end, r0 - 3); ex1, ey1 = pt(end, r1 + 3)
        body.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" '
                    'stroke="%s" stroke-width="2.5"/>' % (ex0, ey0, ex1, ey1, col))
        mid = f["s"] + f["n"] / 2.0
        placed.append((mid, r1 if out else r0, f, col))

    #  ---- labels outside the ring, on leader lines ---------------------------
    #  Inside the band the text is small, clipped by the arc ends and sitting on
    #  a saturated fill; outside it is page ink on the page background, which is
    #  legible in both themes and needs no contrast compromise from the fill.
    items = []
    R_ELB = R + nout * LANE_STEP + 14
    for mid, rr, f, col in placed:
        t = 2 * math.pi * ((mid - 1) % L) / L - math.pi / 2
        ax, ay = pt(mid, rr + 3)
        ex_, ey_ = pt(mid, R_ELB)
        right = math.cos(t) >= 0
        #  The label's y is a COLUMN position, not the elbow's. Taking it
        #  radially puts a gene at six o'clock a full R_LAB below centre --
        #  on duck HBV that threw preS and P to y=553 and 572 on a 552-high
        #  canvas. Clamping to the ring's own vertical extent keeps both
        #  columns beside the circle, and the leader's last leg takes up the
        #  slack as a short diagonal.
        ky = min(max(CY + (R_ELB + 6) * math.sin(t), CY - R + 4), CY + R - 4)
        items.append([ax, ay, ex_, ey_, ky, right, f, col])

    #  Each label is TWO lines, the symbol and the length, so the gap has
    #  to clear the whole block and not just one baseline. At 20 the
    #  symbol of the next label landed on the previous one's aa count.
    MINGAP = 38
    for right in (True, False):
        col_items = sorted([i for i in items if i[5] == right], key=lambda i: i[4])
        for i in range(1, len(col_items)):
            if col_items[i][4] - col_items[i - 1][4] < MINGAP:
                col_items[i][4] = col_items[i - 1][4] + MINGAP
        #  Pushing down can run the last one off the bottom; lift the column
        #  back and re-space upward so the gap survives the correction.
        over = col_items[-1][4] - (CY + R + 36) if col_items else 0
        if over > 0:
            for it in col_items:
                it[4] -= over
            for i in range(len(col_items) - 2, -1, -1):
                if col_items[i + 1][4] - col_items[i][4] < MINGAP:
                    col_items[i][4] = col_items[i + 1][4] - MINGAP

    lab_lo = min(i[4] for i in items) if items else CY
    lab_hi = max(i[4] for i in items) if items else CY
    for ax, ay, ex_, ey_, ky, rt, f, col in items:
        cx_lab = CX + (R_ELB + 30) * (1 if rt else -1)
        tx = cx_lab + (7 if rt else -7)
        body.append('<polyline points="%.1f,%.1f %.1f,%.1f %.1f,%.1f" '
                    'fill="none" stroke="%s" stroke-width="1.4" opacity="0.9"/>'
                    % (ax, ay, ex_, ey_, cx_lab, ky, col))
        body.append('<circle cx="%.1f" cy="%.1f" r="2.8" fill="%s"/>' % (ax, ay, col))
        #  The two lines straddle the leader's endpoint rather than hanging
        #  below it, so the leader points at the middle of the block. 18 px
        #  between baselines clears 13.5 px type with room to read.
        body.append('<text x="%.1f" y="%.1f" dy="-3" text-anchor="%s" '
                    'font-size="13.5" font-weight="600" fill="currentColor">%s</text>'
                    % (tx, ky, "start" if rt else "end", html.escape(f["sym"])))
        body.append('<text x="%.1f" y="%.1f" dy="15" text-anchor="%s" '
                    'font-size="10.5" fill="currentColor" opacity="0.65">%d aa</text>'
                    % (tx, ky, "start" if rt else "end", (f["n"] // 3) - 1))

    body.append('<text x="%d" y="%d" text-anchor="middle" font-size="15" '
                'font-weight="600" fill="currentColor">%s nt</text>'
                % (CX, CY - 2, format(L, ",")))
    body.append('<text x="%d" y="%d" text-anchor="middle" font-size="11.5" '
                'fill="currentColor" opacity="0.7">circular</text>'
                % (CX, CY + 16))

    #  ---- legend: one row per reading frame ----------------------------------
    leg, ly = [], 20
    for gi, members in enumerate(groups):
        slot = gi if gi < 3 else 3
        names = "/".join(m["sym"] for m in sorted(members, key=lambda m: -m["n"]))
        lx = 4 + (gi % 2) * 330
        yy = ly + (gi // 2) * 20
        leg.append('<rect x="%d" y="%d" width="15" height="11" rx="2" '
                   'style="fill:var(--orf%d)"/>' % (lx, yy - 9, slot))
        leg.append('<text x="%d" y="%d" font-size="12" fill="currentColor">'
                   '%s</text>' % (lx + 21, yy, html.escape(names)))
    leg.append('<text x="4" y="%d" font-size="11" fill="currentColor" '
               'opacity="0.7">one colour per reading frame &middot; '
               'outside the ring = forward strand</text>'
               % (ly + ((len(groups) + 1) // 2) * 20 + 2))

    H = max(CY + R + max(nin, 1) * LANE_STEP + 46, lab_hi + 40)
    return leg + body, 720, int(H)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gto", help="annotated .gto for ONE genome")
    ap.add_argument("--tbl", help="annotate_by_viral_pssm.pl .feature.tbl instead "
                                  "of a GTO; needs --fasta for contig lengths")
    ap.add_argument("--fasta", help="the contigs that --tbl was called on")
    ap.add_argument("--module-json", required=True)
    ap.add_argument("--module", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--standalone", action="store_true")
    ap.add_argument("--caption", default="",
                    help="extra sentence appended to the figcaption")
    a = ap.parse_args()

    if not a.gto and not a.tbl:
        ap.error("one of --gto or --tbl is required")
    if a.tbl and not a.fasta:
        ap.error("--tbl needs --fasta for the contig lengths")
    G = load_gto(a.gto) if a.gto else gto_from_tbl(a.tbl, a.fasta)
    F = json.load(open(a.module_json))[a.module]["features"]
    SEG = json.load(open(a.module_json))[a.module].get("segments") or {}

    contigs = [(c["id"], len(c["dna"])) for c in G["contigs"]]
    clen = dict(contigs)
    LMAX = max(clen.values())

    X0, XW = 72, 880
    CH, BOX, LANE = 6.2, 10, 25

    def sx(nt):
        return X0 + XW * (nt - 1) / LMAX

    #  ---- collect features, grouped by contig -------------------------------
    per = {cid: [] for cid, _ in contigs}
    for f in G.get("features", []):
        fa = f.get("family_assignments")
        if not fa:
            continue
        key = feature_key(fa, a.module)
        spec = F.get(key, {})
        exons, cid = [], None
        for c, s, st, n in f["location"]:
            cid = c
            s, n = int(s), int(n)
            exons.append((s, s + n - 1) if st == "+" else (s - n + 1, s))
        if cid not in per:
            continue
        kind = ("mat" if f["type"] == "mat_peptide"
                else ("ess" if spec.get("copy_num") else "oth"))
        per[cid].append(dict(
            key=key, sym=spec.get("gene_symbol") or key, exons=sorted(exons),
            strand=f["location"][0][2], kind=kind, loc=f["location"],
            lo=min(x for x, _ in exons), hi=max(y for _, y in exons)))

    #  ---- is this a single circular replicon? --------------------------------
    #  Only the one-contig case draws as a ring. A multi-segment circular family
    #  would need one ring per segment on a shared scale, and nothing in the kit
    #  has one yet; falling back to the linear axis is wrong but visible, where
    #  a half-written ring would be wrong and look finished.
    def seg_for(cid):
        Lc = clen[cid]
        hits = [w for w in SEG.values()
                if isinstance(w, dict)
                and w.get("min_len", 0) <= Lc <= w.get("max_len", 10 ** 9)]
        return hits[0] if len(hits) == 1 else {}

    is_circ = (len(contigs) == 1
               and str(seg_for(contigs[0][0]).get("replicon_geometry", "")
                       ).lower() == "circular")

    #  ---- lane packing, per contig per strand -------------------------------
    for cid, feats in per.items():
        lanes = {"+": [], "-": []}
        for f in sorted(feats, key=lambda f: (f["kind"] == "mat", f["lo"])):
            x1, x2 = sx(f["lo"]), sx(f["hi"])
            cx, lw = (x1 + x2) / 2, len(f["sym"]) * CH
            occ = (min(x1, cx - lw / 2) - 3, max(x2, cx + lw / 2) + 3)
            ls = lanes[f["strand"]]
            for i, used in enumerate(ls):
                if all(occ[1] < u0 or occ[0] > u1 for u0, u1 in used):
                    used.append(occ)
                    f["lane"] = i
                    break
            else:
                ls.append([occ])
                f["lane"] = len(ls) - 1
        per[cid] = (feats, len(lanes["+"]), len(lanes["-"]))

    COL = {"ess": "var(--built)", "oth": "var(--muted)", "mat": "var(--warn)"}

    def arrow(x1, x2, y, strand, col, tip):
        h = BOX if col != COL["mat"] else BOX - 4
        y0 = y - h / 2
        t = min(5, (x2 - x1) * 0.6) if tip else 0
        if strand == "+":
            pts = [(x1, y0), (x2 - t, y0), (x2, y), (x2 - t, y0 + h), (x1, y0 + h)]
        else:
            pts = [(x1 + t, y0), (x2, y0), (x2, y0 + h), (x1 + t, y0 + h), (x1, y)]
        return ('<polygon points="%s" style="fill:%s"/>'
                % (" ".join("%.1f,%.1f" % p for p in pts), col))

    #  ---- segment naming from the module JSON -------------------------------
    def seg_label(cid):
        L = clen[cid]
        hits = [s for s, w in SEG.items()
                if isinstance(w, dict)
                and w.get("min_len", 0) <= L <= w.get("max_len", 10 ** 9)]
        return hits[0] if len(hits) == 1 else ""

    #  ---- draw --------------------------------------------------------------
    TOP = 50                       # legend band
    body, y = [], TOP
    if is_circ:
        cid, L = contigs[0]
        body, VBW, H = draw_circular(G, F, a, cid, L, per[cid][0], COL)
    else:
        VBW = 980
    for cid, L in ([] if is_circ else contigs):
        feats, npl, nml = per[cid]
        ax_y = y + npl * LANE + 20
        for f in feats:
            fy = (ax_y - 18 - f["lane"] * LANE) if f["strand"] == "+" \
                else (ax_y + 30 + f["lane"] * LANE)
            col, ex = COL[f["kind"]], f["exons"]
            if len(ex) > 1:
                body.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" '
                            'style="stroke:%s" stroke-width="1"/>'
                            % (sx(ex[0][1]), fy, sx(ex[-1][0]), fy, col))
            last = ex[-1] if f["strand"] == "+" else ex[0]
            for p, q in ex:
                body.append(arrow(sx(p), sx(q) + 1, fy, f["strand"], col,
                                  (p, q) == last))
            cx = (sx(f["lo"]) + sx(f["hi"])) / 2
            ty = (fy - BOX / 2 - 3) if f["strand"] == "+" else (fy + BOX / 2 + 11)
            body.append('<text x="%.1f" y="%.1f" text-anchor="middle" '
                        'font-size="10.5" fill="currentColor">%s</text>'
                        % (cx, ty, html.escape(f["sym"])))

        body.append('<line x1="%d" y1="%.1f" x2="%.1f" y2="%.1f" '
                    'stroke="currentColor" stroke-width="1.5"/>'
                    % (X0, ax_y, sx(L), ax_y))
        step = tick_step(LMAX)
        k = 0
        while k <= L:
            x = sx(max(k, 1))
            body.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" '
                        'stroke="currentColor"/>' % (x, ax_y - 4, x, ax_y + 4))
            body.append('<text x="%.1f" y="%.1f" text-anchor="middle" '
                        'font-size="10.5" fill="currentColor">%s</text>'
                        % (x, ax_y + 17, "%.1f kb" % (k / 1000.0) if step < 1000
                           else "%d kb" % (k // 1000)))
            k += step
        lab = seg_label(cid)
        body.append('<text x="4" y="%.1f" font-size="11" fill="currentColor">'
                    '<tspan font-weight="600">%s</tspan></text>'
                    % (ax_y - 2, html.escape(lab or cid)))
        body.append('<text x="4" y="%.1f" font-size="10" fill="currentColor" '
                    'opacity="0.75">%s</text>'
                    % (ax_y + 12, html.escape("%s nt" % format(L, ","))
                       if lab else html.escape("%s nt" % format(L, ","))))
        y = ax_y + 34 + nml * LANE + 16

    if not is_circ:
        H = y

    leg, lx = [], 4
    for kind, label in ([] if is_circ else (("ess", "essential CDS"), ("oth", "other called CDS"),
                        ("mat", "mature peptide"))):
        leg.append('<rect x="%d" y="14" width="18" height="10" style="fill:%s"/>'
                   % (lx, COL[kind]))
        leg.append('<text x="%d" y="23" font-size="11" fill="currentColor">%s</text>'
                   % (lx + 24, label))
        lx += 24 + len(label) * 6.6 + 26

    name = G.get("scientific_name") or (a.gto or a.tbl)
    gid = G.get("id", "")
    nf = sum(len(per[c][0]) for c, _ in contigs)
    ne = sum(1 for c, _ in contigs for f in per[c][0] if f["kind"] == "ess")
    nm = sum(1 for c, _ in contigs for f in per[c][0] if f["kind"] == "mat")
    total = sum(clen.values())
    unit = ("%d segments, %s nt total" % (len(contigs), format(total, ","))
            if len(contigs) > 1 else "%s nt" % format(total, ","))
    claim = ("%s (%s, %s) as the %s module annotates it: %d features "
             "(%d essential, %d mature peptides)."
             % (name, gid, unit, a.module, nf, ne, nm))

    tail = (" All segments are drawn on one scale, so their relative sizes are "
            "real." if len(contigs) > 1 else "")
    if a.caption:
        tail += " " + a.caption

    if is_circ:
        nw = sum(1 for f in per[contigs[0][0]][0]
                 if circ_pieces(f["loc"], contigs[0][1])[1])
        where = ("Drawn to scale from this genome&rsquo;s own calls; "
                 "forward-strand genes outside the ring, reverse-strand inside, "
                 "a tick at each 3&prime; end. Colour marks the reading frame, "
                 "not the gene &mdash; genes sharing a stop codon share a hue, "
                 "so the nested forms are visible as one frame. ")
        where += (("%d gene%s cross%s the origin, each drawn as one unbroken arc "
                   "through it. " % (nw, "" if nw == 1 else "s",
                                     "es" if nw == 1 else ""))
                  if nw else "")
        where += ("The origin is not biology: it is wherever the submitter cut "
                  "the circle.")
    else:
        where = ("Drawn to scale from this genome&rsquo;s own calls; "
                 "forward-strand genes above the axis, reverse-strand below.")

    #  The four reading-frame colours travel INSIDE the <figure>. The page
    #  that embeds this defines --built/--muted/--warn and nothing else, and a
    #  fill referencing an undefined custom property paints nothing -- the
    #  figure is present in the HTML and invisible on screen, the failure this
    #  file's docstring already records once. A <style> placed BEFORE the
    #  <figure> does not survive: gen_coverage_audit.py lifts the figure
    #  element alone, which drops the block and leaves 32 fills pointing at
    #  nothing. Measured, not assumed -- the first attempt shipped that way.
    orfcss = ""
    if is_circ:
        L4 = ";".join("--orf%d:%s" % (i, c) for i, c in enumerate(ORF_LIGHT))
        D4 = ";".join("--orf%d:%s" % (i, c) for i, c in enumerate(ORF_DARK))
        orfcss = ('<style>.lvmap{%s}'
                  '@media (prefers-color-scheme:dark){'
                  ':root:not([data-theme="light"]) .lvmap{%s}}'
                  ':root[data-theme="dark"] .lvmap{%s}</style>' % (L4, D4, D4))

    fig = ('<figure class="lvmap">' + orfcss + '<svg viewBox="0 0 %d %d" role="img" aria-label="%s">'
           '%s%s</svg><figcaption>%s %s%s</figcaption></figure>'
           % (VBW, H, html.escape(claim), "".join(leg), "".join(body),
              html.escape(claim), where, tail))

    if a.standalone:
        fig = ('<!doctype html><meta charset="utf-8"><title>%s genome map</title>'
               '<style>:root{--built:#2a78d6;--muted:#898781;--warn:#b07d18;'
               'color-scheme:light dark}'
               '@media (prefers-color-scheme:dark){:root{--built:#3987e5;'
               '--muted:#898781;--warn:#d3a445}}'
               'body{margin:0;padding:24px;font:14px/1.5 system-ui,sans-serif;'
               'background:#f9f9f7;color:#0b0b0b}'
               '@media (prefers-color-scheme:dark){body{background:#0d0d0d;'
               'color:#fff}}'
               'figure{margin:0}figure svg{width:100%%;height:auto;display:block}'
               'figcaption{margin-top:10px;color:#52514e}</style>'
               % html.escape(a.module)) + fig

    open(a.out, "w").write(fig)
    print("%s  %d features on %d contig(s), height %d"
          % (a.out, nf, len(contigs), H))


if __name__ == "__main__":
    main()
