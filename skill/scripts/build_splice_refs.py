#!/usr/bin/env python3
"""Build Splice-Variants references for a one-intron CDS, from a query protein.

For a protein the literature describes but the dump does not call -- the ones
SKILL.md step 2b turns up -- there is no collection to cluster, so there is no
PSSM to build. If its CDS crosses exactly ONE intron it is instead a splice
feature, and get_splice_variant_features.pl can call it from references of the
shape Splice-Variants/README.md defines.

"Doubly spliced" in a paper is about the TRANSCRIPT. Count only the junctions
the CODING sequence crosses: a leader junction upstream of the start codon
does not count, and the caller models one intron, not one splice event.

    MMTV Rem   301 aa  = SP + C-term            one CDS intron   PMID 16282474
    PFV  Bet   482 aa  = Bel-1[1..88] + Bel-2   one CDS intron   PMID 1846194

Method, which is also the verification: tblastn the query protein against each
genome, take the two plus-strand HSPs that run in query order without
overlapping on the genome, nudge the exon boundary within +/-3 nt onto a
canonical GT..AG, then TRANSLATE the spliced product. A reference is kept only
if it starts at M, carries no internal stop, and lands within 10% of the query
length. Everything else is dropped and counted -- the drops are the evidence
that the junction is real, not a BLAST artefact.

    python3 build_splice_refs.py --query bet.faa --contigs Contigs.Spuma \
            --out Spumaretrovirinae/BET.fasta --name BET
"""
import argparse, glob, os, subprocess, sys, tempfile

B = 'TCAG'
AAs = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
TAB = {a + b + c: AAs[i] for i, (a, b, c) in enumerate(
    (x, y, z) for x in B for y in B for z in B)}


#  The near-cognate initiators seen in this genus. ACG, CTG, GTG and TTG all
#  initiate at reduced efficiency; the submitter writes the product as if it
#  began at Met, which is why a deposit and its own genome disagree on
#  residue 1 and on nothing else.
ALT_STARTS = {"ATG", "CTG", "TTG", "GTG", "ACG"}


def tr(s):
    return "".join(TAB.get(s[i:i + 3], 'X') for i in range(0, len(s) - 2, 3))


def tblastn(query, subject, evalue, seg="no"):
    r = subprocess.run(["tblastn", "-query", query, "-subject", subject,
                        "-outfmt", "6 qstart qend sstart send bitscore",
                        "-seg", seg, "-evalue", str(evalue)],
                       capture_output=True, text=True)
    return [list(map(float, l.split("\t")))
            for l in r.stdout.strip().split("\n") if l]


def head_hsp(qseq, seq, upto, hi, tmpd):
    """Place the FIRST exon when tblastn cannot see it.

    Exon A can be short -- PFV Bet opens with 88 residues, feline foamy Bet
    with about 30 -- and a 30-residue exon scores below any sane evalue in a
    whole-genome tblastn, so the protein looks unspliced and is dropped. Here
    the N-terminal residues the main HSP did not cover are re-searched, alone,
    against only the sequence upstream of exon B, where a permissive threshold
    is safe because the window is small and the query is anchored.
    """
    if upto < 8 or hi < 60:
        return None
    qf = os.path.join(tmpd, "head.faa")
    sf = os.path.join(tmpd, "up.fna")
    open(qf, "w").write(">h\n%s\n" % qseq[:upto])
    open(sf, "w").write(">u\n%s\n" % seq[:hi])
    rows = [r for r in tblastn(qf, sf, 100) if r[3] > r[2]]
    if not rows:
        return None
    rows.sort(key=lambda x: -x[4])
    return rows[0]


def product_is_homologous(prot, qfile, qlen, tmpd, min_id=30.0, min_cov=0.60):
    """Reject a junction whose spliced product is not really the query protein.

    A canonical GT..AG, a clean translation and a plausible length are still
    not enough. On equine foamy virus this builder found a junction that
    spliced 37 residues of upstream sequence onto the whole of Bel-2: the
    product began with M, ran to its own stop and sat inside the length
    window, yet residues 82-410 were 100% identical to that genome's own
    called BEL2 and nothing in the first 81 matched Bel-1. It was not Bet.

    Alignment REGISTER does not separate these -- the artefact still aligns
    from near its N-terminus. Identity does, and with room to spare. Measured
    over the 78 candidate BET references: the seven artefacts (two equine
    foamy, one equine, four simian) score 43-56 bits at 20-27% identity,
    while every genuine reference scores 273-1030 bits at 35.8-99.8%. There
    is nothing in between, so the floors below sit in open space rather than
    on top of the distribution.
    """
    pf = os.path.join(tmpd, "prod.faa")
    open(pf, "w").write(">p\n%s\n" % prot)
    r = subprocess.run(["blastp", "-query", pf, "-subject", qfile,
                        "-outfmt", "6 pident length bitscore", "-evalue", "1e-3"],
                       capture_output=True, text=True)
    rows = [l.split("\t") for l in r.stdout.strip().split("\n") if l]
    if not rows:
        return False
    b = max(rows, key=lambda x: float(x[2]))
    return float(b[0]) >= min_id and int(b[1]) >= min_cov * qlen


def exon_a_starts(seq, exA_s, back=300):
    """Candidate start codons for exon A, nearest first.

    tblastn reports where the query first MATCHES, not where the CDS begins.
    On the divergent simian foamy genomes only the last ~30 residues of Bet's
    88-residue first exon align, so the HSP opens a couple of dozen codons
    inside the protein and the spliced translation does not start with M --
    the genome was then dropped as "internal stop" even though its junction
    was right. This is the mirror of the rule below that ends the CDS at its
    own stop codon: the CDS also STARTS at its own start codon. Walk upstream
    in frame and offer every ATG reached before an in-frame terminator, the
    HSP start included. The translation and length checks, not this walk,
    decide which one is the real initiator.
    """
    out = [exA_s]
    i = exA_s - 3
    while i >= 1 and exA_s - i <= back:
        cod = seq[i - 1:i + 2]
        if cod in ("TAA", "TAG", "TGA"):
            break
        if cod == "ATG":
            out.append(i)
        i -= 3
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--query", required=True, help="protein FASTA, one sequence")
    p.add_argument("--contigs", required=True, help="dir of <genome_id>.fna")
    p.add_argument("--out", required=True)
    p.add_argument("--name", required=True, help="feature key, for the report")
    p.add_argument("--min-contig", type=int, default=3000)
    p.add_argument("--len-tol", type=float, default=0.10,
                   help="fractional length window around the query (default 0.10)")
    p.add_argument("--max-pairs", type=int, default=6,
                   help="how many HSP pairs to try per genome, best exon B "
                        "first (default 6). More than one is needed because a "
                        "spurious short HSP at the query's N-terminus "
                        "otherwise claims the exon-A slot.")
    p.add_argument("--alt-start", action="store_true",
                   help="accept a near-cognate initiator (CTG/TTG/GTG/ACG) as "
                        "well as ATG. Needed for BLV G4; the product still has "
                        "to be homologous to the query.")
    p.add_argument("--noncanonical", action="store_true",
                   help="if no GT..AG pair sits in the window, choose the "
                        "boundary by the spliced translation instead. The "
                        "caller checks the SD/SA windows, never the "
                        "dinucleotide.")
    p.add_argument("--nudge", type=int, default=12,
                   help="how far, in nt, the exon boundary may be moved to reach "
                        "a canonical GT..AG (default 12). tblastn trims an HSP to "
                        "the last well-matching codon, so the reported boundary "
                        "sits a few codons short of the true donor; +/-3 was too "
                        "tight and silently lost most foamy virus genomes. The "
                        "translation check below, not this window, is the guard.")
    a = p.parse_args()

    qseq = "".join(l.strip() for l in open(a.query) if not l.startswith(">"))
    qlen = len(qseq)
    tmpd = tempfile.mkdtemp(prefix="splicerefs.")
    refs = []
    stats = {"no2hsp": 0, "nojunction": 0, "stop": 0, "len": 0, "nterm": 0}
    for fa in sorted(glob.glob(os.path.join(a.contigs, "*.fna"))):
        gid = os.path.basename(fa)[:-4]
        seq = "".join(l.strip() for l in open(fa) if not l.startswith(">")).upper()
        if len(seq) < a.min_contig:
            continue
        plus = sorted([x for x in tblastn(a.query, fa, 1e-5) if x[3] > x[2]],
                      key=lambda x: x[0])
        if not plus:
            stats["no2hsp"] += 1; continue
        #  EVERY PLAUSIBLE HSP PAIR, NOT THE FIRST ONE IN QUERY ORDER.
        #  `plus` is sorted by query start, so a spurious short HSP that
        #  happens to open at residue 1 is paired with the real exon B and the
        #  genuine exon A is never looked at. On HTLV-1 11908.2880 that put
        #  Rex's first exon at position 1044, inside gag, four kilobases from
        #  the real one at 5125 which a GT..AG junction and a 20-residue exon
        #  fit exactly. The reference translated cleanly and passed the
        #  homology check -- exon B is 163 of its 183 residues, so a wrong
        #  20-residue N-terminus cannot fail it -- and the caller then
        #  reproduced the same wrong join on that genome. 14 of 103 PTLV-1
        #  references carried a mismatched N-terminus, one of them with stop
        #  codons in it.
        pairs = [(plus[i], plus[j])
                 for i in range(len(plus)) for j in range(len(plus))
                 if i != j and plus[i][1] < plus[j][0] + 10
                 and plus[i][3] < plus[j][2]]
        #  best exon B first: it carries most of the protein and is the piece
        #  the junction is measured against.
        pairs.sort(key=lambda t: (-t[1][4], -t[0][4]))
        if not pairs:
            #  one long HSP that does not start at residue 1: the missing
            #  N-terminus is exon A, too short for the genome-wide search.
            y = max(plus, key=lambda v: v[4])
            x = head_hsp(qseq, seq, int(y[0]) - 1, int(y[2]) - 1, tmpd)
            if x is None or not (x[1] < y[0] + 10):
                stats["no2hsp"] += 1; continue
            pairs = [(x, y)]
        hits, why = [], "stop"
        for x, y in pairs[:a.max_pairs]:
            exA_s, exA_e, exB_s = int(x[2]), int(x[3]), int(y[2])

            #  Every GT..AG the window reaches, nearest boundary first, and the
            #  first one whose SPLICED TRANSLATION is a real protein wins. Taking
            #  the first canonical dinucleotide instead picks a junction that is
            #  merely possible over the one that codes.
            cand = sorted(((d, e) for d in range(-a.nudge, a.nudge + 1)
                           for e in range(-a.nudge, a.nudge + 1)
                           if seq[exA_e + d:exA_e + d + 2] == "GT"
                           and seq[exB_s + e - 3:exB_s + e - 1] == "AG"),
                          key=lambda t: (abs(t[0]) + abs(t[1]), abs(t[0])))
            if not cand and a.noncanonical:
                #  GT..AG IS THIS SCRIPT'S RULE, NOT THE CALLER'S.
                #  get_splice_variant_features.pl never looks at the dinucleotide:
                #  it takes the SD and SA windows from the reference header and
                #  requires the subject to match them, and that conservation is
                #  the whole test. So a junction the annotator can use is not
                #  obliged to be canonical, and demanding it here threw away an
                #  entire feature. BLV G4 is annotated in GenBank as
                #  join(431..503, 7071..7315) (BHU32669 on LC927011) -- donor TC,
                #  acceptor CC -- and no GT..AG pair anywhere in that intron
                #  reproduces the protein. With this flag the boundary is chosen
                #  by the SPLICED TRANSLATION instead. product_is_homologous() is
                #  then the only thing standing between a real junction and an
                #  artefact, which it already was: a canonical dinucleotide plus a
                #  clean translation passed 7 artefacts on an earlier feature set,
                #  and identity to the query was what separated them.
                cand = sorted(((d, e) for d in range(-a.nudge, a.nudge + 1)
                               for e in range(-a.nudge, a.nudge + 1)),
                              key=lambda t: (abs(t[0]) + abs(t[1]), abs(t[0])))
            if not cand:
                stats["nojunction"] += 1; continue

            #  COLLECT EVERY WORKING COMBINATION AND KEEP THE BEST, instead of
            #  keeping the first. The start walk offers upstream ATGs as well as
            #  the HSP start, and with --alt-start the HSP start may itself be a
            #  near-cognate initiator -- so a genome can satisfy the checks at two
            #  different starts, giving proteins of different lengths. Taking
            #  whichever came first put a 113 aa G4 next to the published 105 aa
            #  one (GenBank BHU32669) in the same reference set. The product
            #  closest to the query length wins, ties going to an ATG and then to
            #  the boundary nearest the HSPs.
            for exA_try in exon_a_starts(seq, exA_s):
              for d, e in cand:
                ie, ib = exA_e + d, exB_s + e
                #  THE CDS ENDS AT ITS OWN STOP CODON, not at the query's length.
                #  Padding exon B out to qlen*3 instead runs the reference past a
                #  real terminator whenever this genome's copy is a few residues
                #  shorter than the query, and the genuine junction is then thrown
                #  out as "internal stop". That lost most of the foamy virus
                #  genomes for BET, whose Bet runs 482-490 aa across isolates.
                head = seq[exA_try - 1:ie]
                tail = seq[ib - 1:]
                off = (3 - len(head) % 3) % 3          # nt of the split codon in B
                prot = tr(head + tail)
                cut = prot.find("*")
                if cut < 0:
                    why = "len"; continue
                nt = (cut + 1) * 3 - len(head)         # nt taken from exon B
                if nt <= off or ib + nt - 1 > len(seq):
                    why = "len"; continue
                ref = seq[exA_try - 1:ib + nt - 1]
                sd_last = ie - exA_try + 1
                sa_first = ib - exA_try + 1
                prot = prot[:cut]
                #  A NON-AUG INITIATOR IS NOT A BROKEN REFERENCE.
                #  Requiring the product to start at M is this script's rule and
                #  nothing downstream shares it: the annotation code handles an
                #  alternative start, and get_splice_variant_features.pl never
                #  inspects the first codon. Near-cognate initiation is described
                #  in this genus (PMID 10710429, HTLV-1 Rex and p21x), and BLV G4
                #  uses it: GenBank BHU32669 gives the protein as MACTRVCFL...
                #  while LC927011 reads CTG at position 431 and translates
                #  LACTRVCFL... -- identical from residue 2. Rejecting it on
                #  residue 1 dropped every G4 in the set before the junction was
                #  even examined. --alt-start admits the near-cognate initiators,
                #  and the product still has to be homologous to the query.
                if not (prot.startswith("M") or
                        (a.alt_start and seq[exA_try - 1:exA_try + 2] in ALT_STARTS)):
                    why = "stop"; continue
                if not ((1 - a.len_tol) * qlen <= len(prot) <= (1 + a.len_tol) * qlen):
                    why = "len"; continue
                #  EXON A FIRST IN THE RANKING. A wrong HSP pair gets exon
                #  B right -- it carries most of the protein -- and exon A
                #  wrong, so product length and overall homology both look
                #  fine. The only term that separates them is how the
                #  exon-A residues compare with the query's own N-terminus.
                nA = sd_last // 3
                agree = sum(1 for u, v in zip(prot[1:nA], qseq[1:nA]) if u == v)
                hits.append(((-agree,
                              abs(len(prot) - qlen),
                              0 if seq[exA_try - 1:exA_try + 2] == "ATG" else 1,
                              abs(d) + abs(e)),
                             (gid, ref, sd_last, sa_first, len(prot)), prot))
        #  HOMOLOGY LAST, ON THE RANKED LIST. product_is_homologous() shells
        #  out to blastp, and --noncanonical offers (2*nudge+1)^2 boundaries
        #  per start instead of the handful that are canonical -- running the
        #  blast inside that loop turned a 13 minute pass into one that had
        #  not finished a twentieth of the contigs in two minutes. The cheap
        #  tests shrink the field first; the blast then confirms the best
        #  candidate and falls through to the next if it fails.
        hit = None
        for _k, cand_hit, prot in sorted(hits, key=lambda t: t[0]):
            if product_is_homologous(prot, a.query, qlen, tmpd):
                hit = cand_hit; break
            why = "nterm"
        if hit is None:
            stats[why] += 1; continue
        refs.append(hit)

    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        for gid, ref, sd, sa, _ in refs:
            #  windows bracket each junction and are matched exactly between
            #  reference and subject -- that is the conservation check.
            fh.write(">%s SD:%d-%d;%d SA:%d-%d;%d\n%s\n"
                     % (gid, sd - 2, sd + 8, sd, sa - 12, sa + 1, sa, ref))
    print("  %s references: %d written -> %s" % (a.name, len(refs), a.out))
    print("  dropped: %d without two usable HSPs, %d without a canonical "
          "junction, %d with an internal stop, %d out of length range, "
          "%d not homologous enough to the query"
          % (stats["no2hsp"], stats["nojunction"], stats["stop"], stats["len"],
             stats["nterm"]))
    if refs:
        L = sorted(r[4] for r in refs)
        print("  query %d aa; spliced product min %d med %d max %d"
              % (qlen, L[0], L[len(L) // 2], L[-1]))
    else:
        #  "more than one intron" was the only explanation this offered, and
        #  it was the wrong one twice in a row. BLV R3 failed on all 20 test
        #  genomes at the default nudge of 12 and succeeded on all 20 at 30:
        #  its first exon is 17 residues and its second 27, and tblastn places
        #  the ends of pieces that short imprecisely, so the real junction sat
        #  outside the window. BLV G4 failed because of this script's own
        #  start-at-M and GT..AG rules, neither of which the caller shares.
        #  Say what to try before blaming the biology.
        sys.exit("  no references built. Before concluding the CDS has more "
                 "than one intron, try:\n"
                 "    --nudge 30      a short exon's HSP ends are imprecise "
                 "and the junction can fall outside a narrow window\n"
                 "    --alt-start     the CDS may begin at a near-cognate "
                 "codon; the deposit still writes it as M\n"
                 "    --noncanonical  the junction need not be GT..AG -- the "
                 "caller checks the SD/SA windows, not the dinucleotide")


main()
