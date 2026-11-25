#!/usr/bin/env python3
"""Build a transcript-editing reference set for a frameshifted polyprotein.

A frameshift product is not collinear with the genome, so it carries no PSSM.
`get_transcript_edited_features.pl` blastn's a NUCLEOTIDE reference carrying the
extra base against the genome, fills the gap in the subject from the query and
translates, so each genome gets its own protein rather than a copy of the
reference. References live at
`<LOWVAN_DATA_DIR>/Transcript-Editing/<Module>/<FEAT>.fasta`.

**The construction is self-validating.** Rather than deriving the insertion
point from the slippery heptamer's phase and hoping, it scans candidate
insertion points and keeps only one whose translation reproduces a protein the
genome itself deposits. A genome that cannot reproduce its own protein donates
no reference. That is the skill's "self round trip" applied at build time.

Two validators, because the pool for the first is usually small:

  --fused     the deposited FUSION protein (nsp1ab, ORF1ab). Strongest: the
              construct must equal it exactly.
  --n-part    the deposited pre-shift protein (nsp1a, ORF1a), and
  --c-part    the deposited post-shift protein (ORF1b, the polymerase).
              Used together: the construct must START with the pre-shift
              protein's first --prefix residues, which are shared because both
              products initiate at the same Met, and END with the post-shift
              protein exactly, which is what pins the insertion point.

The second validator matters. Building Astroviridae, Mamastrovirus deposits 217
fused proteins and gains 49 from pairs, but **Avastrovirus deposits 24 and gains
131** -- a six-fold denser reference set, against a correction gate that demands
95% identity.

Locating the protein in the genome is done by exact match in a three-frame
translation, not by tblastn: it is faster, needs no BLAST database per genome,
and an exact match is a stronger anchor than an alignment.

The trailing `*`: a product that really terminates at a stop keeps the stop
codon in the reference (V, W, ORF1ab, ssGP, TF). One that ends at a cleavage
does not (NSP12, Gn). This builder is for the first kind and keeps it.
"""
import argparse, collections, os, sys

CODON = {}
_b = "TCAG"
_a = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_i = 0
for _x in _b:
    for _y in _b:
        for _z in _b:
            CODON[_x + _y + _z] = _a[_i]; _i += 1


def translate(nt):
    return "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def read_fasta(path):
    h, s, out = None, [], []
    for line in open(path):
        if line.startswith(">"):
            if h is not None:
                out.append((h.split()[0], "".join(s)))
            h = line[1:].strip(); s = []
        else:
            s.append(line.strip())
    if h is not None:
        out.append((h.split()[0], "".join(s)))
    return out


def contig(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def locate(genome, protein, probe=30):
    """0-based nt index of the protein's first codon, by exact 3-frame match."""
    frag = protein[:probe]
    for fr in range(3):
        aa = translate(genome[fr:])
        i = aa.find(frag)
        if i >= 0:
            return fr + i * 3
    return -1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--fused", help="FASTA of deposited fusion proteins, id = genome id")
    ap.add_argument("--n-part", help="FASTA of deposited pre-shift proteins")
    ap.add_argument("--c-part", help="FASTA of deposited post-shift proteins")
    ap.add_argument("--out", required=True)
    ap.add_argument("--window", type=int, default=400,  # CODONS back from the seed
                    help="nt downstream of the start to search for the slip (default 400 "
                         "before the pre-shift protein's own stop)")
    ap.add_argument("--c-suffix", type=int, default=300,
                    help="residues of the POST-shift protein the construct must end with "
                         "(0 = all of it). The default avoids depending on where the "
                         "submitter chose to start that CDS.")
    ap.add_argument("--prefix", type=int, default=200,
                    help="residues of the pre-shift protein the construct must start with")
    ap.add_argument("--heptamer", default="",
                    help="comma-separated slippery sequence(s), e.g. TTTAAAC,GTTAAAC. "
                         "When given this is a CONSTRAINT: the duplicated base must lie "
                         "inside one of them, and a genome whose only valid insertion "
                         "sits elsewhere donates no reference. Omit it and the builder "
                         "takes the insertion closest to the pre-shift stop instead.")
    ap.add_argument("--insert", type=int, default=1, choices=(1, 2),
                    help="how many bases the construct duplicates at the slip "
                         "site. 1 = a -1 frameshift (ORF1ab, V/W, TF, ssGP), "
                         "which is what every reference set built before "
                         "Arteriviridae needed. 2 = a -2 frameshift: the "
                         "arterivirus nsp2TF product shifts -2 at G_GUU_UUU "
                         "(PMID 23043113), and a -2 shift reads the same as "
                         "duplicating two bases. Nothing else changes -- the "
                         "construct is still validated against a deposited "
                         "protein, so a wrong insert length simply yields no "
                         "references rather than bad ones.")
    ap.add_argument("--allow-non-met", action="store_true",
                    help="keep constructs whose protein does not start at Met. Only for "
                         "a product whose N-terminus is a cleavage rather than an "
                         "initiation codon (see the CDS-or-mat_peptide note in "
                         "references/special-features.md)")
    a = ap.parse_args()

    if not a.fused and not (a.n_part and a.c_part):
        sys.exit("need --fused, or both --n-part and --c-part")

    heps = [h.strip().upper() for h in a.heptamer.split(",") if h.strip()]
    fused = dict(read_fasta(a.fused)) if a.fused else {}
    npart = dict(read_fasta(a.n_part)) if a.n_part else {}
    cpart = dict(read_fasta(a.c_part)) if a.c_part else {}
    gids = list(fused) + [g for g in sorted(set(npart) & set(cpart)) if g not in fused]

    kept, stats = [], collections.Counter()
    for gid in gids:
        cp = os.path.join(a.contigs, gid + ".fna")
        if not os.path.exists(cp):
            stats["no contig"] += 1; continue
        g = contig(cp)
        anchor = fused.get(gid) or npart.get(gid)
        start = locate(g, anchor)
        if start < 0:
            stats["protein not found in a 3-frame translation"] += 1; continue

        target = fused.get(gid, "").rstrip("*")
        want_c = cpart.get(gid, "").rstrip("*")
        #  Match the TAIL of the post-shift protein, not all of it. Requiring
        #  the whole thing assumes the submitter's post-shift CDS begins
        #  exactly at the junction, and for the arterivirus ORF1b it often
        #  does not -- GenBank ORF1b starts sit a few codons either side of
        #  the slippery site by convention. The tail pins the READING FRAME,
        #  which is what this validator is for; the slip site pins the
        #  POSITION. Together they are stricter than the old test, not looser.
        if want_c and a.c_suffix and len(want_c) > a.c_suffix:
            want_c = want_c[-a.c_suffix:]
        want_n = (npart.get(gid, "") or "")[:a.prefix]

        #  SEED THE SCAN. The original loop walked every insertion point from
        #  the ORF start, translating the whole candidate at each one: O(n^2)
        #  in the ORF length, and bounded by `--window * 10` nucleotides from
        #  the start. That is fine for a 1-2 kb paramyxovirus V/W and useless
        #  for a 12 kb nidovirus pp1ab, where the slip site sits ~7,500 nt in
        #  -- past the default bound, so the build returned no references at
        #  all, and raising the bound made it take hours per hundred genomes.
        #
        #  The insertion point is not unknown. In `--fused` mode the construct
        #  must equal a deposited protein, so the 0-frame translation agrees
        #  with that protein up to the slip and diverges after it: the first
        #  mismatch IS the slip site. In `--n+c` mode the pre-shift product
        #  ends at the 0-frame stop, and the slip is immediately upstream of
        #  it. Either way one cheap translation locates it, and the scan only
        #  has to cover a window around that seed.
        #  An ambiguous base anywhere in the construct disqualifies the
        #  genome, and testing that FIRST keeps the statistics honest. Left
        #  until after the scan it is reported as "no insertion sits in the
        #  slip site", because an N inside the product makes every candidate
        #  translation differ from the deposit by one X -- which is a true
        #  statement about the scan and a misleading one about the genome.
        span = len(target or "") or (len(want_n or "") + len(want_c or ""))
        if span and set(g[start:start + 3 * span + 6]) - set("ACGT"):
            stats["ambiguous base in the genome's own ORF"] += 1; continue
        tr0 = translate(g[start:])
        stop0 = tr0.find("*")
        seed = stop0 if stop0 >= 0 else len(tr0)
        if target:
            for _i in range(min(len(tr0), len(target))):
                if tr0[_i] != target[_i]:
                    seed = min(seed, _i); break
        lo_p = max(start + 3 + a.insert, start + 3 * (seed - a.window))
        hi_p = min(len(g), start + 3 * (seed + 2) + 3)
        #  Collect EVERY valid insertion point, then choose between them.
        #  Taking the first one is wrong in `--n+c` mode, and the failure is
        #  invisible: that mode pins the READING FRAME (the construct must end
        #  with the post-shift protein) but not the POSITION of the shift,
        #  because nothing constrains the residues between the prefix it
        #  checks and the C-terminal part. On Arteriviridae the first valid
        #  point was routinely 16 nt upstream of the real slippery heptamer,
        #  and 811 of 994 PRRSV references were built with the duplication in
        #  the wrong place while reproducing a perfectly good ORF1b. `--fused`
        #  does not have this problem: comparing the whole protein pins the
        #  position exactly.
        cands = []
        for p in range(lo_p, hi_p):
            cand = g[start:p] + g[p - a.insert:p] + g[p:]
            tr = translate(cand)
            stop = tr.find("*")
            if stop < 0:
                continue
            prot = tr[:stop]
            if target:
                ok = prot == target
            else:
                ok = bool(want_c) and prot.endswith(want_c) and \
                     (not want_n or prot.startswith(want_n))
            if ok:
                cands.append((p, cand[:(stop + 1) * 3]))
        hit = None
        if cands:
            if heps:
                #  The slip site is known, so it is a CONSTRAINT, not a
                #  statistic. Keep only insertions whose duplicated base lies
                #  inside the heptamer (or within two bases of its 3' end,
                #  which is the difference between the two join conventions in
                #  circulation, not a difference in biology).
                keep = [(p, c) for p, c in cands
                        if any(h in g[p - len(h) - 2:p + 2] for h in heps)]
                if keep:
                    cands = keep
                else:
                    stats["no insertion sits in the slip site"] += 1
                    continue
            #  Otherwise take the point CLOSEST TO THE PRE-SHIFT STOP: a
            #  ribosomal slip happens just before the 0-frame terminator.
            hit = max(cands, key=lambda x: x[0])
        if hit is None:
            stats["no insertion reproduces the deposited protein"] += 1; continue
        p, ref = hit
        if set(ref) - set("ACGT"):
            stats["ambiguous base in construct"] += 1; continue
        #  A CDS reference whose protein does not begin at Met is wrong, and
        #  the cause is upstream: the deposited protein it was anchored on was
        #  itself N-terminally truncated, so the construct inherits a start
        #  that is not a start. 5-7% of the Astroviridae references were like
        #  this. Refusing them here is the same rule check_met_start.py applies
        #  to profiles -- be suspicious of any non-mat_peptide that does not
        #  start with Met -- applied to a reference instead of an alignment.
        if not a.allow_non_met and translate(ref)[:1] != "M":
            stats["construct does not begin at Met"] += 1; continue
        if heps:
            stats["slip site confirmed at the insertion"] += 1
        stats["kept (%s)" % ("fused" if target else "n+c pair")] += 1
        #  House convention: genome-derived bases lowercase, the bases this
        #  construct ADDS uppercase. Nothing in the pipeline reads case --
        #  blastn returns qseq/sseq uppercase whatever the file holds, and the
        #  caller builds the edited sequence from those -- so the mark is for the
        #  curator, which is precisely why a
        #  program has to write it. Once the file is on disk the duplication is
        #  otherwise unrecoverable without the genome, and every hand-curated
        #  set in the kit (the paramyxovirus V/W/P/I sets, Orthoebolavirus
        #  GPC/Gn/ssGP) carries exactly one or two uppercase characters per
        #  sequence for this reason.
        #
        #  The mark is also a check on this script. In the Astroviridae sets,
        #  built before --heptamer existed, the uppercase base lands 9-18 nt
        #  from the slippery site in most references -- the "valid insertion
        #  point that is not the slip site" failure described above. That is
        #  obvious at a glance in a marked file and invisible in an unmarked
        #  one. skill/scripts/mark_transcript_edit_inserts.py marks sets
        #  written before this.
        i = p - start
        kept.append((gid, ref[:i].lower() + ref[i:i + a.insert].upper()
                     + ref[i + a.insert:].lower()))

    d = os.path.dirname(a.out)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(a.out, "w") as fh:
        for gid, ref in kept:
            fh.write(">%s\n%s\n" % (gid, ref))
    print("  references written: %d -> %s" % (len(kept), a.out))
    for k, v in stats.most_common():
        print("     %-46s %d" % (k, v))
    if kept:
        L = sorted(len(r) // 3 - 1 for _g, r in kept)
        print("  protein length from the references: median %d, range %d-%d"
              % (L[len(L) // 2], L[0], L[-1]))


main()
