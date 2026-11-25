#!/usr/bin/env python3
"""Mark the inserted bases in a transcript-editing reference set: lowercase
genome, UPPERCASE for the bases the construct adds.

WHY THIS EXISTS
---------------
A `special: transcript_edit` reference is not collinear with the genome. It is
the genome with one or two bases DUPLICATED at the ribosomal slip site (or
inserted at an RNA-editing site), so that the two reading frames can be
translated as one open frame. Which bases those are is the entire content of
the curation, and once the file is written they are invisible: the construct is
just a few thousand indistinguishable nucleotides.

The house convention makes them visible. The genome-derived bases are
lowercase; the added bases are UPPERCASE. Every hand-curated set in the kit
follows it -- the paramyxovirus V/W/P/I sets carry exactly one or two uppercase
characters per sequence, and Orthoebolavirus GPC/Gn/ssGP the same -- so a
curator can find the edit in any reference by eye, or with `grep -o '[ACGT]'`.
Nothing in the pipeline reads case: `get_transcript_edited_features.pl` builds
the edited sequence from blastn's own qseq/sseq, and blastn returns those
uppercase whatever the input file holds (measured, not assumed). So this is
purely for the human, which is exactly why a program has to put it there: nobody will
notice it missing until they need it.

`build_transcript_edit_refs.py` read the genome with `.upper()` and wrote the
construct unchanged, so every set built with it is uniformly uppercase and the
insertion point is unrecoverable from the file alone. That builder now applies
the convention itself. This script repairs the sets written before it did.

HOW THE INSERTION POINT IS RECOVERED
------------------------------------
Not guessed -- re-derived, by one of two exact methods, and verified before
anything is written.

  --contigs   The genome itself. The construct is genome[start:start+q] +
              genome[start+q-n:start+q] + genome[start+q:]. Walk the construct
              and the genome together from the anchored start; the first
              position where they disagree is q. Then the identity above is
              asserted in full: if the remainder of the construct is not the
              remainder of the genome, this is not an n-base duplication and
              the sequence is REFUSED rather than marked.

  --preshift  The deposited pre-shift protein (ORF1a, nsp1a) for the same
              genome, for a set whose contigs are no longer on disk. The fusion
              translation agrees with the pre-shift protein up to the slip and
              diverges after it, which localises q to one codon; the few
              candidates in that window are tested by the duplication identity
              and by requiring that removing them restores the pre-shift
              protein's own terminator.

A sequence that no method confirms is left exactly as it was and counted. A set
is only rewritten when every sequence in it is confirmed, unless --partial is
given: a file where some sequences carry the mark and others do not is worse
than one where none do, because the reader cannot tell which is which.
"""
import argparse, collections, glob, os, re, sys

CODON = {}
_b = "TCAG"
_a = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_i = 0
for _x in _b:
    for _y in _b:
        for _z in _b:
            CODON[_x + _y + _z] = _a[_i]; _i += 1


def tr(nt):
    nt = nt.upper()
    return "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def read_fasta(path):
    h, s = None, []
    for line in open(path):
        if line.startswith(">"):
            if h is not None:
                yield h, "".join(s)
            h, s = line[1:].rstrip("\n"), []
        else:
            s.append(line.strip())
    if h is not None:
        yield h, "".join(s)


def rc(s):
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def gid_of(header):
    """the genome this reference came from, however the header spells it"""
    tok = header.split()[0]
    m = re.search(r"(\d+\.\d+)", tok)
    return m.group(1) if m else tok


def load_contigs(paths):
    """gid -> [contig sequences]. Accepts directories of <gid>.fna and files."""
    idx = collections.defaultdict(list)
    files = []
    for p in paths:
        if os.path.isdir(p):
            files += sorted(glob.glob(os.path.join(p, "*.fna")) +
                            glob.glob(os.path.join(p, "*.fasta")) +
                            glob.glob(os.path.join(p, "*.fa")))
        else:
            files.append(p)
    for f in files:
        base = os.path.basename(f)
        m = re.match(r"(\d+\.\d+)\.", base)
        for h, s in read_fasta(f):
            g = m.group(1) if m else gid_of(h)
            idx[g].append(s.upper())
    return idx


def _check(u, tail, max_n):
    """tail is the genome from the construct's start. Return (q, n) or None."""
    q, lim = 0, min(len(u), len(tail))
    while q < lim and u[q] == tail[q]:
        q += 1
    if q >= len(u):              # collinear -- nothing was added
        return None
    for n in range(1, max_n + 1):
        if q < n:
            continue
        rest = u[q + n:]
        #  the construct must BE the genome with u[q-n:q] repeated
        if u[q:q + n] == u[q - n:q] and rest == tail[q:q + len(rest)]:
            return q, n
    return None


def find_insert_by_genome(ref, contigs, max_n):
    """Return (q, n) or None. q is the offset of the first added base.

    Anchoring is tried several ways because the insertion can sit anywhere,
    including inside the first few codons: the Cardiovirus 2B* construct opens
    at the start of 2B and duplicates a base 30 nt later, so a 40 nt prefix
    anchor spans the insertion and matches nothing. Short prefixes are tried in
    turn, and failing all of them the SUFFIX is anchored instead -- everything
    after the insertion is genome, so the far end always places exactly.
    """
    u = ref.upper()
    for g in contigs:
        for strand in (g, rc(g)):
            for k in (40, 24, 16, 12):
                if k >= len(u):
                    continue
                start = strand.find(u[:k])
                if start < 0:
                    continue
                hit = _check(u, strand[start:], max_n)
                if hit:
                    return hit
            #  suffix anchor: the construct's tail is pure genome
            k = min(40, len(u) // 2)
            end = strand.rfind(u[-k:])
            if end < 0:
                continue
            for n in range(1, max_n + 1):
                start = end + k - (len(u) - n)
                if start < 0:
                    continue
                hit = _check(u, strand[start:], max_n)
                if hit:
                    return hit
    return None


def find_insert_by_preshift(ref, pre, max_n):
    """Place the slip with the deposited pre-shift protein, for a set whose
    genomes are no longer on disk.

    The test is the duplication's defining consequence: delete n bases at q and
    the construct becomes the genome again, so frame 0 must now run into the
    pre-shift gene's OWN terminator, and the peptide it reads up to that stop
    must be part of the deposited pre-shift protein.

    `part of`, not `equal to`: a reference does not have to start where the
    pre-shift gene starts. The Closteroviridae RdRp references open 21 codons
    before the polymerase, deep inside ORF1a, so what frame 0 yields is the
    pre-shift protein's SUFFIX. Requiring equality rejected all 201 of them.

    Every candidate in the sequence is tested rather than a window around a
    guess, and a sequence where more than one position passes is refused as
    ambiguous instead of resolved by preference.
    """
    u = ref.upper()
    pre = pre.rstrip("*")
    hits = []
    for n in range(1, max_n + 1):
        for q in range(n, len(u)):
            if u[q:q + n] != u[q - n:q]:
                continue
            p = tr(u[:q] + u[q + n:])
            stop = p.find("*")
            #  the restored frame must terminate, at or after the slip, on a
            #  peptide the pre-shift deposit actually contains at its end
            if stop < 0 or stop * 3 < q - 3:
                continue
            if stop and pre.endswith(p[:stop]):
                hits.append((q, n))
    return hits[0] if len(hits) == 1 else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refs", required=True, help="reference fasta to mark")
    ap.add_argument("--contigs", nargs="*", default=[],
                    help="directories of <gid>.fna, or fasta files, holding the genomes")
    ap.add_argument("--preshift", default="",
                    help="deposited pre-shift protein fasta, keyed by genome id")
    ap.add_argument("--max-insert", type=int, default=2)
    ap.add_argument("--out", default="", help="default: rewrite --refs in place")
    ap.add_argument("--partial", action="store_true",
                    help="write even when some sequences could not be confirmed")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    refs = list(read_fasta(a.refs))
    contigs = load_contigs(a.contigs) if a.contigs else {}
    pre = {}
    if a.preshift:
        for h, s in read_fasta(a.preshift):
            pre.setdefault(gid_of(h), s.rstrip("*"))

    out, stats, sites, already = [], collections.Counter(), collections.Counter(), 0
    for h, s in refs:
        if any(c.islower() for c in s) and any(c.isupper() for c in s):
            already += 1; out.append((h, s)); stats["already marked"] += 1; continue
        g = gid_of(h)
        hit = None
        if g in contigs:
            hit = find_insert_by_genome(s, contigs[g], a.max_insert)
            stats["confirmed against the genome" if hit else
                  "genome on disk but no duplication explains the construct"] += 1
        elif g in pre:
            hit = find_insert_by_preshift(s, pre[g], a.max_insert)
            stats["confirmed against the pre-shift protein" if hit else
                  "pre-shift protein on disk but the slip could not be placed"] += 1
        else:
            stats["no genome and no pre-shift protein for this id"] += 1
        if hit:
            q, n = hit
            out.append((h, s[:q].lower() + s[q:q + n].upper() + s[q + n:].lower()))
            sites[n] += 1
        else:
            out.append((h, s))

    ok = sum(sites.values()) + already
    print("%s: %d sequence(s)" % (a.refs, len(refs)))
    for k, v in stats.most_common():
        print("   %-58s %d" % (k, v))
    if sites:
        print("   insert sizes: %s" % ", ".join("%d base(s) x%d" % (n, c)
                                                for n, c in sorted(sites.items())))
    if a.dry_run:
        print("   dry run -- nothing written"); return
    if ok < len(refs) and not a.partial:
        print("   NOT WRITTEN: %d of %d unconfirmed. A half-marked file is worse than an "
              "unmarked one -- the reader cannot tell an unmarked sequence from one with "
              "no insertion. Pass --partial to override." % (len(refs) - ok, len(refs)))
        sys.exit(1)
    dest = a.out or a.refs
    with open(dest, "w") as f:
        for h, s in out:
            f.write(">%s\n%s\n" % (h, s))
    print("   wrote %s (%d marked, %d already marked)" % (dest, sum(sites.values()), already))


main()
