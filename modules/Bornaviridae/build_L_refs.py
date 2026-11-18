#!/usr/bin/env python3
"""Build the Bornaviridae L splice references, and round-trip every one.

A reference is the UNSPLICED genomic region from the exon-1 ATG to the stop,
with the donor and acceptor marked in the header exactly as
get_splice_variant_features.pl parses them:

    >genome_id SD:<win_start>-<win_end>;<last nt of exon 1> SA:<win_start>-<win_end>;<first nt of exon 2>

All coordinates in the header are 1-based into the reference sequence.
Internally everything here is 0-based and converted once, at write time,
because the first version of this script mixed the two and produced
references whose SD window did not contain the donor.

Every reference is round-tripped before it is written: cut it at its own
declared coordinates, translate, and require a clean Met-initiated ORF of
the expected length. A reference that cannot reproduce its own protein is
dropped rather than shipped.
"""
import os, re, sys, collections

CODE = {}
_B = 'TCAG'
_A = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
_i = 0
for _a in _B:
    for _b in _B:
        for _c in _B:
            CODE[_a + _b + _c] = _A[_i]; _i += 1

def tr(s):
    return ''.join(CODE.get(s[i:i+3], 'X') for i in range(0, len(s) - 2, 3))

def rc(s):
    return s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]

def genome(gid):
    return ''.join(l.strip() for l in open('Contigs/%s.fna' % gid) if l[0] != '>').upper()

SD_BEFORE, SD_AFTER = 8, 2     # window around the donor:  8 nt of exon 1 + the GT
SA_BEFORE, SA_AFTER = 12, 2    # window around the acceptor: the AG + 10 nt, + 2 nt of exon 2

def make_ref(seq, s0, d0, a0):
    """s0 = first nt of exon 1, d0 = last nt of exon 1, a0 = first nt of exon 2 (all 0-based).
    Returns (reference_seq, header_fields) or None."""
    spliced = seq[s0:d0+1] + seq[a0:]
    prot = tr(spliced)
    stop = prot.find('*')
    if stop < 1500 or stop > 1800:
        return None
    nt_total = 3 * (stop + 1)
    ex1 = d0 - s0 + 1
    gen = seq[s0:a0 + (nt_total - ex1)]
    sd_last = d0 - s0 + 1          # 1-based within gen
    sa_first = a0 - s0 + 1         # 1-based within gen
    sd_s, sd_e = sd_last - SD_BEFORE, sd_last + SD_AFTER
    sa_s, sa_e = sa_first - SA_BEFORE, sa_first + SA_AFTER - 1
    if sd_s < 1 or sa_e > len(gen) or sa_s < 1:
        return None
    #  ---- round trip, exactly as the caller does it ----
    left = gen[:sd_last]
    right = gen[sa_first-1:]
    rt = tr(left + right)
    rt = rt[:rt.find('*')] if '*' in rt else rt
    if not rt.startswith('M') or len(rt) != stop or '*' in rt or 'X' in rt:
        return None
    hdr = 'SD:%d-%d;%d SA:%d-%d;%d' % (sd_s, sd_e, sd_last, sa_s, sa_e, sa_first)
    return gen, hdr, rt

def solve_spliced(seq, approx_e2, want_start='MSFH'):
    """Find (s0,d0,a0) with a canonical GT..AG that yields a clean L."""
    sols = []
    for a0 in range(max(0, approx_e2 - 150), approx_e2 + 150):
        if seq[a0-2:a0] != 'AG':
            continue
        for d0 in range(max(0, a0 - 1500), a0 - 1000):
            if seq[d0+1:d0+3] != 'GT':
                continue
            for n1 in (17, 18, 16, 15, 19, 20, 14, 21):
                s0 = d0 + 1 - n1
                if s0 < 0 or seq[s0:s0+3] != 'ATG':
                    continue
                r = make_ref(seq, s0, d0, a0)
                if r and r[2].startswith(want_start):
                    sols.append((abs((a0 - d0 - 1) - 1294), abs(len(r[2]) - 1711), s0, d0, a0, r))
    if not sols:
        return None
    sols.sort(key=lambda x: (x[0], x[1]))
    return sols[0][2], sols[0][3], sols[0][4], sols[0][5]


def main():
    meta = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open('borna.full.tsv')}
    seqs = {}; h = None
    for l in open('collections/Bornaviridae/L.fasta'):
        if l[0] == '>': h = l[1:].strip().split()[0]; seqs[h] = ''
        else: seqs[h] += l.strip()

    rows = [l.rstrip('\n').split('\t') for l in open('splice2/junctions.tsv')][1:]
    refs = []; stats = collections.Counter()

    #  1. the orthobornavirus junctions solved exactly against their own deposit
    for r in rows:
        fid, genus, gid, strand = r[0], r[1], r[2], r[3]
        if r[4] in ('COLLINEAR', 'NOT_SOLVED'):
            continue
        s0, d0, a0 = int(r[4]) - 1, int(r[5]) - 1, int(r[6]) - 1
        seq = genome(gid)
        if strand == '-': seq = rc(seq)
        out = make_ref(seq, s0, d0, a0)
        if not out: stats['ortho round-trip failed'] += 1; continue
        gen, hdr, prot = out
        if prot != seqs[fid]: stats['ortho product != deposit'] += 1; continue
        refs.append((gid, hdr, gen, 'spliced', genus)); stats['spliced from own deposit'] += 1

    #  2. the genera where L is collinear: a zero-length intron inside the ORF
    for r in rows:
        fid, genus, gid, strand = r[0], r[1], r[2], r[3]
        if r[4] != 'COLLINEAR' or genus == 'Orthobornavirus':
            continue
        p = seqs[fid]; seq = genome(gid)
        if strand == '-': seq = rc(seq)
        st = None
        for fr in range(3):
            t = tr(seq[fr:]); idx = t.find(p[:40])
            if idx < 0: continue
            c = fr + 3 * idx
            if tr(seq[c:c + 3 * (len(p) + 1)]) == p + '*': st = c; break
        if st is None: stats['collinear not placed'] += 1; continue
        gen = seq[st:st + 3 * (len(p) + 1)]
        sd_last, sa_first = 111, 112                       # pseudo-junction, zero-length intron
        hdr = 'SD:%d-%d;%d SA:%d-%d;%d' % (sd_last - SD_BEFORE, sd_last + SD_AFTER, sd_last,
                                           sa_first - SA_BEFORE, sa_first + SA_AFTER - 1, sa_first)
        rt = tr(gen[:sd_last] + gen[sa_first-1:]); rt = rt[:rt.find('*')]
        if rt != p: stats['collinear round-trip failed'] += 1; continue
        refs.append((gid, hdr, gen, 'collinear', genus)); stats['collinear (zero-length intron)'] += 1

    #  3. lineages with no usable deposit: solve from the genome against an anchor
    have = {g for g, _h, _s, _k, _gn in refs}
    extra = [g for g in sys.argv[1:] if g not in have]
    anchor = max((p for p in seqs.values() if len(p) > 1600), key=len)
    for gid in extra:
        seq0 = genome(gid)
        done = False
        for seq in (seq0, rc(seq0)):
            idx = -1
            for fr in range(3):
                t = tr(seq[fr:]); i2 = t.find(anchor[200:240])
                if i2 >= 0: idx = fr + 3 * (i2 - 200 + 6); break
            if idx < 0: continue
            sol = solve_spliced(seq, idx)
            if not sol: continue
            s0, d0, a0, (gen, hdr, prot) = sol
            refs.append((gid, hdr, gen, 'spliced', meta.get(gid, ['', '', '', ''])[3]))
            stats['spliced, solved from the genome'] += 1
            print('   %-12s solved: exon1 %d nt, intron %d nt, L = %d aa, starts %s'
                  % (gid, d0 - s0 + 1, a0 - d0 - 1, len(prot), prot[:8]))
            done = True; break
        if not done:
            print('   %-12s NO SOLUTION' % gid); stats['unsolved'] += 1

    seen = set(); out = []
    for gid, hdr, gen, kind, genus in refs:
        if gen in seen: continue
        seen.add(gen); out.append((gid, hdr, gen, kind, genus))
    os.makedirs('Splice-Variants/Bornaviridae', exist_ok=True)
    with open('Splice-Variants/Bornaviridae/L.fasta', 'w') as f:
        for gid, hdr, gen, kind, genus in out:
            f.write('>%s %s\n%s\n' % (gid, hdr, gen))
    for k, v in sorted(stats.items()): print('  %-34s %d' % (k, v))
    print('  %-34s %d (%d unique sequences)' % ('references written', len(out), len(seen)))
    sd = collections.Counter(); sa = collections.Counter()
    for gid, hdr, gen, kind, genus in out:
        m = re.search(r'SD:(\d+)-(\d+);(\d+) SA:(\d+)-(\d+);(\d+)', hdr)
        a, b, c, d, e, g2 = map(int, m.groups())
        sd[gen[a-1:b]] += 1; sa[gen[d-1:e]] += 1
    print('  distinct SD motifs %d, distinct SA motifs %d' % (len(sd), len(sa)))
    print('  commonest donor   ', sd.most_common(3))
    print('  commonest acceptor', sa.most_common(3))

if __name__ == '__main__':
    main()
