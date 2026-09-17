#!/usr/bin/env python3
"""Bin the Matonaviridae BV-BRC dump into per-feature collections.

Seven features, from the literature (not from database annotation):

  p200  ->  p150 + p90        one cleavage, Gly1301/Gly1302
            PMID 9557742   single site, two products P150/P90 (Liu 1998)
            PMID 10823845  "two mature products ... NH2-p150-p90-COOH" (Liang 2000)
            PMID 10799588  scissile bond by mutagenesis, G1301S abolishes (Liang 2000)
  p110  ->  C + E2 + E1       host signal peptidase
            PMID 2273395   p110 -> C/E2/E1, signalase cleaves E2-E1 (Oker-Blom 1990)
            PMID 2214022   E2 signal seq REMAINS on capsid C-terminus (Suomalainen 1990)
            PMID 11160697  same, and required for assembly (Law 2001)

Because both signal sequences are retained on the upstream product, the three
structural products tile p110 with no residues removed -- which is why no _SP
feature is declared and why SignalP is not needed anywhere in this build.

Binning is STRING-FIRST, length-second.  Length cannot separate p150 from the
structural polyprotein: rubella p150 is 1301 aa but Rustrela p150 is 1082 and
Ruhugu p110 is 1098, so the windows overlap in 1000-1200.  Only the annotation
string distinguishes them, so the generic strings are resolved by explicit
'structural' vs 'protease/p150' tests rather than by length alone.

Each feature gets three files:
  <FEAT>.fasta           length inside the window -- the only file the pipeline reads
  <FEAT>.partial.fasta   too short (dominated by the WHO E1 genotyping amplicon)
  <FEAT>.outliers.fasta  too long
"""
import collections, os, re, sys

WORK = os.path.dirname(os.path.abspath(__file__))
MODULE = 'Matonaviridae'
OUT = os.path.join(WORK, 'collections', MODULE)

# ---------------------------------------------------------------- length windows
# Covering all three rubivirus species:
#            rubella   Rustrela     Ruhugu
# p200         2116    1910/1911      2048
# p150         1301    1082/1083      1237
# p90           815         828        811
# p110         1063    1146/1147      1098
# C             300     332/333        317
# E2            282         324        295
# E1            481         490        486
WINDOW = {
    'P200': (1700, 2250),
    'P150': (1000, 1400),
    'P90':  ( 700,  900),
    'P110':  (1000, 1250),
    'C':    ( 270,  360),
    'E2':   ( 250,  345),
    'E1':   ( 440,  510),
}

# ------------------------------------------------------- species not in the module
# Divergent rubi-like viruses from metatranscriptomic surveys.  One or two
# records each, mostly single fragmentary contigs, and no two of them share a
# measurable gene layout with each other or with the rubiviruses.  A profile
# cannot be built from a singleton even at the leftover floor of -m 2, so these
# are recorded rather than modelled.  Step 10 checks that no rep contig routes
# them: an excluded taxon must be REJECTED, not captured and dead-ended.
NOT_MODELLED = {
    'Neolamprologus multifasciatus matonavirus': 'cichlid; NS 1607 + ST 874, no C/E2/E1 cassette',
    'Trematocara marginatum matonavirus':        'cichlid; NS 385 + ST 293, fragmentary',
    'Aulonocranus dewindti matonavirus':         'cichlid; single 475 aa fragment',
    'Callochromis macrops matonavirus':          'cichlid; single 122 aa fragment',
    'Lamprologus lemairii matonavirus':          'cichlid; single 166 aa fragment',
    'Tetronarce matonavirus':                    'torpedo ray; ST 1001 + RdRp 399',
    'Tiger flathead matonavirus':                'fish; single 1590 aa RdRp',
    'Eaulepmac virus':                           'fish; RdRp 1075 + capsid 316',
    'Mislepmac virus':                           'fish; RdRp 1702 + 2 hypothetical',
    'Rumple rubivirus':                          'ORF2 947 + polyprotein 943, layout unresolved',
    'Ruche rubivirus':                           'single 155 aa RdRp fragment',
    'Ruffle rubivirus':                          'NS 1171 + ST 1054, 1 genome',
    'Guangdong chinese water snake matonavirus': 'reptile; NS 1710 + ST 687',
    'Matonaviridae sp.':                         '5 records, 97-126 aa fragments only',
}

def excluded(name):
    for k in NOT_MODELLED:
        if name.startswith(k.rsplit(' ', 1)[0]) and k.split()[0] in name:
            return k
    for k in NOT_MODELLED:
        if name.startswith(k):
            return k
    return None

# ----------------------------------------------------------------------- the rules
# (feature, compiled test on the lowercased annotation string)
# Order matters: the first match wins, so specific tests precede generic ones.
def rules():
    R = []
    A = lambda f, pat: R.append((f, re.compile(pat, re.I)))

    # --- p90: name the polymerase/helicase product before any generic polymerase
    A('P90',  r'\bp90\b')
    A('P90',  r'rna-(directed|dependent) rna polymerase.*(triphosphatase|helicase)')

    # --- p150: the protease product
    A('P150', r'\bp150\b')
    A('P150', r'nonstructural protease|non-structural protease')

    # --- p200: the nonstructural precursor
    A('P200', r'\bp200\b')
    A('P200', r'non[- ]?structural (poly)?protein')
    A('P200', r'nonstructural poly')
    A('P200', r'nonstructural proteins? (nsp1|p150)')
    A('P200', r'^orf1')

    # --- structural precursor
    A('P110',  r'structural polyprotein|structural protein precursor|polyprotein precursor')
    A('P110',  r'structural proteins c, e2 and e1')
    A('P110',  r'rubella virus polyprotein')
    A('P110',  r'virion protein')
    A('P110',  r'^orf2')

    # --- the three structural products
    A('C',    r'capsid')
    A('E2',   r'\be2\b')
    A('E1',   r'\be1\b|glycoprotien e1')

    # --- generic strings, resolved by length below
    A('?GEN', r'envelope|glycoprotein|hemagglutinin|spike|polyprotein|'
              r'structural protein|nonstructural protein|rna-dependent rna polymerase|'
              r'\brdrp\b|hypothetical|unnamed')
    return R

RULES = rules()

# --------------------------------------------------- resolving the generic strings
# A quarter of the dump's strings name no protein: 'envelope glycoprotein',
# 'structural protein', 'polyprotein', 'hemagglutinin', 'NS4 protein'.  Length
# cannot resolve them either, because the 200-299 band holds BOTH mature E2
# (282 aa) and the 246-aa WHO E1 genotyping amplicon, which is the single most
# common record in the whole dump.
#
# So they are resolved by HOMOLOGY: blastp against the sequences that the
# specific strings already binned confidently, and take the best hit.  Every
# call is written to GENERIC_RESOLVED.tsv with its bitscore so the decision is
# auditable rather than silent.  'NS4 protein' is the U-number trap in this
# taxon -- rubella encodes no NS4, and the label is positional noise.
MIN_BITS = 60.0

# ------------------------------------------------------------- corrupt records
# BV-BRC holds a block of rubella records whose protein translation is in the
# wrong frame.  Two shapes, and only the first is self-announcing:
#
#   1. internal stop codons -- 116 sequences, up to 8 stops each
#   2. no stops, full length, but the sequence is not the protein.  E1 cluster 3
#      was seven such records beginning GFHLPLYCTGVRHADT where every genuine
#      rubella E1 begins EEAFTYLCTAPGCATQ.
#
# Both built PSSMs in the first pass: E1.lo1 from five 8-stop sequences and
# E1.3 from the frameshifted seven.  Each profile recognised its own members at
# ~1000 bits and NOTHING else -- the signature of a junk cluster -- and they
# were precisely the twelve sequences the dominant E1 profile could not see.
#
# So shape 2 is caught by homology against the feature's own modal-length core:
# genuine divergent species score far above the floor (Rustrela E1 513 bits,
# Ruhugu 596), while a wrong-frame translation scores nothing at all.
MIN_CORE_BITS = 100.0

def has_internal_stop(s):
    return '*' in s.rstrip('*')

def drop_corrupt(recs, workdir):
    """Remove wrong-frame translations. Returns (kept, dropped_rows)."""
    import subprocess, tempfile, statistics
    dropped = []
    # shape 1: internal stops
    keep1 = []
    for row in recs:
        if row[1] and has_internal_stop(row[6]):
            dropped.append((row[2], row[3], row[1], len(row[6]), 'internal stop codon'))
        else:
            keep1.append(row)
    # shape 2: no similarity to the feature's modal-length core
    byfeat = collections.defaultdict(list)
    for row in keep1:
        if row[1]:
            byfeat[row[1]].append(row)
    keepset = {id(r) for r in keep1}
    tmp = tempfile.mkdtemp(prefix='corrupt.')
    for feat, rows in byfeat.items():
        lo, hi = WINDOW[feat]
        full = [r for r in rows if lo <= len(r[6]) <= hi]
        if len(full) < 20:
            continue
        mode = statistics.mode([len(r[6]) for r in full])
        core = [r for r in full if len(r[6]) == mode]
        if len(core) < 5:
            continue
        ref = f'{tmp}/{feat}.ref'
        with open(ref, 'w') as fh:
            for j, r in enumerate(core[:40]):
                fh.write(f">c{j}\n{r[6]}\n")
        qry = f'{tmp}/{feat}.q'
        with open(qry, 'w') as fh:
            for j, r in enumerate(rows):
                fh.write(f">q{j}\n{r[6]}\n")
        subprocess.run(['makeblastdb', '-in', ref, '-dbtype', 'prot'],
                       check=True, capture_output=True)
        out = subprocess.run(
            ['blastp', '-query', qry, '-db', ref, '-outfmt', '6 qseqid bitscore',
             '-evalue', '10', '-num_threads', '8', '-comp_based_stats', '0',
             '-max_target_seqs', '5'],
            check=True, capture_output=True, text=True).stdout
        best = {}
        for line in out.splitlines():
            q, b = line.split('\t'); b = float(b)
            if b > best.get(q, 0):
                best[q] = b
        for j, r in enumerate(rows):
            if best.get(f'q{j}', 0.0) < MIN_CORE_BITS:
                keepset.discard(id(r))
                dropped.append((r[2], r[3], r[1], len(r[6]),
                                f'no similarity to {feat} core '
                                f'({best.get(f"q{j}",0):.0f} bits, floor {MIN_CORE_BITS:.0f})'))
    kept = [r for r in keep1 if id(r) in keepset]
    with open(f'{workdir}/collections/CORRUPT_DROPPED.tsv', 'w') as fh:
        fh.write("representative_id\tgenome_name\tfeature\taa_len\treason\n")
        for row in sorted(dropped, key=lambda r: (r[2] or '', r[4])):
            fh.write('\t'.join(str(x) for x in row) + '\n')
    return kept, dropped

def resolve_by_homology(pending, confident, workdir):
    """pending: [(idx, rep, name, ann, seq)]   confident: {feat: [(hdr,seq)]}"""
    import subprocess, tempfile
    if not pending:
        return {}, []
    tmp = tempfile.mkdtemp(prefix='genres.')
    ref = os.path.join(tmp, 'ref.faa')
    with open(ref, 'w') as fh:
        for feat, recs in confident.items():
            for j, (_h, s) in enumerate(recs):
                fh.write(f">{feat}__{j}\n{s}\n")
    qry = os.path.join(tmp, 'q.faa')
    with open(qry, 'w') as fh:
        for i, rep, nm, a, s in pending:
            fh.write(f">q{i}\n{s}\n")
    subprocess.run(['makeblastdb', '-in', ref, '-dbtype', 'prot'],
                   check=True, capture_output=True)
    out = subprocess.run(
        ['blastp', '-query', qry, '-db', ref, '-max_target_seqs', '5',
         '-outfmt', '6 qseqid sseqid bitscore pident length', '-evalue', '1e-5',
         '-num_threads', '8', '-comp_based_stats', '0'],
        check=True, capture_output=True, text=True).stdout
    best = {}
    for line in out.splitlines():
        q, s, bits, pid, alen = line.split('\t')
        bits = float(bits)
        if bits < MIN_BITS:
            continue
        if q not in best or bits > best[q][1]:
            best[q] = (s.split('__')[0], bits, float(pid), int(alen))
    assigned, log = {}, []
    for i, rep, nm, a, s in pending:
        hit = best.get(f'q{i}')
        if hit:
            assigned[i] = hit[0]
            log.append((rep, nm, a, len(s), hit[0], f"{hit[1]:.0f}",
                        f"{hit[2]:.1f}", str(hit[3])))
        else:
            log.append((rep, nm, a, len(s), '-', '-', '-', '-'))
    with open(f'{workdir}/collections/GENERIC_RESOLVED.tsv', 'w') as fh:
        fh.write("representative_id\tgenome_name\tannotation\taa_len\t"
                 "assigned_feature\tbitscore\tpident\tali_len\n")
        for row in sorted(log, key=lambda r: (r[4], -int(r[3]))):
            fh.write('\t'.join(str(x) for x in row) + '\n')
    return assigned, log

# --------------------------------------------------------------------------- load
def load():
    md5 = [l.rstrip('\n') for l in open(f'{WORK}/Matona.uniq.md5')]
    ann, rep = [], []
    for l in open(f'{WORK}/Matona.uniq.id_ann'):
        p = l.rstrip('\n').split('\t')
        rep.append(p[0]); ann.append(p[1] if len(p) > 1 else '')
    seq = []
    for l in open(f'{WORK}/Matona.uniq.seq'):
        p = l.rstrip('\n').split('\t')
        seq.append(p[1] if len(p) > 1 else '')
    gname, ggenus = {}, {}
    for l in open(f'{WORK}/Matona.id_name_gs'):
        p = l.rstrip('\n').split('\t')
        gname[p[0]] = p[1]
        ggenus[p[0]] = p[3] if len(p) > 3 else ''
    # every genome carrying each md5
    carriers = collections.defaultdict(set)
    for l in open(f'{WORK}/Matona.id_md5'):
        fid, m = l.rstrip('\n').split('\t')
        g = re.match(r'fig\|(\d+\.\d+)', fid)
        if g:
            carriers[m].add(g.group(1))
    return md5, ann, rep, seq, gname, ggenus, carriers

def main():
    md5, ann, rep, seq, gname, ggenus, carriers = load()
    os.makedirs(OUT, exist_ok=True)

    notmod = collections.Counter()
    recs = []      # (i, feat_or_None, rep, name, genus, ann, seq)

    # ------------------------------------------------- pass 1: the specific strings
    for i, m in enumerate(md5):
        a, s, r = ann[i], seq[i], rep[i]
        if not s:
            continue
        gids = sorted(carriers.get(m, []))
        names = [gname.get(g, '') for g in gids]
        # a sequence is excluded only if EVERY genome carrying it is excluded
        exc = [excluded(n) for n in names]
        if names and all(exc):
            notmod[exc[0]] += 1
            continue
        nm = next((n for n, e in zip(names, exc) if not e), names[0] if names else '')
        genus = next((ggenus.get(g, '') for g in gids), '')

        feat = None
        for f, pat in RULES:
            if pat.search(a):
                feat = f
                break
        if feat == '?GEN':
            feat = None                       # pass 2 decides
        recs.append([i, feat, r, nm, genus, a, s])

    # ------------------------------------------------ pass 2: homology for the rest
    confident = collections.defaultdict(list)
    for i, feat, r, nm, genus, a, s in recs:
        if feat:
            lo, hi = WINDOW[feat]
            if lo <= len(s) <= hi:
                confident[feat].append((r, s))
    pending = [(i, r, nm, a, s) for i, feat, r, nm, genus, a, s in recs if not feat]
    assigned, _log = resolve_by_homology(pending, confident, WORK)
    n_resolved = len(assigned)
    for row in recs:
        if row[1] is None:
            row[1] = assigned.get(row[0])

    # ------------------------------------------------- drop wrong-frame records
    recs, corrupt = drop_corrupt(recs, WORK)

    # ---------------------------------------------------------------- classify/bin
    bins = collections.defaultdict(list)
    unassigned = []
    syn = collections.defaultdict(collections.Counter)
    species_seen = collections.defaultdict(collections.Counter)
    for i, feat, r, nm, genus, a, s in recs:
        L = len(s)
        if not feat:
            unassigned.append((r, nm, a, L))
            continue
        lo, hi = WINDOW[feat]
        cls = 'full' if lo <= L <= hi else ('partial' if L < lo else 'outliers')
        hdr = f">{r} [{nm}] {genus or 'Rubivirus'} | {a}"
        bins[(feat, cls)].append((hdr, s))
        syn[feat][a] += 1
        if cls == 'full':
            species_seen[feat][' '.join(nm.split()[:2])] += 1

    # ------------------------------------------------------------------- write out
    for (feat, cls), recs_out in sorted(bins.items()):
        name = f"{feat}.fasta" if cls == 'full' else f"{feat}.{cls}.fasta"
        with open(f'{OUT}/{name}', 'w') as fh:
            for h, s in recs_out:
                fh.write(h + '\n')
                for j in range(0, len(s), 60):
                    fh.write(s[j:j+60] + '\n')

    with open(f'{WORK}/collections/synonyms.tsv', 'w') as fh:
        fh.write("feature\tn_uniq_seqs\tbvbrc_annotation_string\n")
        for f in WINDOW:
            for a, n in syn[f].most_common():
                fh.write(f"{f}\t{n}\t{a}\n")

    # The string-collapse artifact needs counts weighted by FEATURE, not by
    # unique sequence -- a protein carried identically by 200 genomes is one
    # unique sequence and 200 features, and the collapse claim is about how
    # many records each controlled string absorbs.  collect_synmap.py reads
    # this from the workdir root with these exact column names.
    md5_feat = {}
    for i, feat, r, nm, genus, a, s in recs:
        md5_feat[md5[i]] = (feat, a)
    per = collections.defaultdict(int)
    for l in open(f'{WORK}/Matona.id_md5'):
        fid, m = l.rstrip('\n').split('\t')
        if m not in md5_feat:
            continue
        feat, a = md5_feat[m]
        g = re.match(r'fig\|(\d+\.\d+)', fid)
        gen = ggenus.get(g.group(1), '') if g else ''
        per[(a, gen or 'Rubivirus', feat or '')] += 1
    with open(f'{WORK}/synonyms.tsv', 'w') as fh:
        fh.write("count\tannotation\tgenus\tbins_to\n")
        for (a, gen, feat), n in sorted(per.items(), key=lambda kv: -kv[1]):
            fh.write(f"{n}\t{a}\t{gen}\t{feat}\n")

    with open(f'{WORK}/collections/UNASSIGNED_TRACKING.tsv', 'w') as fh:
        fh.write("representative_id\tgenome_name\tannotation\taa_len\n")
        for r, nm, a, L in sorted(unassigned, key=lambda x: -x[3]):
            fh.write(f"{r}\t{nm}\t{a}\t{L}\n")

    with open(f'{WORK}/collections/NOT_MODELLED.tsv', 'w') as fh:
        fh.write("species_prefix\tn_uniq_seqs\treason\n")
        for k, why in NOT_MODELLED.items():
            fh.write(f"{k}\t{notmod.get(k,0)}\t{why}\n")

    # ---------------------------------------------------------------------- report
    print(f"{'feat':<6}{'full':>7}{'partial':>9}{'outlier':>9}   window      species in the full set")
    print("-" * 96)
    tot = 0
    for f in WINDOW:
        nf = len(bins.get((f, 'full'), []))
        np_ = len(bins.get((f, 'partial'), []))
        no = len(bins.get((f, 'outliers'), []))
        tot += nf
        sp = ', '.join(f"{k.split()[-1] if k.startswith('Rubivirus') else k.split()[0]}:{v}"
                       for k, v in species_seen[f].most_common(4))
        lo, hi = WINDOW[f]
        print(f"{f:<6}{nf:>7}{np_:>9}{no:>9}   {lo}-{hi:<6}  {sp}")
    print("-" * 96)
    print(f"{'TOTAL':<6}{tot:>7}")
    print(f"\ngeneric strings resolved by homology: {n_resolved}/{len(pending)}"
          f"  -> collections/GENERIC_RESOLVED.tsv")
    nstop = sum(1 for c in corrupt if c[4] == 'internal stop codon')
    print(f"corrupt records dropped: {len(corrupt)} "
          f"({nstop} internal stop, {len(corrupt)-nstop} no core similarity)"
          f"  -> collections/CORRUPT_DROPPED.tsv")
    print(f"unassigned : {len(unassigned)}  -> collections/UNASSIGNED_TRACKING.tsv")
    print(f"not modelled: {sum(notmod.values())} uniq seqs over {len(notmod)} excluded species"
          f"  -> collections/NOT_MODELLED.tsv")
    print(f"synonyms    : {sum(len(v) for v in syn.values())} string->feature mappings"
          f"  -> collections/synonyms.tsv")

if __name__ == '__main__':
    main()
