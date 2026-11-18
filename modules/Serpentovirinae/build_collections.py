#!/usr/bin/env python3
"""Serpentovirinae collections -- binned from SEQUENCE, not from the names.

Genome organisation, ball python nidovirus (Stenglein 2014, PMID 25205093,
NC_024709, 33,452 nt): 5'-ORF1a-ORF1b-S(ORF2)-ORF3-M(ORF4)-N(ORF5a)-ORF5b-
ORF6-3'. ORF3, ORF5b and ORF6 carry transmembrane segments and glycosylation
sites and have no known function -- that paper says so. ORF1b is reached by a
-1 frameshift at the A-rich heptamer AAAAAAC (NC_024709 misc_feature
18405..18411, CDS join(1018..18411,18411..25343)), so pp1ab is declared
`transcript_edit` and ORF1b is never a feature of its own.

## Why this module does not use the name-driven engine

ORF numbering is not portable across this family -- serpentovirus ORF3 is an
unknown accessory protein while bafinivirus and bostovirus ORF3 is M -- and
within the subfamily one protein arrives as `hypothetical protein` (x10),
`membrane protein` (x6), `putative membrane protein` (x4), `M protein` (x2)
and `putative matrix protein`, while `minor M protein` names a DIFFERENT
protein. A rule set over those strings bins the subfamily wrongly; measured,
it put the real M (13 unique sequences) into an accessory bag and left the M
collection with 4.

So the bins come from homology. All 138 deposited proteins under 2,500 aa are
clustered at 25% identity / 50% coverage, and each cluster is tested against
the three proteins the literature NAMES in NC_024709 -- ORF2 = S, ORF4 = M,
ORF5a = N. Measured (blastp, best hit per cluster):

    S   C1                      e = 0
    M   C2 e=1e-163, C9 e=2e-17 over 163 aa
    N   C3 e=2e-115, C11 e=3e-27 over 157 aa, C13 e=9e-10 over 97 aa,
        C29 e=3e-07 over 106 aa, C33 e=9e-09 over 57 aa
    every other cluster: best e >= 0.01 over <= 72 aa -- i.e. no hit

Genome position was read afterwards and agrees (C13 is the ORF after M in
Sectovirus, C29 the ORF after M in Vebetovirus), but nothing here is called by
position: position only separates one orthogroup from another.

## Accessories

Every cluster with no hit to S/M/N and at least two DISTINCT sequences becomes
`UNC1..UNCn`, an `Uncharacterized lineage-specific protein`. No function is
published for any serpentovirus accessory ORF. UNC7 is the ORF reported to
resemble the paramyxovirus haemagglutinin-neuraminidase (PMID 25063552) --
a similarity, not a demonstrated function, so it keeps the uncharacterized
string like the rest.

Clusters whose members reduce to ONE distinct sequence (one isolate, usually
deposited twice as GenBank + RefSeq) are recorded in NOT_MODELLED.tsv rather
than given a profile.
"""
import collections, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
W = os.path.dirname(HERE)
MODULE, T = "Serpentovirinae", "Serp"
ANCHOR_GENOME = "1986118.12"                 # NC_024709, ball python nidovirus 1
ANCHORS = {"S": (25340, 28219), "M": (28977, 29624), "N": (29637, 30095)}
EVAL, MINALN = 1e-5, 50                      # a cluster joins S/M/N only on a real alignment
MIN_SEQ_ID, COV = "0.25", "0.5"
PARTIAL = ("fragment", "truncated", "n-terminal", "c-terminal", "partial")
REPLICASE = ("orf1a", "orf1ab", "orf1b", "pp1a", "pp1ab", "1ab", "1b protein",
             "polyprotein 1a", "polyprotein 1ab", "replicase 1a", "replicase 1b",
             "pp1a polyprotein", "pp1a replicase polyprotein",
             "pp1ab replicase polyprotein")
ORF1A_NAMES = ("orf1a", "pp1a", "polyprotein 1a", "replicase 1a", "pp1a polyprotein",
               "pp1a replicase polyprotein")


def main():
    os.chdir(W)
    gs = {l.split('\t')[0]: l.rstrip('\n').split('\t')
          for l in open('dumps/%s/%s.id_name_gs' % (MODULE, T))}
    locs = {}
    for l in open('cls/feature_locs.tsv'):
        p = l.rstrip('\n').split('\t')
        if p[1] in gs: locs[p[0]] = p
    seq = {l.split('\t')[0]: l.rstrip('\n').split('\t')[1] for l in open('Tob.uniq.seq')}
    d = os.path.join(HERE, 'bins'); os.makedirs(d, exist_ok=True)

    #  ---- 1. cluster every structural-block protein by homology
    small = {f: p for f, p in locs.items() if 0 < int(p[6] or 0) < 2500}
    with open('%s/in.faa' % d, 'w') as o:
        for f, p in sorted(small.items()): o.write('>%s\n%s\n' % (f, seq[p[8]]))
    subprocess.run('mmseqs easy-cluster %s/in.faa %s/c %s/tmp --min-seq-id %s -c %s '
                   '--cov-mode 0 -e 1e-3 >/dev/null 2>&1' % (d, d, d, MIN_SEQ_ID, COV),
                   shell=True, check=True)
    cl = collections.defaultdict(list)
    for l in open('%s/c_cluster.tsv' % d):
        a, b = l.split(); cl[a].append(b)
    clusters = [m for _, m in sorted(cl.items(), key=lambda x: -len(x[1]))]

    #  ---- 2. anchors the literature names, from NC_024709
    anch = {}
    for f, p in locs.items():
        if p[1] == ANCHOR_GENOME:
            for k, (s, e) in ANCHORS.items():
                if int(p[3]) == s and int(p[4]) == e: anch[k] = seq[p[8]]
    assert len(anch) == 3, 'anchors not found in %s' % ANCHOR_GENOME
    with open('%s/anchors.faa' % d, 'w') as o:
        for k, s in anch.items(): o.write('>%s\n%s\n' % (k, s))
    with open('%s/q.faa' % d, 'w') as o:
        for i, mem in enumerate(clusters, 1):
            for m in mem: o.write('>C%d|%s\n%s\n' % (i, m, seq[locs[m][8]]))
    subprocess.run('makeblastdb -in %s/anchors.faa -dbtype prot >/dev/null && '
                   'blastp -query %s/q.faa -db %s/anchors.faa -evalue 10 -outfmt '
                   '"6 qseqid sseqid pident length evalue bitscore" -max_target_seqs 5 '
                   '-num_threads 8 -out %s/hits.tsv' % (d, d, d, d), shell=True, check=True)
    best = {}
    for l in open('%s/hits.tsv' % d):
        q, s, pid, ln, ev, bit = l.rstrip('\n').split('\t')
        c = q.split('|')[0]
        if c not in best or float(bit) > best[c][3]:
            best[c] = (s, float(ev), int(ln), float(bit))

    #  ---- 3. assign each cluster
    assign, unc, notmod = {}, [], collections.Counter()
    for i, mem in enumerate(clusters, 1):
        c = 'C%d' % i
        names = {locs[m][7].lower() for m in mem}
        if any(n in REPLICASE for n in names):
            #  Under 2,500 aa and named for the replicase: a fragment of pp1a
            #  or pp1ab, never the protein. pp1a itself is bound by name below.
            assign[c] = None
            notmod[('replicase fragment: under 2,500 aa and deposited under a pp1a / '
                    'pp1ab / ORF1b name, so it is part of a polyprotein, not a protein',
                    '; '.join(sorted(names))[:60])] += len(mem)
            continue
        b = best.get(c)
        if b and b[1] <= EVAL and b[2] >= MINALN:
            assign[c] = b[0]
        else:
            uniq = {seq[locs[m][8]] for m in mem}
            if len(uniq) >= 2: unc.append(c)
            else:
                assign[c] = None
                notmod[('one distinct sequence only (a single isolate, usually deposited '
                        'twice); a profile would describe one genome',
                        '; '.join(sorted(names))[:60])] += len(mem)
    #  Two accessory orthogroups are ONE feature, measured on a panel run:
    #  the 3'-terminal ~500 aa accessory of Pregotovirus and the ORF reported
    #  to resemble paramyxovirus HN (PMID 25063552) called each other's loci at
    #  identical coordinates in Xinzhou nematode virus 6 and Bellinger River
    #  virus. blastp finds no HSP between the two collections even at e=10, so
    #  this is not visible before a profile exists -- and a profile that scores
    #  on the other bin's locus cannot be separated by raising a cutoff. They
    #  are merged, with one profile per orthogroup inside the one feature.
    MERGE = {}
    for a, b in (('C16', 'C5'),):          # HN-like -> the Pregotovirus 3' accessory
        if a in unc and b in unc: MERGE[a] = b
    unc = [c for c in unc if c not in MERGE]
    for n, c in enumerate(unc, 1): assign[c] = 'UNC%d' % n
    for a, b in MERGE.items(): assign[a] = assign[b]

    #  ---- 4. write the collections
    cdir = os.path.join(HERE, 'collections', MODULE)
    os.makedirs(cdir, exist_ok=True)
    for f in os.listdir(cdir):
        if f.endswith('.fasta'): os.remove(os.path.join(cdir, f))
    out = collections.defaultdict(list)
    syn = collections.Counter()
    dropped = collections.Counter()
    for i, mem in enumerate(clusters, 1):
        key = assign.get('C%d' % i)
        for m in mem:
            p = locs[m]; s = seq[p[8]]
            g = gs[p[1]][3] or 'unclassified'
            if any(x in p[7].lower() for x in PARTIAL):
                dropped['partial'] += 1
                syn[(p[7], g, 'DISCARDED')] += 1
                continue
            syn[(p[7], g, key or 'NOT_MODELLED')] += 1
            if key: out[key].append((m, s))
    #  ORF1A is bound by name, not by the clustering above: the clustering only
    #  sees proteins under 2,500 aa, and pp1a is 3,840-6,570. The window is
    #  what separates a real pp1a from the 66 aa and 294 aa "polyprotein 1a" /
    #  "polyprotein 1ab" fragments, which the clustering does see.
    for f, p in sorted(locs.items()):
        n, s = p[7].lower().strip(), seq[p[8]]
        if len(s) < 2500: continue
        g = gs[p[1]][3] or 'unclassified'
        if any(x in n for x in PARTIAL):
            dropped['partial'] += 1; syn[(p[7], g, 'DISCARDED')] += 1; continue
        if n in ORF1A_NAMES and 3800 <= len(s) <= 6700:
            out['ORF1A'].append((f, s)); syn[(p[7], g, 'ORF1A')] += 1
        else:
            syn[(p[7], g, 'NOT_MODELLED')] += 1
            notmod[('pp1ab / ORF1b: the -1 frameshift product and its 3\' half are '
                    'called from transcript_edit references, never from a PSSM', p[7])] += 1
    for k, rows in sorted(out.items()):
        seen = set()
        with open(os.path.join(cdir, k + '.fasta'), 'w') as fh:
            for f, s in rows:
                if s in seen: continue
                seen.add(s); fh.write('>%s\n%s\n' % (f, s))
    #  ---- 5. the tracking files every downstream tool expects
    td = os.path.join(HERE, 'collections')
    with open(os.path.join(td, 'synonyms.tsv'), 'w') as fh:
        fh.write('count\tannotation\tgenus\tbins_to\n')
        for (a, g, k), c in sorted(syn.items(), key=lambda x: -x[1]):
            fh.write('%d\t%s\t%s\t%s\n' % (c, a, g, k))
    with open(os.path.join(td, 'NOT_MODELLED.tsv'), 'w') as fh:
        fh.write('features\tannotation\treason\n')
        for (why, names), c in sorted(notmod.items(), key=lambda x: -x[1]):
            fh.write('%d\t%s\t%s\n' % (c, names, why))
    with open(os.path.join(td, 'UNASSIGNED_TRACKING.tsv'), 'w') as fh:
        fh.write('features\tuniq_seqs\tnot_already_called\tannotation\n')
    print('%-8s %-6s %s' % ('feature', 'seqs', 'length range'))
    for k in sorted(out):
        lens = sorted(len(s) for _f, s in out[k])
        n = len({s for _f, s in out[k]})
        print('%-8s %-6d %d-%d' % (k, n, lens[0], lens[-1]))
    print('dropped:', dict(dropped))
    print('clusters with one distinct sequence, recorded NOT_MODELLED:',
          sum(notmod.values()))


if __name__ == '__main__':
    main()
