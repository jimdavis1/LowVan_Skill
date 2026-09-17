#!/usr/bin/env python3
"""Generate the Matonaviridae block of Viral_PSSM.json.

Seven features: two polyprotein precursors and the five mature peptides that
tile them exactly.  Every number here was measured; the provenance is below and
in BUILD_NOTES.md.

  bit_cutoff       measure_cutoffs.py -- above the noise ceiling, below the
                   weakest true signal, both measured with -comp_based_stats 0
  min_len/max_len  the collection's own length distribution, widened ~5%
  up/downstream    the literature's cleavage chemistry; see EXT_WHY

Extension flags follow from what each terminus IS.  A protease or signalase cut
site is marked by no codon, so the annotator must not scan outward for one (0);
a genuine ORF start or stop must be found (1).  Getting this backwards lets a
terminus drift with whatever tblastn happened to align -- the failure
check_cleaved_ends.py exists to catch.
"""
import json, os, glob, re, sys

WORK = os.path.dirname(os.path.abspath(__file__))
MODULE = 'Matonaviridae'
OUT = os.path.join(WORK, f'{MODULE}_Viral_PSSM.json')

# --------------------------------------------------------------------- citations
# EVERY citation in this module was found and read by a language model, not
# chosen by a curator, so every one is flagged in PMID_claude_generated as the
# skill requires.  Delete an entry from this registry once a human has read the
# paper and confirmed it supports the claim; the flag then disappears on the
# next regeneration.  That makes clearing the flag an explicit act.
CLAUDE_PMIDS = {
    '9557742':  'Liu, Ropp, Jackson & Frey, J Virol 1998;72:4463. RUB NS protease '
                'cleaves the NSP-ORF product at a SINGLE site giving P150 + P90.',
    '10799588': 'Liang & Gillam, J Virol 2000;74:5133. Scissile bond by mutagenesis '
                '(G1301S abolishes cleavage); p150+p90 needed for +strand synthesis.',
    '10823845': 'Liang, Yao & Gillam, J Virol 2000;74:5412. "two mature products '
                '(p150 and p90) in the order NH2-p150-p90-COOH"; protease in the '
                'C-terminal region of p150; X domain at residues 834-940.',
    '2273395':  'Oker-Blom, Jarvis & Summers, J Gen Virol 1990;71:3047. p110 is '
                'processed to C + E2 + E1; signal peptidase separates the E2 signal '
                'sequence and cleaves E1 from E2.',
    '2214022':  'Suomalainen, Garoff & Baron, J Virol 1990;64:5500. The E2 signal '
                'sequence REMAINS part of the capsid protein in the mature virion '
                '(antipeptide antiserum) and confers membrane association.',
    '11160697': 'Law, Duncan, Esmaili, Nakhasi & Hobman, J Virol 2001;75:1978. E2 '
                'signal peptide stays on the capsid C-terminus after signalase and '
                'is required for perinuclear localisation and assembly.',
    '39207105': 'Das & Kielian, mBio 2024;15:e0196524. 278-PPAY-281 late domain at '
                'the C-terminus of E2; supports the close relationship of human and '
                'animal rubiviruses.',
}

EXT_WHY = {
    'P200': 'CDS: its own initiating Met and its own terminator, so scan for both.',
    'P150': 'N-terminus IS p200 initiating Met -> 1. C-terminus is the protease cut '
            'SRGG/GTCA (PMID 10799588), which no codon marks -> 0.',
    'P90':  'N-terminus is that same protease cut -> 0. C-terminus is p200 own '
            'terminator -> 1.',
    'P110': 'CDS on the subgenomic RNA: own Met, own terminator.',
    'C':    'N-terminus IS p110 initiating Met -> 1. C-terminus is the signalase cut '
            'before E2 -> 0. The E2 signal sequence is retained here (PMID 2214022), '
            'so C runs to the cut with nothing removed.',
    'E2':   'Both termini are signalase cuts (PMID 2273395) -> 0/0.',
    'E1':   'N-terminus is the signalase cut after E2 -> 0. C-terminus is p110 own '
            'terminator -> 1.',
}

#      key      anno                               sym     type          bit  up dn  min   max
FEATS = [
    ('P200', 'Nonstructural polyprotein',          'p200', 'CDS',         400, 1, 1, 1810, 2225,
     ['9557742', '10799588', '10823845']),
    ('P150', 'Protease p150 protein',              'p150', 'mat_peptide', 250, 1, 0, 1020, 1370,
     ['9557742', '10799588', '10823845']),
    ('P90',  'RNA-dependent RNA polymerase',       'p90',  'mat_peptide', 200, 0, 1,  770,  875,
     ['10799588', '10823845']),
    ('P110', 'Structural polyprotein',             'p110', 'CDS',         250, 1, 1,  995, 1205,
     ['2273395']),
    ('C',    'Nucleocapsid protein',               'C',    'mat_peptide', 100, 1, 0,  280,  350,
     ['2214022', '11160697']),
    ('E2',   'Mature envelope glycoprotein E2',    'E2',   'mat_peptide',  90, 0, 0,  245,  345,
     ['2273395', '39207105']),
    ('E1',   'Mature envelope glycoprotein E1',    'E1',   'mat_peptide', 150, 0, 1,  445,  515,
     ['2273395']),
]

# Rep contigs, in the order the .dna files are numbered.
CLOSE = [
    ('Matonaviridae.1.dna', 'NC_001545.2',
     'Rubella virus RVi/1a reference (clade 1)'),
    ('Matonaviridae.2.dna', 'NC_076948.1',
     'Rubella virus RVi/Bismarck.ND.USA/23.08/2B (clade 2)'),
    ('Matonaviridae.3.dna', 'NC_076451.1',
     'Rustrela virus Donkey/19_041-1/2019/Germany'),
    ('Matonaviridae.4.dna', 'NC_076450.1',
     'Ruhugu virus Cyclops leaf-nosed bat/2017/Uganda'),
]

SEGMENT = 'Single RNA Segment'


def main():
    # sanity: every declared feature must have at least one PSSM on disk
    missing = []
    for key, *_ in FEATS:
        n = len(glob.glob(f'{WORK}/Alignments/{MODULE}/{key}/pssms/*.pssm'))
        if n == 0:
            missing.append(key)
    if missing:
        sys.exit(f"refusing to write: no PSSMs for {', '.join(missing)}")

    for fn, _acc, _nm in CLOSE:
        if not os.path.exists(f'{WORK}/Rep-Contigs/{fn}'):
            sys.exit(f"refusing to write: declared rep contig {fn} is not on disk")

    feats = {}
    for key, anno, sym, ftype, bit, up, dn, lo, hi, pmids in FEATS:
        d = {
            'anno': anno,
            'gene_symbol': sym,
            'feature_type': ftype,
            'segment': SEGMENT,
            'bit_cutoff': bit,
            'coverage_cutoff': 0.65,
            'upstream_ext': up,
            'downstream_ext': dn,
            'copy_num': 1,
            'kmers': 1,
            'min_len': lo,
            'max_len': hi,
        }
        # The two lists are DISJOINT, not nested: PMID holds citations a curator
        # has read, PMID_claude_generated holds proposals awaiting a read.  That
        # is the convention Togaviridae shipped and the one check_dlits.py
        # enforces; references/json-schema.md still documents the older nested
        # form and is wrong.  Every Matonaviridae citation is currently a
        # proposal, so PMID is absent on every feature until someone reads them.
        ids = [str(p) for p in pmids]
        curated = [p for p in ids if p not in CLAUDE_PMIDS]
        proposed = [p for p in ids if p in CLAUDE_PMIDS]
        if curated:
            d['PMID'] = curated
        if proposed:
            d['PMID_claude_generated'] = proposed
        feats[key] = d

    block = {
        'close_genomes': {fn: {'genome_ids': acc, 'genome_name': nm}
                          for fn, acc, nm in CLOSE},
        # Rubivirus genomes measured here: min 9429, median 9761, max 9778 nt
        # over 217 near-complete records.  Widened so the quality tool does not
        # flag a healthy genome with ragged ends.
        'segments': {SEGMENT: {'min_len': 9000, 'max_len': 10500,
                               'replicon_geometry': 'linear'}},
        'features': feats,
    }
    with open(OUT, 'w') as fh:
        json.dump({MODULE: block}, fh, indent=3, separators=(',', ' : '),
                  sort_keys=True)
        fh.write('\n')

    npssm = sum(len(glob.glob(f'{WORK}/Alignments/{MODULE}/{k}/pssms/*.pssm'))
                for k, *_ in FEATS)
    print(f"wrote {OUT}")
    print(f"  {len(feats)} features, {npssm} PSSMs, {len(CLOSE)} rep contigs")
    flagged = {p for f in feats.values() for p in f.get('PMID_claude_generated', [])}
    print(f"  {len(flagged)} distinct citations, ALL flagged as model-proposed")
    for p in sorted(flagged, key=int):
        print(f"    PMID {p}: {CLAUDE_PMIDS[p][:72]}")


if __name__ == '__main__':
    main()
