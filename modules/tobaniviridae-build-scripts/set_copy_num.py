#!/usr/bin/env python3
"""Apply copy_num, filtered by the project's standing rules.

apply_copy_num.py measured these rates on each module's held-out panel; this
script records which of its ADDs are taken and why the rest are not, because
copy_num is what marks a feature essential and a wrong one flags every genome
in the taxon.

Taken: a core protein whose measured rate is >= 0.85 and which has no
copy-number or cutoff struggle.

Refused, with the reason:
  HE (Torovirinae)      95.8% measured, but the equine torovirus HE gene is a
                        deleted pseudogene (PMID 29950409) and cell adaptation
                        truncates it in bovine torovirus (PMID 34817201), so
                        copy_num would flag genomes whose HE is genuinely gone.
  ORF1AB (all four)     85.7-87.5% measured. Every miss on the panels is a
                        genome carrying ambiguous bases, where the annotator
                        demotes the call to partial_cds by design; marking it
                        essential would flag exactly those genomes twice.
  MPRO (two modules)    100% raw but 15-20% standalone -- it is a mature
                        peptide lying inside ORF1A, and apply_copy_num blocks
                        that case rather than granting essentiality on another
                        feature's locus.
  ORF1A (Piscanivirinae) 84.6%, below the threshold: two panel contigs are
                        5'-truncated and carry no pp1a.
  G2, ORF6, UNC1-9      uncharacterized or lineage-specific, and in several
                        cases genus-specific accessories. Never marked
                        essential, whatever the rate.
"""
import json, os
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
TAKE = {
    'Torovirinae':     ['M', 'N', 'S', 'ORF1A'],
    'Piscanivirinae':  ['M', 'N', 'S'],
    'Remotovirinae':   ['M', 'N', 'S', 'ORF1A'],
    'Serpentovirinae': ['M', 'N', 'S', 'ORF1A'],
}
for M, keys in TAKE.items():
    p = '%s/%s_Viral_PSSM.json' % (M, M)
    d = json.load(open(p))
    feats = d[M]['features']
    for k, v in feats.items():
        if k in keys: v['copy_num'] = 1
        else: v.pop('copy_num', None)
    json.dump(d, open(p, 'w'), indent=3, sort_keys=True)
    print('%-16s copy_num on %s' % (M, ', '.join(sorted(keys))))
