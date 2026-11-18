#!/usr/bin/env python3
"""min_len/max_len for the ORF1AB transcript_edit feature, measured from the
references this build actually wrote (amino acids, 10% slack either side)."""
import glob, json, os
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
for M in ('Torovirinae', 'Piscanivirinae', 'Remotovirinae', 'Serpentovirinae'):
    lens = []
    cur = 0
    for l in open('%s/Transcript-Editing/%s/ORF1AB.fasta' % (M, M)):
        if l.startswith('>'):
            if cur: lens.append(cur // 3 - 1)
            cur = 0
        else: cur += len(l.strip())
    if cur: lens.append(cur // 3 - 1)
    lo, hi = min(lens), max(lens)
    f = json.load(open('%s/features.json' % M))
    f['ORF1AB']['min_len'] = int(lo * 0.9)
    f['ORF1AB']['max_len'] = int(hi * 1.1)
    f['ORF1AB']['build_note'] = ('min_len/max_len in amino acids, from the %d '
        'transcript_edit references this build wrote: observed %d-%d aa, with '
        '10%% slack either side.' % (len(lens), lo, hi))
    json.dump(f, open('%s/features.json' % M, 'w'), indent=1, sort_keys=True)
    print('%-16s %d refs, pp1ab %d-%d aa -> min_len %d max_len %d'
          % (M, len(lens), lo, hi, f['ORF1AB']['min_len'], f['ORF1AB']['max_len']))
