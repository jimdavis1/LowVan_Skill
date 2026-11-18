#!/usr/bin/env python3
"""Why did the transcript-edit call fall back to 'Uncorrected ... region'?

The caller blastn's each nucleotide reference against the genome with
-word_size 28 and requires >=95% identity AND >=95% query coverage in ONE HSP
pair before it fills the frameshift gap; between 80% and 95% it emits an
uncorrected partial_cds instead. This prints, for one genome, the best
reference and what it scored on each gate.
usage: te_gate_probe.py <Module> <genome_id>"""
import os, re, subprocess, sys, tempfile
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
M, G = sys.argv[1], sys.argv[2]
ref = '%s/Transcript-Editing/%s/ORF1AB.fasta' % (M, M)
gen = 'Contigs/%s.fna' % G
out = subprocess.run(
    'blastn -query %s -subject %s -word_size 28 -outfmt '
    '"6 qseqid qlen qstart qend sstart send pident length gaps bitscore"'
    % (ref, gen), shell=True, capture_output=True, text=True).stdout
rows = [l.split('\t') for l in out.splitlines()]
if not rows:
    print('%s %s: no blastn hit at all to any reference' % (M, G)); sys.exit()
best = {}
for r in rows:
    q = r[0]
    if q not in best or float(r[9]) > float(best[q][9]): best[q] = r
top = sorted(best.values(), key=lambda r: -float(r[9]))[:4]
print('%s %s: %d reference(s) hit; best four by bitscore' % (M, G, len(best)))
for r in top:
    qlen, qs, qe, pid, ln, gaps = int(r[1]), int(r[2]), int(r[3]), float(r[6]), int(r[7]), int(r[8])
    cov = 100.0 * (qe - qs + 1) / qlen
    print('   %-46s %5.1f%% identity, best HSP covers %5.1f%% of the reference '
          '(%d/%d nt), %d gaps' % (r[0].split('|')[0][:46], pid, cov, qe - qs + 1, qlen, gaps))
print('   gate: >=95%% identity AND >=95%% coverage in one HSP -> corrected call;'
      ' 80-95%% -> "Uncorrected ... region"')
