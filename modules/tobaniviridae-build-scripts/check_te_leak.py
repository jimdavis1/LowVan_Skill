#!/usr/bin/env python3
"""Does the panel overlap the transcript-edit reference set?

A reference set built from every genome includes the panel genomes, and then
ORF1AB's call rate on that panel is self-recall wearing different clothes. This
reports the overlap and the call rate on the genomes that are NOT references.
usage: check_te_leak.py <Module>"""
import glob, os, re, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
M = sys.argv[1]
refs = set()
for l in open('%s/Transcript-Editing/%s/ORF1AB.fasta' % (M, M)):
    if l.startswith('>'):
        m = re.search(r'\|(\d+\.\d+)/', l)
        if m: refs.add(m.group(1))
called, panel = set(), set()
for f in sorted(glob.glob('tbl/%s/*.feature.tbl' % M)):
    g = os.path.basename(f)[:-len('.feature.tbl')]
    panel.add(g)
    if any(l.split('\t')[6] == 'ORF1ab' for l in open(f) if l.count('\t') > 6):
        called.add(g)
held = panel - refs
print('%s: %d transcript-edit references, %d panel genomes, %d of the panel are '
      'also references' % (M, len(refs), len(panel), len(panel & refs)))
print('  ORF1AB called on %d/%d reference genomes in the panel'
      % (len(called & refs), len(panel & refs)))
print('  ORF1AB called on %d/%d genomes that are NOT references%s'
      % (len(called & held), len(held),
         (' -- ' + ', '.join(sorted(held - called))) if held - called else ''))
