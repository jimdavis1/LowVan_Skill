#!/usr/bin/env python3
"""Renumber the single-sequence profiles from <Module>.<KEY>.ssN to plain
numbers, continuing from the feature's highest existing profile id.

apply_copy_num.py resolves a feature from the GTO's family_assignments by
stripping a TRAILING NUMERIC profile id -- "Bocaparvovirus.NP1.3_1" -> NP1. An
"ss3" suffix is not numeric, so the key stayed "M.ss3" and every profile was
counted as a feature of its own: Torovirinae M measured 33.3% from one profile
while the feature is called on 100% of the panel. Nothing errors; the rates are
simply wrong, and a wrong rate here adds or drops copy_num.
usage: renumber_singletons.py <Module> [--write]"""
import glob, os, re, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
M = sys.argv[1]; WRITE = '--write' in sys.argv
base = '%s/Alignments/%s' % (M, M)
for feat in sorted(os.listdir(base)):
    pd = os.path.join(base, feat, 'pssms')
    if not os.path.isdir(pd): continue
    nums = []
    for p in glob.glob('%s/*.pssm' % pd):
        m = re.match(r'.*\.([0-9]+(?:_[0-9]+)?)\.pssm$', os.path.basename(p))
        if m: nums.append(int(m.group(1).split('_')[0]))
    nxt = max(nums) + 1 if nums else 1
    for p in sorted(glob.glob('%s/*.ss[0-9]*.pssm' % pd)):
        b = os.path.basename(p)
        tag = re.search(r'\.(ss[0-9]+)\.pssm$', b).group(1)
        new = b.replace('.%s.pssm' % tag, '.%d.pssm' % nxt)
        ali = os.path.join(base, feat, 'reclustered_alis', '%s.fa' % tag)
        nali = os.path.join(base, feat, 'reclustered_alis', '%d.fa' % nxt)
        print('   %s -> %s' % (b, new))
        if WRITE:
            os.rename(p, os.path.join(pd, new))
            if os.path.exists(ali): os.rename(ali, nali)
        nxt += 1
print('%s: %s' % (M, 'renamed' if WRITE else 'dry run'))
