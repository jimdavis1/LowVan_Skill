#!/usr/bin/env python3
"""One profile per sequence, for a feature whose members do not cluster.

Serpentovirinae and Piscanivirinae hold two to seventeen unique sequences per
protein spread over several genera, 25-60% identical to each other. At the
documented -mi 0.8 every sequence is its own cluster, so the main pass and the
leftover pass both return nothing and the feature ships with no profile at all.
The identity floor is NOT lowered to force clustering -- that would redefine
"the same protein" -- so each surviving sequence gets its own profile instead.

psiblast refuses an MSA of one row, so the row is duplicated: the profile is
then the sequence's own residue frequencies, which is exactly what a
single-sequence profile is. A cutoff for it cannot be read off a score
distribution, so it must be MEASURED against a foreign sequence -- that is
measure_cutoffs.py, which runs after this.

usage: build_singleton_pssms.py <Module> <KEY> [--write]
"""
import glob, os, re, shutil, subprocess, sys, tempfile
W = os.path.expanduser('~/lowvan_scratch/tobani')
KIT = os.path.join(W, 'kit/skill/scripts')


def read_fasta(p):
    h, s = None, []
    for line in open(p, errors='replace'):
        if line.startswith('>'):
            if h is not None: yield h, ''.join(s)
            h, s = line[1:].strip(), []
        else: s.append(line.strip())
    if h is not None: yield h, ''.join(s)


def build(seq, title, out, tmp):
    """psiblast over a two-row MSA holding the same sequence twice."""
    msa = os.path.join(tmp, 'g.msa')
    with open(msa, 'w') as fh:
        fh.write('>1 %s\n%s\n>2\n%s\n' % (title, seq, seq))
    subj = os.path.join(tmp, 'g.subj')
    with open(subj, 'w') as fh:
        fh.write('>1 %s\n%s\n' % (title, seq))
    raw = os.path.join(tmp, 'g.raw')
    if os.path.exists(raw): os.remove(raw)
    subprocess.run(['psiblast', '-subject', subj, '-in_msa', msa, '-out_pssm', raw],
                   capture_output=True)
    if not os.path.exists(raw): return False
    txt = subprocess.run([sys.executable, os.path.join(KIT, 'norm_pssm.py'), raw],
                         capture_output=True, text=True).stdout
    os.remove(raw)
    open(out, 'w').write(txt)
    return True


def main():
    M, KEY = sys.argv[1], sys.argv[2]
    write = '--write' in sys.argv
    os.chdir(W)
    import json
    anno = json.load(open('%s/features.json' % M))[KEY]['anno']
    fd = '%s/Alignments/%s/%s' % (M, M, KEY)
    coll = '%s/collections/%s/%s.fasta' % (M, M, KEY)
    rows = list(read_fasta(coll))
    #  sequences a profile already covers, by exact match to an existing
    #  alignment's ungapped rows
    covered = set()
    for ali in glob.glob('%s/corrected_alis/*.fa' % fd) + glob.glob('%s/reclustered_alis/*.fa' % fd):
        for _h, s in read_fasta(ali): covered.add(s.replace('-', '').upper())
    todo = [(h, s) for h, s in rows if s.upper() not in covered]
    print('%s %s: %d sequences, %d already covered by a profile, %d to build'
          % (M, KEY, len(rows), len(rows) - len(todo), len(todo)))
    if not write:
        print('  (dry run; rerun with --write)'); return
    pd, rd = os.path.join(fd, 'pssms'), os.path.join(fd, 'reclustered_alis')
    os.makedirs(pd, exist_ok=True); os.makedirs(rd, exist_ok=True)
    nums = []
    for q in glob.glob('%s/*.pssm' % pd):
        m = re.match(r'.*\.([0-9]+)(?:_[0-9]+)?\.pssm$', os.path.basename(q))
        if m: nums.append(int(m.group(1)))
    nstart = max(nums) if nums else 0
    tmp = tempfile.mkdtemp(prefix='ss.')
    built = 0
    try:
        for i, (h, s) in enumerate(todo, 1):
            tag = '%d' % (nstart + i)
            out = os.path.join(pd, '%s.%s.%s.pssm' % (M, KEY, tag))
            if build(s.replace('-', '').upper(), anno, out, tmp):
                with open(os.path.join(rd, '%s.fa' % tag), 'w') as fh:
                    fh.write('>%s\n%s\n>%s_copy\n%s\n' % (h, s, h.split()[0], s))
                built += 1
            else:
                print('  psiblast failed on', h.split()[0])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print('  built %d single-sequence profile(s); %d profiles total for %s'
          % (built, len(glob.glob('%s/*.pssm' % pd)), KEY))


if __name__ == '__main__':
    main()
