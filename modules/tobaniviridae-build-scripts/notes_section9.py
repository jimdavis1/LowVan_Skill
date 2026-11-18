#!/usr/bin/env python3
"""Fill section 9 of the build notes from the measurements on disk."""
import os, subprocess, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
cov = subprocess.run([sys.executable, 'scripts/coverage_summary.py'],
                     capture_output=True, text=True).stdout
panel = {}
for M in ('Torovirinae', 'Piscanivirinae', 'Remotovirinae', 'Serpentovirinae'):
    panel[M] = subprocess.run([sys.executable, 'scripts/per_genome_calls.py', M],
                              capture_output=True, text=True).stdout.split(
        'per-feature call rate from the feature tables:')[-1].strip()
body = """## 9. What the measurements say

### Profiles

| module | features | profiles | single-sequence | collections |
|---|---|---|---|---|
| Torovirinae | 7 (1 `transcript_edit`) | 57 | 42 | ORF1A 32, S 62, M 37, HE 87, N 52, MPRO 32 |
| Piscanivirinae | 6 (1) | 24 | 19 | ORF1A 6, S 13, M 9, N 7, MPRO 5 |
| Remotovirinae | 7 (1) | 10 | 5 | ORF1A 16, S 12, M 8, N 5, G2 10, ORF6 9 |
| Serpentovirinae | 14 (1) | 85 | 83 | S 17, M 17, N 14, ORF1A 12, UNC1-9 2-7 each |

`-m` was lowered to 2 (recorded in every `BUILD_PARAMS`) because at the
documented floor of 5 most of these features build **zero** profiles: 6 to 17
unique sequences per protein spread over several genera never reach five per
cluster. `-mi` was **not** lowered — the identity floor is what makes a profile
a model of one protein. The cost of holding it is that most sequences end up
in a cluster of one, so they get a single-sequence profile built from the
sequence duplicated into a two-row alignment, which psiblast will accept. A
single-sequence profile is a lookup table with a bit score: it calls a close
relative and does not generalise.

### Held-out panels

Drawn from the build's own contigs with the reference genomes excluded by
accession. Rates are computed from the feature tables, by gene symbol — NOT
from `run_gto_eval.py`'s own per-feature figure, which counts a feature as
called when its annotation STRING appears in the GTO and therefore reports one
identical rate for the nine Serpentovirinae features that share
"Uncharacterized lineage-specific protein".

%s

0 duplicate calls, and `check_cleaved_ends` / the junction chemistry in §7
confirm the Mpro boundaries.

### Whole-taxon run

Every Tobaniviridae record BV-BRC holds between 1 and 50 kb (353 genomes),
clustered at 95%% nucleotide identity, one exemplar per cluster, routed against
a repository holding this family and nothing else.

```
%s```

Quality here is a property of how the family is deposited, not of the modules:
the median poor-quality exemplar is 6,550 nt against genomes of 20-36 kb, and
every flag on it is "missing essential feature" — ORF1a on 30, N on 26, M on
24, S on 22. Those are single-gene partial submissions.

The Serpentovirinae column is worth reading twice. Most of its near-complete
exemplars are the metagenomic `sp.` contigs dropped at intake, so they are in
no collection and are not transcript-edit references — `ORF1AB` falls to 30%%
there because the correction gate needs a reference within 95%% identity, while
the structural proteins are still called on 94-100%% of them. Those profiles
are reaching genomes that contributed nothing to them.

### Where ORF1AB misses, classified

Each panel miss was traced rather than counted:

| cause | n | what it means |
|---|---|---|
| ambiguous bases in the genome | 7 | the annotator demotes to `partial_cds` by design, because the gap-fill draws its bases from the subject; the same ambiguity is why those genomes are not references |
| no pp1a on the contig | 2 | 5'-truncated submissions |
| no reference inside the 95%% gate | 2 | divergent, or a partial best-HSP |

`check_te_leak.py` states the other half of this plainly: the reference sets
are built from every genome that yields a clean construct, so most panel
genomes are also references and ORF1AB's rate on them is self-recall. On the
genomes that are NOT references it is 2 of 13 — and 11 of those 13 are
precisely the genomes rejected from the reference set for ambiguity or
partiality. The two populations coincide, which is the honest reading.
""" % (''.join('**%s** (%d genomes)\n\n```%s\n```\n\n' %
               (M, sum(1 for _ in open('panel/%s.tsv' % M)), panel[M])
               for M in ('Torovirinae', 'Piscanivirinae', 'Remotovirinae', 'Serpentovirinae')),
       cov)
p = 'TOBANIVIRIDAE_BUILD_NOTES.md'
txt = open(p).read()
a = txt.index('## 9. What the measurements say')
b = txt.index('## 10. Open questions')
open(p, 'w').write(txt[:a] + body + '\n---\n\n' + txt[b:])
print('section 9 written, %d chars' % len(body))
