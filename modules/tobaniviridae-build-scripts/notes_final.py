#!/usr/bin/env python3
"""Rewrite section 9 of the build notes from the final measurements."""
import os, subprocess, sys
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
run = lambda *a: subprocess.run([sys.executable] + list(a), capture_output=True, text=True).stdout
panel = {M: run('scripts/per_genome_calls.py', M).split(
    'per-feature call rate from the feature tables:')[-1].strip()
    for M in ('Torovirinae', 'Piscanivirinae', 'Remotovirinae', 'Serpentovirinae')}
nc = run('scripts/nc_summary.py')
flat = run('scripts/flat_routing.py').split('\n\n')[0]
cuts = run('scripts/measure_cutoffs.py')
body = """## 9. What the measurements say

### Profiles

| module | features | profiles | single-sequence | collections |
|---|---|---|---|---|
| Torovirinae | 7 (1 `transcript_edit`) | 57 | 42 | ORF1A 32, S 62, M 37, HE 87, N 52, MPRO 32 |
| Piscanivirinae | 6 (1) | 24 | 19 | ORF1A 6, S 13, M 9, N 7, MPRO 5 |
| Remotovirinae | 7 (1) | 10 | 5 | ORF1A 16, S 12, M 8, N 5, G2 10, ORF6 9 |
| Serpentovirinae | 14 (1) | 85 | 83 | S 17, M 17, N 14, ORF1A 12, UNC1-9 2-7 each |

`-m` was lowered to 2, recorded in every `BUILD_PARAMS`, because at the
documented floor of 5 most of these features build **zero** profiles: 6 to 17
unique sequences per protein spread over several genera never reach five per
cluster. `-mi` was **not** lowered — the identity floor is what makes a profile
a model of one protein. The cost of holding it is that most sequences end up in
a cluster of one and get a single-sequence profile, built by duplicating the
sequence into a two-row alignment, which is what psiblast will accept. A
single-sequence profile is a lookup table with a bit score: it calls a close
relative and does not generalise.

### Held-out panels

Drawn from the build's own contigs with the reference genomes excluded by
accession. Rates are computed from the feature tables, by gene symbol — NOT
from `run_gto_eval.py`'s per-feature figure, which counts a feature as called
when its annotation STRING appears in the GTO and so reports one identical rate
for the nine Serpentovirinae features that share "Uncharacterized
lineage-specific protein".

%s
### Bit cutoffs, reset from measured scores

The generator derives a cutoff from the profile, and for a short protein that
lands at 30 bits. On Berne virus a second *Membrane glycoprotein* was called at
position 11,322 — inside the replicase ORF — scoring **30.34 against a cutoff
of 30**, while the genuine M on that genome scored 476.8. Every feature's
cutoff was reset to half the lowest score its own panel calls reached, rounded
down, never below the generator's value, and **capped at 500**: half of
Piscanivirinae ORF1A's lowest panel score is 3,000, and at 3,000 the divergent
fish tobaniviruses whose pp1a matches partially stop being called at all
(measured: 70 per cent falling to 33). Features with fewer than five observed calls — UNC4 to
UNC9 — keep the generator's value, because one observation says nothing about
how a divergent member of the same orthogroup would score.

The measured score distributions behind those cutoffs:

```
%s```

After the change the false call is gone, Berne annotates ORF1a, Mpro, S, M and
N (correctly no HE: its HE is a pseudogene), and every panel rate above is
unchanged.

### Whole-taxon run

Every Tobaniviridae record BV-BRC holds between 1 and 50 kb (353 genomes),
clustered at 95%% nucleotide identity, one exemplar per cluster, routed against
a repository holding this family and nothing else.

**Routing.** 91 of the 117 exemplars were annotated through the GTO pipeline
and **every one of them routed**. The other 26 never reached a module:
`rast-create-genome` returned HTTP 403 from the BV-BRC ID server, which the
coverage report then counts as UNROUTED — an infrastructure failure wearing a
routing failure's name. Re-run through the flat annotator, which needs no ID
server:

```
%s
```

So 115 of 117 distinct genome types route, 98.3%%. The two that do not are a
1,625 nt and a 1,316 nt fragment.

**Quality, with an explicit completeness cutoff.** "Near-complete" is measured
against the genome each module expects, not a family-wide constant: an
exemplar counts as near-complete when it reaches 90%% of the median length of
its own module's reference genomes (Torovirinae 28,511 nt, Piscanivirinae
27,004, Remotovirinae 20,386, Serpentovirinae 30,104). 49 of the 117 clear
that; the other 68 have a median length of 4,731 nt and **90.5%% of them are
correctly flagged incomplete**.

The 49 were re-run against the module as it now ships:

```
%s```

**The six poor ones are all profile coverage, and all visible.** Four carry
`ORF1a polyprotein: Feature is too short` — a divergent genome whose
polyprotein matches partially, which is what a thin profile set does; one is
missing N; and Bangali torovirus carries `Feature is too long` on ORF1a. None
is a silent failure, and the fix for all of them is more sequences rather than
different thresholds.

### Where ORF1AB misses, classified

| cause | n | what it means |
|---|---|---|
| ambiguous bases in the genome | 7 | the annotator demotes to `partial_cds` by design, because the gap-fill draws its bases from the subject; the same ambiguity is why those genomes are not references |
| no pp1a on the contig | 2 | 5'-truncated submissions |
| no reference inside the 95%% gate | 2 | divergent, or a partial best HSP |

`check_te_leak.py` states the other half plainly: the reference sets are built
from every genome that yields a clean construct, so most panel genomes are
also references and ORF1AB's rate on them is self-recall. On the genomes that
are NOT references it is 2 of 13 — and 11 of those 13 are exactly the genomes
rejected from the reference set for ambiguity or partiality. The two
populations coincide, which is the honest reading.
""" % (''.join('**%s** (%d panel genomes)\n\n```%s\n```\n\n' %
               (M, sum(1 for _ in open('panel/%s.tsv' % M)), panel[M])
               for M in ('Torovirinae', 'Piscanivirinae', 'Remotovirinae', 'Serpentovirinae')),
       cuts, flat, nc)
p = 'TOBANIVIRIDAE_BUILD_NOTES.md'
txt = open(p).read()
a = txt.index('## 9. What the measurements say')
b = txt.index('## 10. Open questions')
open(p, 'w').write(txt[:a] + body + '\n---\n\n' + txt[b:])
print('section 9 rewritten (%d chars)' % len(body))
