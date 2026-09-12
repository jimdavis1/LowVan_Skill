# The BV-BRC input dump

Every module starts here. Five tab-separated files, named for the taxon (`<T>`).

## The contract

```
<T>.id_md5        feature_id \t md5
<T>.id_name_gs    genome_id \t genome_name \t family \t genus
<T>.uniq.md5      md5
<T>.uniq.id_ann   representative_feature_id \t annotation
<T>.uniq.seq      md5 \t sequence
```

The three `uniq.*` files are **line-index aligned**: line *n* of `uniq.md5`,
line *n* of `uniq.id_ann` and line *n* of `uniq.seq` describe the same unique
sequence. That index is the join key. Sorting any one of them independently
destroys the dump. Check first:

```bash
wc -l <T>.uniq.md5 <T>.uniq.id_ann <T>.uniq.seq     # must be equal
```

`id_md5` is the many-to-one map from every feature in every genome to its unique
sequence. Use it when you need per-genome context — gene order, copy number,
what else that genome encodes. `uniq.id_ann` gives one *representative* feature
id per unique sequence, which is fine for binning and useless for counting: 9
genomes carrying an identical protein appear once.

Counting rule: **count features via `id_md5`, bin sequences via `uniq.*`.**

## Producing the dump

These five files come out of SOLR via `query_PATRIC_bob.pl`
(`assets/query_PATRIC_bob.pl`), which takes a core, an input field and a list of
return fields, and reads identifiers on stdin. The fields are the SOLR schema
fields documented at <https://github.com/BV-BRC/BV-BRC-Solr>.

```bash
# 1. every genome in the family
echo 'Rhabdoviridae' \
  | query_PATRIC_bob.pl -c genome -i family \
      -r "genome_id genome_name family genus"            > Rhabdo.id_name_gs

# 2. every feature in those genomes -> id + md5
cut -f1 Rhabdo.id_name_gs \
  | query_PATRIC_bob.pl -c genome_feature -i genome_id \
      -r "patric_id aa_sequence_md5" | grep 'fig|'       > Rhabdo.id_md5

# 3. the distinct protein sequences, as md5s
grep "CDS\|mat\|peg" Rhabdo.id_md5 | cut -f2 | sort -u \
  | perl -e 'while (<>){chomp; print "$_\n" if /\w/;}'   > Rhabdo.uniq.md5

# 4. one representative feature id + product per distinct sequence
perl -e 'open(IN,"<Rhabdo.id_md5"); my %ids = map{chomp; ($id,$md5)=split /\t/; $md5,$id}(<IN>);
         close IN; while(<>){chomp; print "$ids{$_}\n" if exists $ids{$_};}' < Rhabdo.uniq.md5 \
  | query_PATRIC_bob.pl -c genome_feature -i patric_id \
      -r "patric_id product"                             > Rhabdo.uniq.id_ann

# 5. the sequences themselves
query_PATRIC_bob.pl -c feature_sequence -i md5 \
      -r "md5 sequence" < Rhabdo.uniq.md5                > Rhabdo.uniq.seq
```

Steps 4 and 5 are both driven from `Rhabdo.uniq.md5`, which is what keeps the
three `uniq.*` files line-aligned — the alignment is a consequence of the query
returning rows in input order, not something enforced afterwards. Verify it
rather than trust it.

### The step that goes wrong

**Step 3 is the dangerous one.** It selects protein-coding features with a
positive filter on the feature id, and a positive filter silently omits anything
you did not think to name. `grep "CDS\|mat"` — without `peg` — drops every
BV-BRC-called gene. On Rhabdoviridae that was **15,069 distinct sequences, 45%
of the protein diversity**, and it produced two modules that looked
declared-but-empty for reasons that had nothing to do with the data. Nothing
downstream errors; the proteins are simply invisible.

Prefer excluding what you know is not protein-coding over enumerating what is:

```bash
cut -f2 Rhabdo.id_md5 | sort -u | grep '\w'              > Rhabdo.uniq.md5
```

Rows for RNA, `misc_feature`, `gap` and the rest carry a blank md5 and drop out
on their own, so the filter needs no list of feature types to keep current.

## Check it before building anything

```bash
python3 scripts/check_dump.py --workdir .
```

Two properties, neither of which anything else in the pipeline verifies:

- **completeness** — every distinct md5 named in `id_md5` has a sequence in
  `uniq.seq`. When some do not, it reports which feature-id shapes are affected,
  which names the filter that dropped them.
- **index alignment** — row *n* of `uniq.md5`, `uniq.id_ann` and `uniq.seq`
  describe the same sequence. They are joined by line number, so a query that
  ever paginates or reorders breaks every downstream join at once, silently.

Run it on arrival and again after any re-export.

## What to check on arrival

```bash
# how many features, genomes, unique sequences
wc -l <T>.id_md5 <T>.id_name_gs <T>.uniq.md5

# genus distribution -- this drives module partitioning
cut -f4 <T>.id_name_gs | sort | uniq -c | sort -rn

# annotation string vocabulary -- this drives the triage rules
cut -f2 <T>.uniq.id_ann | sort | uniq -c | sort -rn | head -60
```

Three things to notice immediately:

**Genome-id parsing.** Feature ids look like `fig|11292.22889.CDS.1` but also
`fig|11292.22889.peg.2`, and occasionally other suffixes. Extract the genome id
with `^fig\|(\d+\.\d+)` — anchoring on `.CDS.` alone silently drops every `peg`
feature. In this project that bug affected 36,389 features before it was found.

**Taxonomic skew.** Viral dumps are wildly unbalanced. Rhabdoviridae was 34,777
of 48,074 genomes Lyssavirus. A family-wide percentage is dominated by one
genus, so report coverage per module and per feature, never only family-wide.

**Empty and zero-length records.** Some md5s in `id_md5` have no entry in
`uniq.md5`; some sequences are zero length. Guard every lookup rather than
indexing directly, or the first join crashes.

## Sequences BV-BRC does not have

The LowVan manuscript is explicit that many proteins encoded by these viruses
are absent from BV-BRC and appear only in the literature and the NCBI NR. Where
a protein has fewer than three representatives, exemplars were taken from BLAST
searches against NR. That is a legitimate and expected part of building a
module — record where each added sequence came from, and cite the paper in the
feature's `PMID` field.
