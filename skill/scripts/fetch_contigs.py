#!/usr/bin/env python3
"""Pull per-genome contigs for the near-complete members of a taxon.

Rep-contig selection and `evaluate_coverage` both need nucleotide contigs, one
FASTA per genome, which the protein dump does not carry. This fetches them from
the BV-BRC API and writes the metadata table the downstream scripts expect.

    python3 fetch_contigs.py --genus Sapovirus --min-len 6500 \
            --out ~/lowvan_scratch/calici/Contigs.Sapovirus

`--min-len` is a floor on genome length, not on contig length, and it is the
one judgement here. Set it just under the shortest complete genome in the
taxon: too high silently drops divergent members with shorter genomes, too low
fills the directory with genotyping fragments that then dominate the
rep-contig clustering. A segmented taxon wants the length of its smallest
segment and `--rank family` or `--rank species` as appropriate.

Writes `<out>/<genome_id>.fna` per genome and `<out>.metadata`, a three-column
TSV that `build_rep_contigs.py --meta` and `repcontig_budget.py` both read.
"""
import argparse, collections, json, os, sys

try:
    import requests
except ImportError:
    sys.exit("needs `requests` -- pip install requests")

API = "https://www.bv-brc.org/api"
HDR = {"Content-Type": "application/rqlquery+x-www-form-urlencoded",
       "Accept": "application/json"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genus", help="taxon value to match")
    ap.add_argument("--rank", default="genus",
                    help="BV-BRC field to match on (genus, family, species)")
    ap.add_argument("--min-len", type=int, default=6500,
                    help="genome length floor in nt (default 6500)")
    ap.add_argument("--out", required=True, help="contig directory to write")
    ap.add_argument("--batch", type=int, default=100,
                    help="genomes per sequence query (default 100)")
    a = ap.parse_args()
    if not a.genus:
        sys.exit("--genus is required")

    out = os.path.expanduser(a.out)
    os.makedirs(out, exist_ok=True)

    #  Page the genome list. The API caps a page at 25,000 rows and a taxon
    #  like Norovirus exceeds that, so this loops rather than asking once.
    rows, got = [], 0
    while True:
        q = ("eq(%s,%s)&gt(genome_length,%d)"
             "&select(genome_id,genome_name,genome_length)"
             "&limit(25000,%d)&sort(+genome_id)" % (a.rank, a.genus, a.min_len, got))
        r = requests.post(API + "/genome/", data=q, headers=HDR, timeout=600).json()
        r = [x for x in r if isinstance(x, dict)]
        if not r:
            break
        rows += r
        got += len(r)
        if len(r) < 25000:
            break
    print("  %d %s genomes >= %d nt" % (len(rows), a.genus, a.min_len), flush=True)
    if not rows:
        sys.exit("  nothing to fetch -- check --rank and --min-len")

    ids = [x["genome_id"] for x in rows]
    name = {x["genome_id"]: x.get("genome_name", "") for x in rows}

    meta = open(out + ".metadata", "w")
    meta.write("genome.genome_id\tgenome.genome_name\tgenome.genome_length\n")
    n = 0
    for i in range(0, len(ids), a.batch):
        q = ("in(genome_id,(%s))&select(genome_id,sequence_id,sequence)&limit(5000)"
             % ",".join(ids[i:i + a.batch]))
        d = requests.post(API + "/genome_sequence/", data=q, headers=HDR,
                          timeout=600).json()
        by = collections.defaultdict(list)
        for x in d:
            if isinstance(x, dict) and x.get("sequence"):
                by[x["genome_id"]].append(x)
        for g, v in by.items():
            with open(os.path.join(out, g + ".fna"), "w") as fh:
                for j, x in enumerate(v, 1):
                    fh.write(">%s\n" % (x.get("sequence_id") or "%s_%d" % (g, j)))
                    s = x["sequence"]
                    for k in range(0, len(s), 70):
                        fh.write(s[k:k + 70] + "\n")
            meta.write("%s\t%s\t%d\n"
                       % (g, name.get(g, ""), sum(len(x["sequence"]) for x in v)))
            n += 1
        if i and i % 5000 == 0:
            print("    %d/%d" % (i, len(ids)), flush=True)
    meta.close()
    print("  wrote %d genome contig files -> %s" % (n, out))
    print("  metadata -> %s.metadata" % out)
    if n < len(ids):
        print("  note: %d genome(s) returned no sequence rows" % (len(ids) - n))


main()
