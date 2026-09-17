#!/usr/bin/env python3
"""Download every Matonaviridae contig in the length window from BV-BRC.

Replaces the kit's `New-annotate-viral-taxon.pl -download-only`, which shells
out to BVBRC_clean_contigs.pl -- a program the kit does not ship -- and so
downloads nothing. This queries the genome_sequence core directly, in batches,
and writes the same two artifacts the rest of the workflow reads:

    Contigs/<genome_id>.contigs      one FASTA per genome
    Matonaviridae.metadata             file \t genome_name \t taxon_id
"""
import os, re, subprocess, sys, collections

MIN, MAX = 1000, 50000
meta = {}
for line in open("Matona.genome_meta.tsv"):
    p = line.rstrip("\n").split("\t")
    if len(p) < 7:
        continue
    L = int(p[3] or 0)
    if MIN < L < MAX:
        meta[p[0]] = (p[2], p[1])          # name, taxon_id

ids = sorted(meta)
print(f"{len(ids)} genomes in [{MIN},{MAX}]")
os.makedirs("Contigs", exist_ok=True)

# one request per 100 genomes; sequence rows come back unordered so they are
# bucketed by genome_id rather than trusted to arrive grouped
proc = subprocess.Popen(
    ["perl", "batch_query.pl", "-c", "genome_sequence", "-i", "genome_id",
     "-r", "genome_id accession description sequence", "-b", "100"],
    stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    text=True)

def feeder():
    for i in ids:
        proc.stdin.write(i + "\n")
    proc.stdin.close()

import threading
threading.Thread(target=feeder, daemon=True).start()

buf = collections.defaultdict(list)
n = 0
for line in proc.stdout:
    p = line.rstrip("\n").split("\t")
    if len(p) < 4 or not p[3].strip():
        continue
    gid, acc, desc, seq = p[0], p[1], p[2], p[3].strip()
    buf[gid].append((acc or gid, desc, seq))
    n += 1
proc.wait()
print(f"{n} contig rows for {len(buf)} genomes")

wrote = 0
with open("Matonaviridae.metadata", "w") as m:
    for gid in ids:
        rows = buf.get(gid)
        if not rows:
            continue
        with open(os.path.join("Contigs", gid + ".contigs"), "w") as fh:
            for acc, desc, seq in rows:
                fh.write(">%s %s\n" % (acc, desc))
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")
        m.write("%s\t%s\t%s\n" % (gid + ".contigs", meta[gid][0], meta[gid][1]))
        wrote += 1
print(f"wrote {wrote} genomes to Contigs/, {len(ids)-wrote} had no sequence")
