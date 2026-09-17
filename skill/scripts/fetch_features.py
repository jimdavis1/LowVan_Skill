#!/usr/bin/env python3
"""Batched BV-BRC genome_feature dump.

query_PATRIC_bob.pl issues one HTTP query per input line and rebuilds the API
object each time; 87k genomes takes ~4 h.  This asks for 400 genome ids at a
time over a thread pool and takes ~3 min for the same input.
"""
import sys, json, time, threading, argparse
import requests
from concurrent.futures import ThreadPoolExecutor

URL = "https://www.bv-brc.org/api/genome_feature/"
HDR = {"Content-Type": "application/rqlquery+x-www-form-urlencoded",
       "Accept": "application/json"}
FIELDS = "genome_id,patric_id,aa_sequence_md5,product,feature_type,length,start,end,strand,accession"
LIMIT  = 25000

lock = threading.Lock()
done = [0]

def chunk(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i+n]

def fetch(ids, retries=5):
    q = "in(genome_id,(%s))&select(%s)&limit(%d)" % (",".join(ids), FIELDS, LIMIT)
    for attempt in range(retries):
        try:
            r = requests.post(URL, data=q, headers=HDR, timeout=300)
            if r.status_code == 200:
                rows = r.json()
                if len(rows) >= LIMIT:
                    sys.stderr.write("WARN batch hit LIMIT (%d) -- splitting\n" % LIMIT)
                    if len(ids) == 1:
                        raise RuntimeError("single genome exceeds LIMIT: %s" % ids[0])
                    h = len(ids)//2
                    return fetch(ids[:h]) + fetch(ids[h:])
                return rows
            sys.stderr.write("HTTP %d attempt %d\n" % (r.status_code, attempt+1))
        except Exception as e:
            sys.stderr.write("ERR %s attempt %d\n" % (e, attempt+1))
        time.sleep(2*(attempt+1))
    raise RuntimeError("batch failed after %d retries" % retries)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True, help="file of genome ids, one per line")
    ap.add_argument("--out", required=True, help="output TSV")
    ap.add_argument("--batch", type=int, default=400)
    ap.add_argument("--jobs",  type=int, default=12)
    a = ap.parse_args()

    ids = [l.strip() for l in open(a.ids) if l.strip()]
    batches = list(chunk(ids, a.batch))
    sys.stderr.write("%d genomes in %d batches\n" % (len(ids), len(batches)))

    cols = FIELDS.split(",")
    with open(a.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        with ThreadPoolExecutor(max_workers=a.jobs) as ex:
            for rows in ex.map(fetch, batches):
                for r in rows:
                    out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
                with lock:
                    done[0] += 1
                    if done[0] % 20 == 0:
                        sys.stderr.write("  %d/%d batches\n" % (done[0], len(batches)))
    sys.stderr.write("done\n")

main()
