#!/usr/bin/env python3
"""Batched BV-BRC fetch, as a drop-in for the per-line query_PATRIC_bob.pl loop.

query_PATRIC_bob.pl opens a new P3DataAPI object and issues one `eq` query for
every input line. On the Hepacivirus dump that was 268 rows/min: 65,735 genome
lookups took 4h17m and the 100,855 annotation lookups were tracking to six
hours more. The API accepts `in(field,(v1,v2,...))`, which is what
repcontig_budget.py already uses, so the fix is batching rather than
parallelism -- ~100x fewer round trips instead of ~6x.

Resumable: reads whatever the output file already holds and asks only for the
keys still missing, so it can take over from a partial serial run.

    python3 fetch_bvbrc.py --core genome_feature --key patric_id \
        --fields patric_id,product --in ids.txt --out out.tsv
    python3 fetch_bvbrc.py --core feature_sequence --key md5 \
        --fields md5,sequence --in md5.txt --out seq.tsv
"""
import argparse, os, sys, time, urllib.parse
from concurrent.futures import ThreadPoolExecutor
import requests

API = "https://www.bv-brc.org/api/%s/"
HDR = {"Content-Type": "application/rqlquery+x-www-form-urlencoded",
       "Accept": "application/json"}


def quoted(v):
    #  Percent-encode, which is what query_PATRIC_bob.pl does with uri_escape.
    #  Measured against the API on patric_ids, which contain a pipe:
    #    bare              HTTP 400
    #    "double-quoted"   HTTP 400
    #    percent-encoded   HTTP 200   <-
    #    "percent-encoded" HTTP 400
    #  md5s are hex and work either way, so one encoded path serves both.
    return urllib.parse.quote(v, safe="")


def fetch(core, key, fields, batch, tries=4):
    q = "in(%s,(%s))&select(%s)&limit(%d)" % (
        key, ",".join(quoted(b) for b in batch), ",".join(fields), len(batch) * 4 + 100)
    for a in range(tries):
        try:
            r = requests.post(API % core, data=q, headers=HDR, timeout=300)
            if r.status_code == 200:
                return r.json()
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(2 ** a)
                continue
            sys.stderr.write("HTTP %d on a batch of %d\n" % (r.status_code, len(batch)))
            return []
        except Exception as e:
            if a == tries - 1:
                sys.stderr.write("failed batch of %d: %s\n" % (len(batch), e))
                return []
            time.sleep(2 ** a)
    return []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--core", required=True)
    ap.add_argument("--key", required=True)
    ap.add_argument("--fields", required=True, help="comma-separated, written in this order")
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--batch", type=int, default=150)
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    fields = a.fields.split(",")

    want = [l.strip() for l in open(a.inp) if l.strip()]
    seen = set()
    if os.path.exists(a.out):
        for line in open(a.out, errors="replace"):
            p = line.split("\t")
            if p and p[0].strip():
                seen.add(p[0].strip())
    todo = [w for w in want if w not in seen]
    print("  %d wanted, %d already present, %d to fetch" % (len(want), len(seen), len(todo)),
          flush=True)
    if not todo:
        return 0

    batches = [todo[i:i + a.batch] for i in range(0, len(todo), a.batch)]
    t0 = time.time()
    done = rows = 0
    with open(a.out, "a") as fh, ThreadPoolExecutor(max_workers=a.jobs) as pool:
        for recs in pool.map(lambda b: fetch(a.core, a.key, fields, b), batches):
            for x in recs:
                fh.write("\t".join(str(x.get(f, "")) for f in fields) + "\n")
                rows += 1
            done += 1
            if done % 20 == 0:
                el = time.time() - t0
                pct = 100.0 * done / len(batches)
                rate = (done * a.batch) / (el / 60.0)
                eta = (len(batches) - done) * a.batch / max(rate, 1)
                print("    %d/%d batches (%.0f%%)  %d rows  %.0f keys/min  eta %.0f min"
                      % (done, len(batches), pct, rows, rate, eta), flush=True)
                fh.flush()
    print("  wrote %d rows in %.1f min" % (rows, (time.time() - t0) / 60.0))
    return 0


sys.exit(main())
