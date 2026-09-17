#!/usr/bin/env python3
"""Fetch protein sequences by md5 from feature_sequence, batched.

Emits md5<TAB>sequence in exactly the order of the input md5 file, which is what
keeps <T>.uniq.md5 / .uniq.id_ann / .uniq.seq line-index aligned.
"""
import sys, time, argparse, threading
import requests
from concurrent.futures import ThreadPoolExecutor

URL = "https://www.bv-brc.org/api/feature_sequence/"
HDR = {"Content-Type": "application/rqlquery+x-www-form-urlencoded",
       "Accept": "application/json"}
LIMIT = 25000

def fetch(md5s, retries=5):
    q = "in(md5,(%s))&select(md5,sequence)&limit(%d)" % (",".join(md5s), LIMIT)
    for attempt in range(retries):
        try:
            r = requests.post(URL, data=q, headers=HDR, timeout=300)
            if r.status_code == 200:
                return {d["md5"]: d.get("sequence","") for d in r.json()}
            sys.stderr.write("HTTP %d attempt %d\n" % (r.status_code, attempt+1))
        except Exception as e:
            sys.stderr.write("ERR %s attempt %d\n" % (e, attempt+1))
        time.sleep(2*(attempt+1))
    raise RuntimeError("batch failed")

ap = argparse.ArgumentParser()
ap.add_argument("--md5", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--batch", type=int, default=300)
ap.add_argument("--jobs",  type=int, default=12)
a = ap.parse_args()

order = [l.strip() for l in open(a.md5) if l.strip()]
batches = [order[i:i+a.batch] for i in range(0, len(order), a.batch)]
sys.stderr.write("%d md5s in %d batches\n" % (len(order), len(batches)))

seqs = {}
n = 0
with ThreadPoolExecutor(max_workers=a.jobs) as ex:
    for d in ex.map(fetch, batches):
        seqs.update(d); n += 1
        if n % 20 == 0: sys.stderr.write("  %d/%d\n" % (n, len(batches)))

miss = [m for m in order if m not in seqs]
sys.stderr.write("missing sequences: %d\n" % len(miss))
if miss: open(a.out + ".missing", "w").write("\n".join(miss) + "\n")
with open(a.out, "w") as out:
    for m in order:
        out.write("%s\t%s\n" % (m, seqs.get(m, "")))
sys.stderr.write("done\n")
