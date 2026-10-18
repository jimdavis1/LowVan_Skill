#!/usr/bin/env python3
"""The BV-BRC input dump for a family, via the REST API.

Same five-file contract as query_PATRIC_bob.pl, which queries one identifier per
request: 190,698 Sedoreoviridae genomes measured out at ~12 hours for the
feature step alone. Batched through the same SOLR cores this runs in minutes.

The three uniq.* files MUST stay line-index aligned -- that index is the join
key for everything downstream. Here alignment is not a side effect of the server
returning rows in input order; uniq.md5 is written first and id_ann and seq are
emitted by walking that list, so they cannot drift.

Step 3 uses a NEGATIVE filter (drop blank md5) rather than naming feature types.
`grep "CDS\\|mat"` without `peg` once dropped 45% of Rhabdoviridae protein
diversity with no error anywhere.
"""
import argparse, collections, os, sys, time, urllib.parse
import requests

API = "https://www.bv-brc.org/api"
HDR = {"Content-Type": "application/rqlquery+x-www-form-urlencoded",
       "Accept": "application/json"}


def post(core, q, retries=5):
    for a in range(retries):
        try:
            r = requests.post("%s/%s/" % (API, core), data=q, headers=HDR, timeout=600)
            if r.status_code == 200:
                return [x for x in r.json() if isinstance(x, dict)]
            time.sleep(2 * (a + 1))
        except Exception:
            time.sleep(2 * (a + 1))
    raise SystemExit("  query failed after %d tries: %s %s" % (retries, core, q[:120]))


def enc(v):
    """patric_id carries a '|', which RQL rejects raw and rejects quoted. Only
    percent-encoding is accepted by the in() list."""
    return urllib.parse.quote(v, safe="")


def batched(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


ap = argparse.ArgumentParser()
ap.add_argument("--family", required=True)
ap.add_argument("--prefix", required=True, help="output file prefix, e.g. Sedo")
ap.add_argument("--workdir", default=".")
ap.add_argument("--batch", type=int, default=200)
#  A family-wide dump is the right default, but some families are dominated by
#  amplicons rather than genomes and a blind dump is then mostly junk.
#  Retroviridae holds 1,212,379 Lentivirus records of which 93.5% are under
#  3 kb -- HIV drug-resistance surveillance -- against 31,854 at >=8 kb. The
#  module only ever uses the whole genomes, so pulling the rest costs hours
#  and buys nothing.
ap.add_argument("--min-len", type=int, default=0,
                help="skip genomes shorter than this (nt). Use for families "
                     "dominated by partial records; 0 disables.")
ap.add_argument("--genus", default="",
                help="restrict to one genus within the family.")
a = ap.parse_args()
W, T = a.workdir, a.prefix


def path(suffix): return os.path.join(W, "%s.%s" % (T, suffix))


def log(msg): print("  %s  %s" % (time.strftime("%H:%M"), msg), flush=True)


#  1. genomes -------------------------------------------------------------
if not os.path.exists(path("id_name_gs")) or os.path.getsize(path("id_name_gs")) == 0:
    log("1. genomes")
    rows, got = [], 0
    while True:
        sel = "eq(family,%s)" % a.family
        if a.genus:
            sel = "and(%s,eq(genus,%s))" % (sel, a.genus)
        if a.min_len:
            sel = "and(%s,gt(genome_length,%d))" % (sel, a.min_len - 1)
        r = post("genome", "%s&select(genome_id,genome_name,family,genus)"
                           "&limit(25000,%d)&sort(+genome_id)" % (sel, got))
        if not r: break
        rows += r; got += len(r)
        if len(r) < 25000: break
    with open(path("id_name_gs"), "w") as f:
        for x in rows:
            f.write("\t".join([x.get("genome_id", ""), x.get("genome_name", ""),
                               x.get("family", ""), x.get("genus", "")]) + "\n")
    log("   %d genomes" % len(rows))
gids = [l.split("\t")[0] for l in open(path("id_name_gs")) if l.strip()]
log("1. genomes: %d" % len(gids))

#  2. features -> patric_id + md5 ------------------------------------------
log("2. features -> id + md5")
n = 0
if os.path.exists(path("id_md5")) and os.path.getsize(path("id_md5")) > 0:
    n = sum(1 for _ in open(path("id_md5")))
    log("   reusing existing id_md5: %d features" % n)
else:
  if True:
   with open(path("id_md5"), "w") as f:
    for i, chunk in enumerate(batched(gids, a.batch)):
        for x in post("genome_feature",
                      "in(genome_id,(%s))&select(patric_id,aa_sequence_md5)&limit(100000)"
                      % ",".join(chunk)):
            pid = x.get("patric_id", "")
            if pid.startswith("fig|"):
                f.write("%s\t%s\n" % (pid, x.get("aa_sequence_md5", "")))
                n += 1
        if i % 100 == 0: log("   batch %d/%d, %d features" % (i, len(gids)//a.batch + 1, n))
log("   %d features" % n)

#  3. distinct sequences ----------------------------------------------------
log("3. distinct sequences (negative filter: drop blank md5)")
md5s = sorted({l.rstrip("\n").split("\t")[1] for l in open(path("id_md5"))
               if len(l.rstrip("\n").split("\t")) > 1 and l.rstrip("\n").split("\t")[1].strip()})
with open(path("uniq.md5"), "w") as f:
    for m in md5s: f.write(m + "\n")
log("   %d distinct sequences" % len(md5s))

#  4. one representative feature id + product per sequence ------------------
log("4. representative id + product")
rep = {}
for l in open(path("id_md5")):
    p = l.rstrip("\n").split("\t")
    if len(p) > 1 and p[1] and p[1] not in rep: rep[p[1]] = p[0]
prod = {}
ids = [rep[m] for m in md5s if m in rep]
for i, chunk in enumerate(batched(ids, a.batch)):
    for x in post("genome_feature",
                  "in(patric_id,(%s))&select(patric_id,product)&limit(100000)"
                  % ",".join(enc(c) for c in chunk)):
        prod[x.get("patric_id", "")] = (x.get("product") or "").replace("\t", " ")
    if i % 100 == 0: log("   batch %d/%d" % (i, len(ids)//a.batch + 1))
with open(path("uniq.id_ann"), "w") as f:          # walk md5s -> alignment cannot drift
    for m in md5s:
        r = rep.get(m, "")
        f.write("%s\t%s\n" % (r, prod.get(r, "")))
log("   %d rows" % len(md5s))

#  5. the sequences ---------------------------------------------------------
log("5. sequences")
seq = {}
for i, chunk in enumerate(batched(md5s, a.batch)):
    for x in post("feature_sequence",
                  "in(md5,(%s))&select(md5,sequence)&limit(100000)" % ",".join(chunk)):
        seq[x.get("md5", "")] = x.get("sequence", "")
    if i % 100 == 0: log("   batch %d/%d" % (i, len(md5s)//a.batch + 1))
with open(path("uniq.seq"), "w") as f:             # same walk, same order
    for m in md5s: f.write("%s\t%s\n" % (m, seq.get(m, "")))
miss = sum(1 for m in md5s if not seq.get(m))
log("   %d rows, %d with no sequence" % (len(md5s), miss))
print("DUMP_COMPLETE", flush=True)
