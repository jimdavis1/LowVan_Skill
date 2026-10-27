#!/usr/bin/env python3
"""Is each deposited protein collinear with its own genome, or spliced?

A module declares a PSSM feature only for proteins that are a single ORF on the
genome. One that needs two exons is not collinear and gets a splice reference
set instead; one made by a frameshift gets transcript-edit references. Both are
cheap to get wrong by analogy -- SRLV *tat* is a single exon and gets a profile
while EIAV *tat*, the same name in a sister module, is spliced and gets none --
so the project rule is to MEASURE this per taxon rather than carry the layout
over from a relative.

The measurement: tblastn each deposited protein against the contigs of the
genome it came from, and count how many HSPs are needed to cover the query.

    one HSP covering most of the query      collinear  -> PSSM feature
    two HSPs, same strand, in query order   spliced    -> splice references
    two HSPs, different frames, overlapping frameshift -> transcript edit

Reports per annotation string, because that is the unit a collection is built
from. A string that is collinear in most genomes and spliced in a few is worth
looking at by hand: it usually means two different proteins share a name.

    python3 check_collinearity.py --contigs Contigs.Alpharetrovirus \
        --workdir . --prefix R --genus Alpharetrovirus
"""
import argparse, collections, glob, json, os, re, subprocess, sys, tempfile


def load_dump(workdir, prefix):
    md5 = [l.strip() for l in open(os.path.join(workdir, "%s.uniq.md5" % prefix))]
    ann = [l.rstrip("\n").split("\t")
           for l in open(os.path.join(workdir, "%s.uniq.id_ann" % prefix), errors="replace")]
    m2a = {m: " ".join((a[1] if len(a) > 1 else "").split()) for m, a in zip(md5, ann)}
    m2s = {}
    for l in open(os.path.join(workdir, "%s.uniq.seq" % prefix), errors="replace"):
        f = l.rstrip("\n").split("\t")
        if len(f) >= 2:
            m2s[f[0]] = f[1]
    g2gen = {}
    p = os.path.join(workdir, "%s.id_name_gs" % prefix)
    if os.path.exists(p):
        for l in open(p, errors="replace"):
            f = l.rstrip("\n").split("\t")
            if len(f) >= 4:
                g2gen[f[0]] = f[3]
    return m2a, m2s, g2gen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs", required=True, help="directory of per-genome FASTA")
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--prefix", required=True, help="dump prefix, e.g. R")
    ap.add_argument("--genus", default=None,
                    help="comma-separated genus filter from <prefix>.id_name_gs")
    ap.add_argument("--min-cov", type=float, default=0.90,
                    help="fraction of the query one HSP must cover to be called "
                         "collinear (default 0.90)")
    ap.add_argument("--max-genomes", type=int, default=0,
                    help="sample at most this many genomes per genus (0 = all)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", default=None, help="write the per-string table as JSON")
    a = ap.parse_args()

    m2a, m2s, g2gen = load_dump(a.workdir, a.prefix)
    want = set(x.strip() for x in a.genus.split(",")) if a.genus else None

    #  genome -> [(annotation, sequence)] from id_md5
    per = collections.defaultdict(list)
    for l in open(os.path.join(a.workdir, "%s.id_md5" % a.prefix), errors="replace"):
        f = l.rstrip("\n").split("\t")
        if len(f) < 2 or not f[1].strip():
            continue
        g = ".".join(f[0].split("|")[-1].split(".")[:2])
        if want and g2gen.get(g) not in want:
            continue
        s, q = m2a.get(f[1], ""), m2s.get(f[1], "")
        if s and q:
            per[g].append((s, q))

    genomes = sorted(per)
    if a.max_genomes:
        genomes = genomes[:a.max_genomes]
    print("  %d genome(s) with protein features%s"
          % (len(genomes), " (sampled)" if a.max_genomes else ""))

    verdicts = collections.defaultdict(collections.Counter)
    nhsp = collections.defaultdict(list)
    tested = missing = 0

    for g in genomes:
        fa = os.path.join(a.contigs, "%s.fna" % g)
        if not os.path.exists(fa):
            for alt in ("%s.fasta" % g, "%s.fa" % g):
                if os.path.exists(os.path.join(a.contigs, alt)):
                    fa = os.path.join(a.contigs, alt); break
        if not os.path.exists(fa) or os.path.getsize(fa) == 0:
            missing += 1
            continue
        with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as fh:
            qf = fh.name
            for i, (s, q) in enumerate(per[g]):
                fh.write(">q%d\n%s\n" % (i, q))
        try:
            r = subprocess.run(
                ["tblastn", "-query", qf, "-subject", fa, "-outfmt", "15",
                 "-seg", "no", "-num_threads", str(a.threads), "-evalue", "1e-5"],
                capture_output=True, text=True)
            if r.returncode != 0:
                missing += 1
                continue
            B = json.loads(r.stdout)
        except Exception:
            missing += 1
            continue
        finally:
            os.unlink(qf)

        #  TWO RESULT SHAPES. With -db the per-query block is results->search,
        #  a dict. With -subject BLAST switches to bl2seq mode and results
        #  holds a LIST under "bl2seq", one entry per query, with no "search"
        #  key at all. Reading only "search" silently finds nothing and the
        #  run reports "0 proteins tested, 0 skipped" -- no error, no hits.
        for bo in B.get("BlastOutput2") or []:
            rep = (bo.get("report") or {}).get("results") or {}
            blocks = rep.get("bl2seq")
            if not isinstance(blocks, list):
                blocks = [rep.get("search") or {}]
            for srch in blocks:
                qid = srch.get("query_title") or srch.get("query_id") or ""
                m = re.match(r"q(\d+)$", qid.strip())
                if not m:
                    continue
                s = per[g][int(m.group(1))][0]
                qlen = srch.get("query_len") or 0
                hsps = []
                for h in srch.get("hits") or []:
                    hsps.extend(h.get("hsps") or [])
                if not hsps or not qlen:
                    verdicts[s]["no hit"] += 1
                    continue
                hsps.sort(key=lambda x: -x.get("bit_score", 0))
                top = hsps[0]
                cov = (top["query_to"] - top["query_from"] + 1) / qlen
                tested += 1
                if cov >= a.min_cov:
                    verdicts[s]["collinear"] += 1
                    nhsp[s].append(1)
                    continue
                #  Not covered by one HSP. Two on the same strand, in query
                #  order and non-overlapping on the query, is the splice
                #  signature. tblastn HSPs carry no hit_strand field; the
                #  sign of hit_frame is the strand.
                def strand(h):
                    return 1 if (h.get("hit_frame") or 1) > 0 else -1
                same = [h for h in hsps if strand(h) == strand(top)]
                same.sort(key=lambda x: x["query_from"])
                span, used, last_q = 0, [], 0
                for h in same:
                    if h["query_from"] > last_q:
                        span += h["query_to"] - h["query_from"] + 1
                        used.append(h); last_q = h["query_to"]
                tot = span / qlen if qlen else 0
                frames = {h.get("hit_frame") for h in used}
                if tot >= a.min_cov and len(used) >= 2:
                    verdicts[s]["spliced/multi-exon" if len(frames) == 1
                                else "frameshift"] += 1
                    nhsp[s].append(len(used))
                else:
                    verdicts[s]["partial"] += 1
                    nhsp[s].append(len(used) or 1)

    print("  %d protein(s) tested, %d genome(s) skipped (no contig or BLAST error)\n"
          % (tested, missing))
    rows = []
    print("  %-44s %5s  %-22s %s" % ("annotation string", "n", "verdict (majority)", "breakdown"))
    for s, c in sorted(verdicts.items(), key=lambda x: -sum(x[1].values())):
        n = sum(c.values())
        top = c.most_common(1)[0]
        rows.append({"annotation": s, "n": n, "verdict": top[0],
                     "pct": round(100.0 * top[1] / n, 1), "breakdown": dict(c)})
        print("  %-44s %5d  %-22s %s"
              % (s[:44], n, "%s %.0f%%" % (top[0], 100.0 * top[1] / n),
                 ", ".join("%s %d" % (k, v) for k, v in c.most_common())))
    if a.out:
        json.dump(rows, open(a.out, "w"), indent=1)
        print("\n  wrote %s" % a.out)


main()
