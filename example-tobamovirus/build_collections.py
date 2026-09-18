#!/usr/bin/env python3
"""Triage Tobamovirus into feature collections.

Four features, and the interesting one is that TWO OF THEM SHARE THEIR
ANNOTATION STRINGS. Tobamovirus replicase is expressed as a 126K protein that
terminates at a leaky amber, and a 183K readthrough product that continues
through it. BV-BRC calls both "replicase", "RNA-dependent RNA polymerase",
"replication protein" and so on, so no text rule can separate them.

Length can, and the split is unambiguous. Over the 1,207 full-length
replicase-like sequences in the dump:

    1100-1199 aa   578 sequences      <- the opal-terminated 126K form
    1200-1399 aa     5 sequences      <- the valley
    1500-1699 aa   607 sequences      <- the 183K readthrough form

    median short 1116, median long 1616, difference exactly 500 aa

So the rules bin on `(string, length)` and the cut is 1400. Anything outside
either window is a fragment and goes to the outliers, as usual.

REP183 carries `internal_stop: 1` in the JSON and its alignments must span the
amber, or the profile has never seen what follows it. REP126 is the terminated
product and is a real protein in its own right, exactly as alphavirus P123 is
alongside nsP1234 -- its C-terminus IS the stop.
"""
import re, os, sys, collections, argparse

MODULE = "Tobamovirus"
DUMP   = "Tobamovirus"

def rx(p): return re.compile(p, re.I)

#  Replicase strings, which apply to BOTH forms. Length decides which.
REPLICASE = rx(r"replicas|replicat|polymeras|\brdrp\b|\brna pol|"
               r"1[0-9]{2} ?k(da)?( protein| replicase)?|p1[0-9]{2}\b|"
               r"helicase|methyltransferase|readthrough|read.through")
MOVEMENT  = rx(r"movement|cell.to.cell|\bmp\b|30 ?kda?( protein)?$|\bp30\b")
COAT      = rx(r"\bcoat\b|capsid|\bcp\b|17.[0-9] ?kda|\bp17\b")

#  Length windows. REP126/REP183 from the bimodal distribution above; MP and CP
#  from their modal bins (MP mode 260-279, 81% within +/-120; CP mode 140-159,
#  100% within +/-120).
#  REP126's top bound sits BELOW the valley, not across it. An earlier build
#  used 1000-1399, which admitted the 5 sequences between the two modes; those
#  are longer than a genuine 126K and one of them produced a leftover profile
#  (lo12, master 1220 aa) that outscored every real 126K profile and pulled the
#  call 570 nt downstream into the readthrough region. The behavioural check on
#  TMV caught it: 126K came out 926 aa starting at 639 instead of 1116 at 69.
EXPECTED = {
    "REP126": (980, 1199),
    "REP183": (1400, 1780),
    "MP":     ( 200,  360),
    "CP":     ( 130,  200),
}
REP_SPLIT = 1300   # anything above this is the readthrough form

FORBIDDEN = [rx(r"^hypothetical protein"), rx(r"^unknown"), rx(r"^orf ?\d+$"),
             rx(r"^unnamed protein product$"), rx(r"^putative protein$")]

MISBINNED = {}
NOT_MODELLED = []

def normalise(a):
    n = re.sub(r"[-_]", " ", a)
    return re.sub(r"\s+", " ", n).strip()

def classify(ann, n):
    """Return a feature key. Replicase strings are split on length."""
    for cand in (ann, normalise(ann)):
        if REPLICASE.search(cand):
            return "REP126" if n < REP_SPLIT else "REP183"
        if MOVEMENT.search(cand):  return "MP"
        if COAT.search(cand):      return "CP"
    return None

def forbidden(a):
    return any(p.search(a) or p.search(normalise(a)) for p in FORBIDDEN)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workdir", default="."); a = ap.parse_args()
    W = a.workdir
    g2g = {}
    for l in open(os.path.join(W, DUMP + ".id_name_gs")):
        p = l.rstrip("\n").split("\t"); g2g[p[0]] = p[3] or "unclassified"
    seq = dict(l.rstrip("\n").split("\t") for l in open(os.path.join(W, DUMP + ".uniq.seq")))
    md5 = [l.strip() for l in open(os.path.join(W, DUMP + ".uniq.md5"))]
    ann = [l.rstrip("\n").split("\t") for l in open(os.path.join(W, DUMP + ".uniq.id_ann"))]

    out = collections.defaultdict(list); outl = collections.defaultdict(list)
    unass = collections.Counter(); syn = collections.Counter(); hit = collections.Counter()
    for m, (fid, aname) in zip(md5, ann):
        gid = fid.split("|", 1)[1].rsplit(".", 2)[0]
        genus = g2g.get(gid, "Tobamovirus"); s = seq.get(m, "")
        if not s: continue
        key = None if forbidden(aname) else classify(aname, len(s))
        if fid in MISBINNED: key = MISBINNED[fid][0]
        syn[(aname, genus, key or "UNASSIGNED")] += 1
        if key is None: unass[(aname, genus)] += 1; continue
        lo, hi = EXPECTED[key]
        rec = ">%s|%s|%s\n%s\n" % (fid, genus, m, s)
        if lo <= len(s) <= hi: out[key].append(rec); hit[key] += 1
        else: outl[key].append(rec)

    cdir = os.path.join(W, "collections", MODULE); os.makedirs(cdir, exist_ok=True)
    for key in EXPECTED:
        open(os.path.join(cdir, key + ".fasta"), "w").writelines(out.get(key, []))
        if outl.get(key):
            open(os.path.join(cdir, key + ".outliers.fasta"), "w").writelines(outl[key])
    cd = os.path.join(W, "collections")
    with open(os.path.join(cd, "synonyms.tsv"), "w") as f:
        f.write("count\tannotation\tgenus\tbins_to\n")
        for (s_, g, k), c in sorted(syn.items(), key=lambda x: -x[1]):
            f.write("%d\t%s\t%s\t%s\n" % (c, s_, g, k))
    with open(os.path.join(cd, "UNASSIGNED_TRACKING.tsv"), "w") as f:
        f.write("count\tannotation\tgenus\n")
        for (s_, g), c in sorted(unass.items(), key=lambda x: -x[1]):
            f.write("%d\t%s\t%s\n" % (c, s_, g))
    with open(os.path.join(cd, "NOT_MODELLED.tsv"), "w") as f:
        f.write("feature\ttaxon\tn_sequences\treason\n")
        for r in NOT_MODELLED: f.write("%s\t%s\t%d\t%s\n" % r)

    print("%-10s %8s %10s   %s" % ("feature", "in-range", "outliers", "window"))
    for key in EXPECTED:
        print("%-10s %8d %10d   %d-%d aa" % (key, hit[key], len(outl.get(key, [])), *EXPECTED[key]))
    print("\nunassigned: %d sequences (%d strings)   total: %d"
          % (sum(unass.values()), len(unass), len(md5)))
main()
