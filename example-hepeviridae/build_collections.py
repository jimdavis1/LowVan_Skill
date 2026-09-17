#!/usr/bin/env python3
"""Triage the Hepeviridae annotation strings into feature collections.

Taxon-specific. Everything it knows about Hepeviridae is in RULES, EXPECTED and
MISBINNED; the mechanics follow references/annotation-triage.md.

Four features. Architecture was measured from the coordinates of every
near-complete genome in the family (see PARTITIONING.md): each genus is
ORF1 + ORF2 + ORF3, and Rocahepevirus adds a fourth small ORF.

  ORF1        nonstructural polyprotein, ~1455-1790 aa in the named genera,
              out to ~2504 aa in the unclassified hepe-like viruses
  ORF2        capsid, ~600-680 aa
  ORF3        viroporin, ~87-125 aa
  ORF4_ROCA   Rocahepevirus only -- see the note below

ORF4 IS NOT ONE PROTEIN.  Paslahepevirus ORF4 and Rocahepevirus ORF4 share a
positional label and nothing else: all-vs-all BLASTp over every ORF4-labelled
sequence in the dump gives a best cross-genus alignment of 5 residues at
e-value 138.  Rocahepevirus ORF4 is internally coherent (90-99% identity over
19 sequences) and is collected.  Paslahepevirus ORF4 has 2 unique sequences,
which cannot train a profile at any threshold, and is retired to NOT_MODELLED.

ORF1 IS NOT CUT.  The seven domains RefSeq annotates on ORF1 descend from
Koonin 1992 (PMID 1518855), a computational assignment of putative domains.
Perttila 2013 (PMID 23255617) and LeDesma 2023 (PMID 36852909) both find no
proteolytic processing.  So the domain strings below -- methyltransferase, Y
domain, papain-like protease, poly-proline hinge, X domain, helicase, RdRp --
route to ORF1 and are length-gated like anything else: full-length ones join
the collection, fragments go to ORF1.outliers.fasta where fragments belong.
"""
import re, os, sys, collections, argparse

MODULE = "Hepeviridae"
DUMP   = "Hepeviridae"

# ---------------------------------------------------------------- rules
# (feature_key, genus or None, regex)  -- matched IN ORDER, first hit wins.
# Genus-specific rules come first.
def rx(p): return re.compile(p, re.I)

# ORF4 is NOT scoped by genus, deliberately.  An audit of every ORF1 in the
# collection against every other found 16 of 1417 genomes whose BV-BRC genus
# label contradicts the sequence -- rat and ferret HEV deposited twice under two
# taxon ids, so the same virus appears as both Paslahepevirus and Rocahepevirus
# at 96-99% identity (genus_audit.txt).  A genus-scoped rule would miss the
# mislabelled copies.  HOMOLOGY_GATE below does the discrimination instead, on
# sequence, which is the thing that is actually reliable here.
RULES = [
    ("ORF4_ROCA", None, rx(r"^orf ?-?4\b")),
    ("ORF4_ROCA", "Rocahepevirus", rx(r"^(putative protein|intrinsic disordered protein|zinc-binding protein)$")),

    # --- Avihepevirus names its ORF3 for the phenotype, not the position.
    ("ORF3", "Avihepevirus", rx(r"cytoskeleton[- ]?(related|associated)")),

    # --- ORF3.  Small viroporin.  "phosphoprotein" is unambiguous inside this
    #     family: nothing else in a hepevirus genome is called one.
    ("ORF3", None, rx(r"^p?orf ?-?3\b")),
    ("ORF3", None, rx(r"^orf ?-?3$")),
    ("ORF3", None, rx(r"phosphoprotein")),
    ("ORF3", None, rx(r"^(putative )?accessory protein$")),
    ("ORF3", None, rx(r"^viral accessory (poly)?protein$")),
    ("ORF3", None, rx(r"^viroporin$")),
    ("ORF3", None, rx(r"^(small )?multi-?functional protein$")),

    # --- ORF2.  Capsid.
    ("ORF2", None, rx(r"^orf ?-?2\b")),
    ("ORF2", None, rx(r"^orf ?-?2$")),
    ("ORF2", None, rx(r"\bcapsid")),
    ("ORF2", None, rx(r"^(caspid|capisd|capside|capsid) ?(protein|region)?$")),  # observed typos
    ("ORF2", None, rx(r"^cp( protein)?$")),
    ("ORF2", None, rx(r"^coat protein$")),
    ("ORF2", None, rx(r"^structural (viral )?protein( 2)?( precursor)?$")),
    ("ORF2", None, rx(r"^putative structural protein$")),
    ("ORF2", None, rx(r"^structural polyprotein$")),
    ("ORF2", None, rx(r"^immunogenic protein$")),
    ("ORF2", None, rx(r"^envelope? glycoprotein$")),  # mislabelled; length gate + QC catch it
    ("ORF2", None, rx(r"^vp1$")),

    # --- ORF1.  The nonstructural polyprotein, plus every domain name, because
    #     the domains are not separate proteins in this family.
    ("ORF1", None, rx(r"^orf ?-?1\b")),
    ("ORF1", None, rx(r"^orf ?-?1$")),
    ("ORF1", None, rx(r"non[- ]?struct")),            # nonstructural, non-structural, nonstructual
    ("ORF1", None, rx(r"^mon-structural polyprotein$")),
    ("ORF1", None, rx(r"^non structure protein$")),
    ("ORF1", None, rx(r"polyprotein")),
    ("ORF1", None, rx(r"replicase")),
    ("ORF1", None, rx(r"rna[- ]?(dependent|directed|dep)[- ]?rna[- ]?polymerase")),
    ("ORF1", None, rx(r"^rna[- ]?polymerase$")),
    ("ORF1", None, rx(r"^(truncated )?rdrp\b")),
    ("ORF1", None, rx(r"^polymerase$")),
    ("ORF1", None, rx(r"methyl ?transferase")),
    ("ORF1", None, rx(r"helicase")),
    ("ORF1", None, rx(r"protease")),
    ("ORF1", None, rx(r"^x[- ]?domain( protein)?$")),
    ("ORF1", None, rx(r"^y[- ]?domain( protein)?$")),
    ("ORF1", None, rx(r"poly[- ]?proline")),
    ("ORF1", None, rx(r"^proline-rich region$")),
    ("ORF1", None, rx(r"^hinge( protein)?$")),
    ("ORF1", None, rx(r"^hypervariable region$")),
    ("ORF1", None, rx(r"^nucleotide binding protein$")),
    ("ORF1", None, rx(r"^zinc-binding protein$")),
]

# Strings that must never be gathered by regex: they name a slot or a shrug,
# not a protein.  Homology rescue may still adopt them; a rule may not.
FORBIDDEN = [rx(r"^hypothetical protein"), rx(r"^unknown$"), rx(r"^undefined cds$"),
             rx(r"^unnamed protein product$"), rx(r"^function-unknown protein$"),
             rx(r"^gp\d+ protein$"), rx(r"^blsv$"), rx(r"^vp\d+$"),
             rx(r"^hypothetical proteins$"), rx(r"^transmembrane protein$")]

# ---------------------------------------------------------------- lengths
# Full-length windows.  Anything outside goes to <FEAT>.outliers.fasta rather
# than being discarded: the median record in this family is 347 nt, so most
# sequences are fragments of partial submissions and the outlier files are
# expected to be large.  That is the data, not a fault in the rules.
EXPECTED = {
    "ORF1":      (1380, 2600),
    "ORF2":      ( 560,  700),
    "ORF3":      (  70,  160),
    "ORF4_ROCA": ( 100,  230),
}

# feature_id -> (correct key or None to drop, the evidence)
MISBINNED = {}

# A feature whose members must also be homologous to a verified seed, not merely
# named alike.  ORF4 is the U-number trap in this family: Paslahepevirus ORF4
# (the HEV-1 ER-stress factor, PMID 27035822) and Rocahepevirus ORF4 share a
# positional label and nothing else -- best cross-genus BLASTp alignment over
# every ORF4-labelled sequence in the dump is 5 residues at e-value 138.
# Members failing the gate are written to <FEAT>.rejected.fasta, never dropped
# silently.
HOMOLOGY_GATE = {
    # key       : (seed feature_id,            min % identity)
    "ORF4_ROCA" : ("fig|1213422.7.CDS.2",      50.0),
}

NOT_MODELLED = [
    ("ORF4_PASLA", "Paslahepevirus", 2,
     "Not homologous to Rocahepevirus ORF4 -- best cross-genus BLASTp alignment "
     "is 5 residues at e-value 138. A distinct protein sharing only the label. "
     "2 unique sequences cannot train a profile at any threshold (PMID 27035822)."),
    ("ORF1_DOMAINS", "all", 0,
     "MeT / Y / PCP / poly-proline hinge / X / helicase / RdRp are computational "
     "domain predictions (Koonin 1992, PMID 1518855), not demonstrated mature "
     "peptides. Perttila 2013 (PMID 23255617) and LeDesma 2023 (PMID 36852909) "
     "find no proteolytic processing of pORF1. Domain-named sequences route to "
     "ORF1 and are length-gated."),
]

# ---------------------------------------------------------------- machinery
# Submitters prepend boilerplate that carries no information but defeats an
# anchored rule: "unnamed protein product; ORF1" is an ORF1.
JUNK_PREFIX = re.compile(r"^(unnamed protein product|undefined cds|putative)\s*[;:,]?\s*", re.I)

def normalise(a):
    n = re.sub(r"[-_]", " ", a)
    return re.sub(r"\s+", " ", n).strip()

def strip_junk(a):
    prev = None
    while prev != a:
        prev = a
        a = JUNK_PREFIX.sub("", a).strip()
    return a

def classify(ann, genus):
    for key, g, prx in RULES:
        if g and g != genus: continue
        if prx.search(ann): return key
    for cand in (normalise(ann), strip_junk(ann), normalise(strip_junk(ann))):
        if cand == ann or not cand: continue
        for key, g, prx in RULES:
            if g and g != genus: continue
            if prx.search(cand): return key
    return None

def forbidden(ann):
    return any(p.search(ann) or p.search(normalise(ann)) for p in FORBIDDEN)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    a = ap.parse_args()
    W = a.workdir

    g2g = {}
    for l in open(os.path.join(W, DUMP + ".id_name_gs")):
        p = l.rstrip("\n").split("\t")
        g2g[p[0]] = (p[3] or "unclassified")

    seq = dict(l.rstrip("\n").split("\t") for l in open(os.path.join(W, DUMP + ".uniq.seq")))
    md5 = [l.strip() for l in open(os.path.join(W, DUMP + ".uniq.md5"))]
    ann = [l.rstrip("\n").split("\t") for l in open(os.path.join(W, DUMP + ".uniq.id_ann"))]

    out   = collections.defaultdict(list)
    outl  = collections.defaultdict(list)
    unass = collections.Counter()
    syn   = collections.Counter()
    hit   = collections.Counter()

    for m, (fid, aname) in zip(md5, ann):
        gid   = fid.split("|", 1)[1].rsplit(".", 2)[0]
        genus = g2g.get(gid, "unclassified")
        s     = seq.get(m, "")
        if not s:
            continue
        key = None if forbidden(aname) else classify(aname, genus)
        if fid in MISBINNED:
            key = MISBINNED[fid][0]
            sys.stderr.write("MISBINNED override %s -> %s (%s)\n" % (fid, key, MISBINNED[fid][1]))
        syn[(aname, genus, key or "UNASSIGNED")] += 1
        if key is None:
            unass[(aname, genus)] += 1
            continue
        lo, hi = EXPECTED[key]
        rec = (">%s|%s|%s\n%s\n" % (fid, genus, m, s))
        if lo <= len(s) <= hi:
            out[key].append(rec); hit[key] += 1
        else:
            outl[key].append(rec)

    # --- homology gate -------------------------------------------------
    gate_rejected = collections.defaultdict(list)
    for key, (seed_fid, min_id) in HOMOLOGY_GATE.items():
        recs = out.get(key, [])
        if not recs: continue
        import subprocess, tempfile
        tmp = tempfile.mkdtemp()
        seed = None
        for r in recs:
            if r.split("|", 2)[0][1:] + "|" + r.split("|")[1] == seed_fid:
                seed = r.split("\n")[1]
        if seed is None:
            sys.stderr.write("GATE %s: seed %s not in collection -- gate skipped\n" % (key, seed_fid))
            continue
        open(tmp + "/seed.faa", "w").write(">seed\n%s\n" % seed)
        open(tmp + "/q.faa", "w").writelines(recs)
        subprocess.run(["makeblastdb", "-in", tmp + "/seed.faa", "-dbtype", "prot",
                        "-out", tmp + "/db"], check=True, stdout=subprocess.DEVNULL)
        o = subprocess.run(["blastp", "-query", tmp + "/q.faa", "-db", tmp + "/db",
                            "-outfmt", "6 qseqid pident", "-evalue", "1e-3"],
                           capture_output=True, text=True).stdout
        ok = {l.split("\t")[0] for l in o.splitlines() if float(l.split("\t")[1]) >= min_id}
        keep, drop = [], []
        for r in recs:
            (keep if r[1:].split("\n")[0] in ok else drop).append(r)
        out[key] = keep
        gate_rejected[key] = drop
        hit[key] = len(keep)
        sys.stderr.write("GATE %s: %d kept, %d rejected below %.0f%% identity to %s\n"
                         % (key, len(keep), len(drop), min_id, seed_fid))

    cdir = os.path.join(W, "collections", MODULE)
    os.makedirs(cdir, exist_ok=True)
    for key, recs in gate_rejected.items():
        if recs:
            open(os.path.join(cdir, key + ".rejected.fasta"), "w").writelines(recs)
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
        for r in NOT_MODELLED:
            f.write("%s\t%s\t%d\t%s\n" % r)

    print("collections written to %s\n" % cdir)
    print("%-12s %8s %10s   %s" % ("feature", "in-range", "outliers", "length window"))
    for key in EXPECTED:
        print("%-12s %8d %10d   %d-%d aa" % (key, hit[key], len(outl.get(key, [])), *EXPECTED[key]))
    print("\nunassigned unique sequences : %d  (%d distinct strings)"
          % (sum(unass.values()), len(unass)))
    print("total unique sequences      : %d" % len(md5))

main()
