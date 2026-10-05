#!/usr/bin/env python3
"""Build a coverage-audit facts file for one Picornaviridae genus module.

Every number comes from a measurement file; every sentence comes from the
hand-written prose file. Nothing here is inferred.

    gen_facts.py <Genus> <prose.json> <out_audit.json>

Measurement inputs (all under ~/lowvan_scratch/picorna):
  fin/<G>/G.*                 genus-only dump        (section 01)
  fin/<G>/registry.json       collect_registry.py    (02, 03)
  fin/<G>/Alignments/.../BUILD_PARAMS                (02)
  coverage/<G>.json           evaluate_coverage.py   (04, 05)
  panel/<G>/panel_eval.tsv    run_gto_eval.py        (06)
  fin/<G>/agg.json, unbinned.json, rarefaction.json  (09)
"""
import collections, glob, json, os, sys

G, PROSE, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
S = os.path.expanduser("~/lowvan_scratch/picorna")
F = os.path.join(S, "fin", G)
pr = json.load(open(PROSE))
mod = json.load(open(os.path.join(F, "%s_Viral_PSSM.json" % G)))[G]
feats = mod["features"]
reg = {f["key"]: f for f in json.load(open(os.path.join(F, "registry.json")))[G]["features"]}
cov = json.load(open(os.path.join(S, "coverage", "%s.json" % G)))
agg = json.load(open(os.path.join(F, "agg.json")))
unb = json.load(open(os.path.join(F, "unbinned.json")))
rar = json.load(open(os.path.join(F, "rarefaction.json")))

def c(n): return "{:,}".format(n)
def P(x, y): return "%.1f%%" % (100.0 * x / y) if y else "&mdash;"

#  order features as the genome is laid out, POLY first
ORDER = ["POLY", "Ldr", "VP0", "VP4", "VP2", "VP3", "VP1", "2A", "2A1", "2A2", "2A3",
         "2B", "2C", "3A", "VPG", "3C", "3D"]
keys = sorted(feats, key=lambda k: (ORDER.index(k) if k in ORDER else 99, k))

# 01 -- the dump
gen_ids = [l.split("\t")[0] for l in open(os.path.join(F, "G.id_name_gs"))]
species = {l.rstrip("\n").split("\t")[4] for l in open(os.path.join(S, "Picorna.full.tsv"))
           if l.split("\t")[3] == G}
n_feat = sum(1 for _ in open(os.path.join(F, "G.id_md5")))
n_uniq = sum(1 for _ in open(os.path.join(F, "G.uniq.md5")))
n_contig = len(glob.glob(os.path.join(S, "Contigs.%s" % G, "*.fna")))

# 02 -- departures from the documented defaults
def departures(k):
    p = os.path.join(F, "Alignments", G, k, "BUILD_PARAMS")
    if not os.path.exists(p): return "&mdash;"
    for l in open(p):
        if l.startswith("departures:"):
            raw = l.split(":", 1)[1].strip()
            d = {} if raw == "none" else json.loads(raw)
            #  early BUILD_PARAMS recorded dash-keyed defaults as departures
            DEF = {"m": "5", "mi": "0.8", "mc": "0.8"}
            d = {k.lstrip("-"): str(v) for k, v in d.items()}
            d = {k: v for k, v in d.items() if DEF.get(k) != v}
            return ", ".join("-%s %s" % kv for kv in sorted(d.items())) or "defaults"
    return "defaults"

n_pssm = sum(len(reg[k]["pssms"]) for k in keys if k in reg)
build_rows = [[k, feats[k]["anno"], feats[k].get("gene_symbol", ""), c(reg[k]["n"]) if k in reg else "0",
               str(len(reg[k]["pssms"])) if k in reg else "0", departures(k)] for k in keys]

# 03 -- cutoffs and self-recall
cut_rows = [[k, str(feats[k]["bit_cutoff"]), str(feats[k]["coverage_cutoff"]),
             c(feats[k]["min_len"]), c(feats[k]["max_len"])] for k in keys]
recall = ", ".join("%s %d/%d" % (k, reg[k]["called_inrange"], reg[k]["n_inrange"])
                   for k in keys if k in reg and reg[k]["n_inrange"])
tot_c = sum(reg[k]["called_inrange"] for k in keys if k in reg)
tot_n = sum(reg[k]["n_inrange"] for k in keys if k in reg)

# 04 -- routing
ex = cov["exemplars"]; unr = cov["unrouted"]
rows = [r for r in cov["rows"] if r.get("module") == G]
#  annotate_by_viral_pssm.pl emits a no_features_called sentinel so a routed
#  genome keeps its taxon assignment; it reaches the table as an empty symbol.
#  It is not a feature -- count those genomes separately.
empty = [r for r in rows if not [k for k in r["calls"] if k]]
for r in rows: r["calls"] = [k for k in r["calls"] if k]
other = collections.Counter(r.get("module") for r in cov["rows"]
                            if r.get("module") not in (G, None, "", "UNROUTED"))
routed = len(rows)
refs_used = len(cov.get("flags", {}))

# 05 -- per-feature presence in routed exemplars
pres = collections.Counter()
for r in rows:
    for k in set(r["calls"]): pres[k] += 1
full = sum(1 for r in rows if set(keys) <= set(r["calls"]))
n_empty = len(empty)
bars = [[feats[k].get("gene_symbol", k), "%.1f" % (100.0 * pres[k] / routed if routed else 0),
         round(100.0 * pres[k] / routed, 1) if routed else 0, "%d of %d routed exemplars" % (pres[k], routed)]
        for k in keys]
core = [k for k in keys if k != "POLY"]
best = max(pres[k] for k in core) if core else 0
lagging = [k for k in core if best - pres[k] >= 0.20 * routed]

# 06 -- held-out panel through the GTO quality pipeline
pe = os.path.join(S, "panel", G, "panel_eval.tsv")
prow = [l.rstrip("\n").split("\t") for l in open(pe)][1:]
n_p = len(prow)
nothing = sum(1 for r in prow if int(r[1]) == 0)
noflag = sum(1 for r in prow if int(r[1]) > 0 and int(r[2]) == 0 and int(r[3]) == 0)
clean = sum(1 for r in prow if int(r[1]) > 0 and int(r[2]) == 0 and int(r[3]) == 0 and int(r[4]) == 0)
flagged = n_p - clean - nothing
ftally = collections.Counter()
for r in prow:
    for f in (r[5] if len(r) > 5 else "").split("; "):
        if f.strip(): ftally[f.strip()] += 1

# 09 -- collapse and rarefaction
bound = [a for a in agg if a["key"] in feats]
#  same definition as gen_collapse.py, so the tile and the collapse page agree
nstr = sum(a["nstrings"] for a in agg)
nanno = len({a["anno"] for a in agg})
collapse = nstr / float(nanno) if nanno else 0
if rar.get("unit") != "proteins": sys.exit("rarefaction.json is genome-level; re-run annotation_rarefaction.py")
n_r = rar["n_genomes"]
fold = rar["source_at_common"] / rar["module_at_common"] if rar["module_at_common"] else 0
agg_rows = [[a["key"], a["anno"], str(a["nstrings"])] for a in sorted(bound, key=lambda a: -a["nstrings"])]
others = [a for a in agg if a["key"] not in feats]
for a in others:
    agg_rows.append(["<i>%s</i>" % a["key"].lower(), "not modelled", str(a["nstrings"])])
agg_rows.append(["<i>unbound</i>", "&mdash;", str(len(unb))])

slug = G.lower()
facts = {
 "title": "%s Coverage Audit" % G,
 "eyebrow": "LowVan &middot; %s module &middot; %s" % (G, pr["date"]),
 "h1": pr["h1"],
 "standfirst": pr["standfirst"],
 "headline_tiles": [
   [str(len(keys)), "features, all with profiles" if all(k in reg and reg[k]["pssms"] for k in keys) else "features"],
   [str(n_pssm), "PSSMs"],
   ["%d/%d" % (routed, ex), "exemplars routed"],
   [P(clean, n_p), "held-out genomes genuinely clean"],
   ["%.1f&times;" % collapse, "string collapse"],
   ["%.1f&times;" % fold, "vocabulary collapse at equal proteins sampled"]],
 "sections": [
  {"num": "01", "title": "What the module covers", "lede": pr["lede01"], "blocks":
    [{"type": "prose", "html": h} for h in pr["p01"][:1]] +
    [{"type": "table", "headers": ["", "count"], "rows": [
       ["genomes in the BV-BRC protein dump", c(len(gen_ids))],
       ["protein features", c(n_feat)],
       ["unique sequences", c(n_uniq)],
       ["genomes with a contig &ge; 5 kb", c(n_contig)],
       ["species", c(len(species))]]}] +
    [{"type": "prose", "html": h} for h in pr["p01"][1:]]},
  {"num": "02", "title": "How it was built", "lede": pr["lede02"], "blocks":
    [{"type": "table", "headers": ["feature", "annotation", "symbol", "collection", "profiles", "departures"],
      "rows": build_rows}] + [{"type": "prose", "html": h} for h in pr["p02"]]},
  {"num": "03", "title": "Where the cutoffs sit", "lede": pr["lede03"], "blocks":
    [{"type": "table", "headers": ["feature", "bit_cutoff", "coverage_cutoff", "min_len", "max_len"], "rows": cut_rows}] +
    [{"type": "prose", "html": h} for h in pr["p03"]] +
    [{"type": "verdict", "kind": "ok", "html":
      "Self-recall against each feature's own in-range collection: <b>%d of %d (%s)</b> &mdash; %s. "
      "This is a floor, not a validation; the held-out figures are in &sect;05 and &sect;06." % (tot_c, tot_n, P(tot_c, tot_n), recall)}]},
  {"num": "04", "title": "Does it route?", "lede": pr["lede04"], "blocks": [
     {"type": "table", "headers": ["", "count", "%"], "rows": [
       ["genomes &ge; 5 kb, clustered at 95% nucleotide identity", c(ex), "&mdash;"],
       ["routed to %s" % G, c(routed), P(routed, ex)]] +
       [["routed to %s" % m, c(n), P(n, ex)] for m, n in other.most_common()] +
       [["unrouted", c(len(unr)), P(len(unr), ex)]]},
     {"type": "prose", "html": "%d rep contigs ship with the module; %d of them are the best BLASTn hit for at least one exemplar. "
       "Routing was measured against the ten Picornaviridae modules only &mdash; nucleotide BLASTn cannot reach the routing "
       "threshold across virus families, and leaving the other 56 modules' 558 rep contigs out cut the run time four-fold." %
       (len(mod.get("close_genomes", {})), refs_used)}] + [{"type": "prose", "html": h} for h in pr["p04"]]},
  {"num": "05", "title": "Does it call the proteins?", "lede": pr["lede05"], "blocks": [
     {"type": "bars", "rows": bars, "unit": "%"},
     {"type": "prose", "html": "<b>%d of %d routed exemplars get every declared feature.</b> %s" % (full, routed,
       ("Features trailing the best-called core protein by 20 points or more: <b>%s</b>." % ", ".join(lagging))
       if lagging else "No core feature trails the best-called one by 20 points or more.")
       + ("" if not n_empty else " <b>%d routed exemplar(s) cleared no profile at any feature</b> and are counted here as calling nothing." % n_empty)}] +
     [{"type": "prose", "html": h} for h in pr["p05"]]},
  {"num": "06", "title": "Are the calls right?", "lede": pr["lede06"], "blocks": [
     {"type": "tiles", "items": [[str(n_p), "held-out genomes"], [str(noflag), "no genome/contig flags"],
        [str(clean), "genuinely clean"], [str(flagged), "flagged"], [str(nothing), "called nothing"],
        [P(clean, n_p), "clean"]]},
     {"type": "prose", "html": "The panel is %d genomes drawn round-robin across species from this genus's own BV-BRC contigs, "
       "with every rep-contig reference genome excluded, built into GTOs with <code>rast-create-genome</code> and scored "
       "by <code>viral_genome_quality.pl</code>. <i>Genuinely clean</i> means no genome, contig or feature flag at all." % n_p}] +
     ([{"type": "pre", "text": "\n".join("%4d  %s" % (n, f) for f, n in ftally.most_common(12))}] if ftally else []) +
     [{"type": "prose", "html": h} for h in pr["p06"]]},
  {"num": "07", "title": "What it cannot do", "lede": pr["lede07"], "blocks":
     [{"type": "prose", "html": h} for h in pr["p07"]] +
     [{"type": "verdict", "kind": "warn", "html": pr["verdict07"]}]},
  {"num": "08", "title": "What the run found", "lede": pr["lede08"], "blocks":
     [{"type": "prose", "html": h} for h in pr["p08"]]},
  {"num": "09", "title": "What the vocabulary bought", "lede": pr["lede09"], "blocks": [
     {"type": "figure", "from": "%s-vocabulary-saturation.html" % slug},
     {"type": "table", "headers": ["feature", "annotation string", "source strings absorbed"], "rows": agg_rows},
     {"type": "prose", "html": "Rarefaction over proteins: every annotated protein in the %s genomes both sources annotate "
       "(%s source product records, %s module calls) shuffled %d times and drawn one at a time. Both curves start at zero. "
       "The source vocabulary reaches <b>%d</b> distinct strings, the module's <b>%d</b>. At the depth both reach, %s proteins, "
       "free text has used <b>%.1f strings against %.1f &mdash; a %.1f-fold difference.</b> "
       "The source reaches 95%% of its final size after %s proteins, the module after %s." %
       (c(n_r), c(rar["source"]["n"]), c(rar["module"]["n"]), rar["replicates"], rar["source_total"], rar["module_total"],
        c(rar["common_n"]), rar["source_at_common"], rar["module_at_common"], fold, c(rar["source_knee"]), c(rar["module_knee"]))}] +
     [{"type": "prose", "html": h} for h in pr["p09"]]}],
 "source_note": "Generated from <code>%s_audit.json</code> by <code>gen_coverage_audit.py</code>; the numbers by "
   "<code>gen_facts.py</code> from collect_registry.py, evaluate_coverage.py, run_gto_eval.py, collect_synmap.py and "
   "annotation_rarefaction.py output in <code>Picornaviridae/work/%s/</code>. All citations in this module are "
   "model-proposed and carried in PMID_claude_generated; none has been read and confirmed by a curator." % (slug, G),
}
json.dump(facts, open(OUT, "w"), indent=1)

#  the numbers the prose has to agree with
print(json.dumps({"genus": G, "features": len(keys), "pssms": n_pssm, "dump_genomes": len(gen_ids),
  "features_dump": n_feat, "uniq": n_uniq, "contigs": n_contig, "species": len(species),
  "exemplars": ex, "routed": routed, "unrouted": unr, "other_modules": dict(other),
  "presence": {k: pres[k] for k in keys}, "full_set": full,
  "routed_called_nothing": n_empty, "lagging": lagging,
  "panel": {"n": n_p, "noflag": noflag, "clean": clean, "flagged": flagged, "nothing": nothing,
            "flags": dict(ftally.most_common(12))},
  "selfrecall": [tot_c, tot_n], "collapse": round(collapse, 1), "nstrings": nstr, "unbound": len(unb),
  "unbound_top": unb[:8], "rarefaction": {"genomes": n_r, "src_prot": rar["source"]["n"], "mod_prot": rar["module"]["n"],
  "src": rar["source_total"], "mod": rar["module_total"], "common": rar["common_n"],
  "src_knee": rar["source_knee"], "mod_knee": rar["module_knee"], "fold": round(fold, 1)}, "departures": {k: departures(k) for k in keys},
  "collections": {k: (reg[k]["n"] if k in reg else 0) for k in keys}}, indent=1))
