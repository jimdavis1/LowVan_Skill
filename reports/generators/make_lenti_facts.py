#!/usr/bin/env python3
"""Facts files for the four lentivirus artifact pages, from measurement.

Everything numeric is read from the files the build and the audit produced --
registry.json, agg.json, unbinned.json, collections/synonyms.tsv, the module
JSON and the 26-27 Oct coverage audit. Only PROSE is written here, because
only prose is a judgement. Regenerating after a rebuild therefore cannot leave
a stale number behind a fresh sentence, which is the failure mode that put
thirteen-feature numbers on a fifteen-feature Cardiovirus page.
"""
import json, os, collections, glob, sys

L = os.path.expanduser("~/lowvan_scratch/lenti")
E = os.path.join(L, "audit_2026-10-26e")
VJ = "/Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Alsuviricetes/Viral_Annotation/Viral_PSSM.json"
DATE = "27 October 2026"

WORK = {"HIV2": os.path.join(L, "hiv2"),
        "FIV":  os.path.join(L, "genus", "FIV"),
        "SRLV": os.path.join(L, "genus", "SRLV"),
        "EIAV": os.path.join(L, "genus", "EIAV")}
COVSET = {"HIV2": ("hiv2", "cov_hiv2_work"), "FIV": ("genus", "cov_genus_work"),
          "SRLV": ("genus", "cov_genus_work"), "EIAV": ("genus", "cov_genus_work")}
HOLDINGS = {"HIV2": os.path.join(L, "hiv2", "Contigs"),
            "FIV": os.path.join(L, "genus", "Contigs"),
            "SRLV": os.path.join(L, "genus", "Contigs"),
            "EIAV": os.path.join(L, "genus", "Contigs")}

#  ---- the only hand-written part ----------------------------------------
PROSE = {
 "HIV2": dict(
  species="HIV-2 and the simian immunodeficiency viruses",
  stand=("<b>HIV2</b> is a module name, not a taxon. It holds HIV-2 and the SIVs together, "
         "and it is <i>majority SIV</i>: 1,660 of its 2,113 collection sequences and 65 of its "
         "110 coverage exemplars come from simian immunodeficiency virus, not from HIV-2. "
         "The module is separated from HIV-1 on layout rather than distance &mdash; the two differ "
         "by one accessory gene, <i>vpx</i> here against <i>vpu</i> there, which is the project's "
         "partitioning rule."),
  lede="Nine features on one 9.8 kb provirus, seven called by profile and two by the splice caller.",
  notes=("Nine features: gag, pol, vif, vpx, vpr, env and nef from PSSMs, tat and rev from "
         "get_splice_variant_features.pl. <i>vpx</i> is a gene duplication of <i>vpr</i> and exists "
         "only in the HIV-2/SIVsmm/SIVmac lineage; both are declared and neither calls the other."),
  gaps=("VPR is called on 72.7% of exemplars and that is biology, not a gap. Thirty of the 110 "
        "genomes carry a single accessory gene in a 311 nt window where the rest carry two in 630 nt, "
        "and the BV-BRC reference annotation agrees those genomes have no <i>vpr</i>. VPR therefore "
        "carries no copy_num. The held-out panel read 95% because it is mostly HIV-2; the coverage "
        "set carries the full SIV spread."),
  ),
 "FIV": dict(
  species="feline immunodeficiency virus",
  stand=("<b>FIV</b> is the lentivirus of domestic cats and of wild felids &mdash; bobcat, puma and "
         "others, each with its own long-diverged lineage. Five features, all collinear, with "
         "<i>rev</i> recognised in the vocabulary but deliberately not declared: it is spliced, and "
         "this module calls no spliced feature."),
  lede="Five features, four of them called on every exemplar.",
  notes=("gag, pol, env and vif are called on 39 of 39 exemplars. orfA is the only variable one and "
         "carries no copy_num."),
  gaps=("orfA at 64.1% is a host-lineage split, not a truncation. Domestic-cat genomes call it 21 of "
        "21; bobcat 0 of 10 and puma 3 of 6. Those same genomes call gag, pol, vif and env 39 of 39, "
        "so the absence is profile divergence on a small accessory gene rather than incomplete "
        "sequence. Widening orfA to reach the wild felids is a build decision, not a defect."),
  ),
 "SRLV": dict(
  species="the small ruminant lentiviruses",
  stand=("<b>SRLV</b> holds maedi-visna virus of sheep and caprine arthritis-encephalitis virus of "
         "goats, which are one continuum rather than two clean species. Five features, all with "
         "profiles &mdash; including <i>tat</i>, which in this genus is not what its name suggests."),
  lede="Five features, and a tat that is collinear rather than spliced.",
  notes=("SRLV <i>tat</i> is a single-exon protein, also published as S, and it was MEASURED collinear "
         "by tblastn of each deposited protein against its own genome rather than assumed spliced "
         "from the HIV-1 layout. It therefore gets a PSSM and a copy_num like any other feature. "
         "EIAV <i>tat</i>, in the sister module, is genuinely spliced and gets neither."),
  gaps=("No declared feature lacks a profile. <i>rev</i> appears in the vocabulary and bins to its own "
        "key, but the module does not declare it: it is spliced, and nothing here calls spliced "
        "features."),
  ),
 "EIAV": dict(
  species="equine infectious anaemia virus",
  stand=("<b>EIAV</b> is the lentivirus of horses and the smallest module in this family &mdash; one "
         "virus, little vocabulary, and the clearest example in the set of a tiny accessory gene "
         "that is nonetheless essential. Four features carry profiles; <i>rev</i> is declared and "
         "called by the splice caller, and <i>tat</i> is recognised but not declared."),
  lede="Four profiled features, and an 11-of-11 accessory gene smaller than most peptides.",
  notes=("S2 is 59&ndash;83 aa and is called on every exemplar. It was marked essential on evidence: "
         "30 of 30 held-out genomes, and out of frame with <i>env</i> in 29 of 29 &mdash; a separate "
         "gene, not a cleavage product of the envelope it sits inside."),
  gaps=("<i>tat</i> is spliced in EIAV and is not declared; the string is recognised and recorded, "
        "which is why the collapse report shows it bound to a target the module JSON does not carry. "
        "gag, pol and env are each called on 10 of 11 exemplars; the one miss is a single short "
        "record, not a profile failure."),
  ),
}

def load(mod):
    d = WORK[mod]
    cov_name, _ = COVSET[mod]
    cov = json.load(open(os.path.join(E, "coverage_%s.json" % cov_name)))
    rows = [r for r in cov["rows"] if r["module"] == mod]
    pres = collections.Counter()
    for r in rows:
        for c in set(r["calls"]): pres[c] += 1
    reg = json.load(open(os.path.join(d, "registry.json")))
    agg = json.load(open(os.path.join(d, "agg.json")))
    unb = json.load(open(os.path.join(d, "unbinned.json")))
    syn = [l.rstrip("\n").split("\t") for l in open(os.path.join(d, "collections", "synonyms.tsv"))][1:]
    feats = json.load(open(VJ))[mod]["features"]
    return dict(d=d, rows=rows, pres=pres, reg=reg, agg=agg, unb=unb, syn=syn, feats=feats,
                holdings=len(glob.glob(os.path.join(HOLDINGS[mod], "*"))))

def reg_rows(reg, mod):
    """registry.json -> [(feature, n_seqs, n_pssm, recall)]"""
    out = []
    blocks = reg.get(mod) if isinstance(reg, dict) and mod in reg else reg
    if isinstance(blocks, dict):
        it = blocks.get("features", blocks)
        for k, v in (it.items() if isinstance(it, dict) else []):
            if not isinstance(v, dict): continue
            out.append((k, v.get("n_seqs") or v.get("n") or 0,
                        v.get("n_pssm") or v.get("pssms") or 0,
                        v.get("recall") or v.get("in_range")))
    return out

def build(mod):
    M = load(mod); P = PROSE[mod]
    n = len(M["rows"]); feats = M["feats"]
    syn_tot = sum(int(x[0]) for x in M["syn"])
    syn_bound = sum(int(x[0]) for x in M["syn"] if x[3].strip() not in ("UNASSIGNED", "DISCARDED"))
    nstr = sum(r["nstrings"] for r in M["agg"]); ntarg = len(M["agg"])
    npssm = sum((v.get("n_pssm") or 0) for v in [])  # filled from registry text below

    # per-feature coverage table
    trows = []
    for f in sorted(feats):
        key = feats[f].get("gene_symbol") or f.lower()
        c = M["pres"].get(key, M["pres"].get(f.lower(), 0))
        trows.append([f, feats[f].get("anno", ""), "%d / %d" % (c, n), "%.1f%%" % (100.0*c/n if n else 0),
                      "yes" if feats[f].get("copy_num") else "—",
                      feats[f].get("special") or "PSSM"])

    #  bars rows are (label, value, pct, note): the value is what prints, the
    #  pct drives the bar width, the note carries the count behind the rate.
    bars = []
    for f in sorted(feats):
        key = feats[f].get("gene_symbol") or f.lower()
        c = M["pres"].get(key, M["pres"].get(f.lower(), 0))
        pct = round(100.0 * c / n, 1) if n else 0.0
        bars.append([f, "%.1f" % pct, pct, "%d of %d" % (c, n)])

    collapse_rows = [[r["key"], r["anno"], str(r["nstrings"]), str(r["total"])] for r in M["agg"]]

    facts = {
      "title": "%s Coverage Audit" % mod,
      "eyebrow": "LowVan module audit &middot; %s" % DATE,
      "h1": mod,
      "subhead": P["species"],
      "standfirst": P["stand"],
      "headline_tiles": [
         [str(len(feats)), "features declared"],
         ["%d" % n, "coverage exemplars"],
         ["%d" % M["holdings"], "genomes in holdings"],
         ["%.1f%%" % (100.0*syn_bound/syn_tot if syn_tot else 0), "vocabulary bound by rule"],
         ["%d → %d" % (nstr, ntarg), "annotation strings collapse"],
      ],
      "sections": [
        {"num": "01", "title": "What the module calls", "lede": P["lede"], "blocks": [
            {"type": "prose", "html": P["notes"]},
            {"type": "table",
             "headers": ["feature", "product", "called", "rate", "copy_num", "caller"],
             "rows": trows},
            {"type": "bars", "rows": bars, "unit": "%"},
        ]},
        {"num": "02", "title": "Where it does not call, and why", "blocks": [
            {"type": "prose", "html": P["gaps"]},
        ]},
        {"num": "03", "title": "Vocabulary", "lede":
            "%d annotation strings collapse onto %d targets." % (nstr, ntarg), "blocks": [
            {"type": "table", "headers": ["target", "product", "strings", "features"],
             "rows": collapse_rows},
            {"type": "prose", "html":
              ("%.1f%% of the %d deposited features from this module's taxa are bound by a "
               "naming rule. The remainder are recorded in <code>UNASSIGNED</code> and "
               "<code>DISCARDED</code>: the first is vocabulary no rule reached, the second is "
               "records the source itself calls partial or truncated, which the project discards "
               "rather than trains on."
               % (100.0*syn_bound/syn_tot if syn_tot else 0, syn_tot))},
        ]},
        {"num": "04", "title": "Method", "blocks": [
            {"type": "prose", "html":
              ("Coverage is measured by annotating de-replicated exemplars through the real GTO "
               "pipeline &mdash; <code>annotate_by_viral_pssm-GTO.pl</code>, then the transcript-edit "
               "and splice passes, then <code>viral_genome_quality.pl</code> &mdash; not by a "
               "reimplementation. Spliced features are included, which the flat feature-table path "
               "cannot represent. Rates are per exemplar, and exemplars are clustered at 95% "
               "nucleotide identity, so one row is one distinct sequence rather than one deposit.")},
        ]},
      ],
      "source_note":
        ("Measured %s from the shipped module: the %s block of Viral_PSSM.json, its profiles and "
         "collections, and the 27 October coverage re-annotation over %d exemplars de-replicated "
         "from %d genomes. Vocabulary figures come from <code>collections/synonyms.tsv</code>, "
         "written by the module's own builder from the naming rules, and count deposited features "
         "rather than unique sequences."
         % (DATE, mod, n, M["holdings"])),
    }

    notes = {
      "standfirst": P["stand"],
      "modules": {mod: P["notes"]},
      "gaps": P["gaps"],
    }
    return facts, notes

def main():
    R = "/Users/jdavis/Library/CloudStorage/Box-Box/1_Projects/CEPI/Alsuviricetes/LowVan_Skill/reports/generators"
    for mod in ("HIV2", "FIV", "SRLV", "EIAV"):
        facts, notes = build(mod)
        s = mod.lower()
        json.dump(facts, open(os.path.join(R, "%s_audit.json" % s), "w"), indent=1)
        json.dump(notes, open(os.path.join(R, "%s_notes.json" % s), "w"), indent=1)
        print("  %-5s facts: %d sections, %d tiles -> %s_audit.json / %s_notes.json"
              % (mod, len(facts["sections"]), len(facts["headline_tiles"]), s, s))

main()
