#!/usr/bin/env python3
"""Write a module's Viral_PSSM.json block from what was actually built.

Earlier modules were assembled by hand, which is how Pestiviridae shipped
with segments written as {"Single RNA Segment": 1} -- annotation worked and
quality scoring was dead on every genome. This derives every number from an
artifact on disk instead:

  min_len / max_len   the built collection's own length range, IN AMINO ACIDS
  close_genomes       Rep-Contigs/close_genomes.json, as written by
                      build_rep_contigs.py
  segments            measured contig lengths, in the dict form the quality
                      script needs

COPY_NUM IS NEVER INVENTED HERE, and that is deliberate. copy_num is what marks a
feature essential, so a feature that carries it flags every genome that lacks
it. Essentiality cannot be known before the module has been run against
genomes -- it is called/routed >= 0.85, not "is this gene famous". Getting it
from the gene's reputation put three shipped modules at 0.0% clean. So the
module ships without copy_num, run_gto_eval.py measures called/routed, and
copy_num is added afterwards to the features that earn it. A missing copy_num
costs nothing; a wrong one fails every genome in the taxon.

A copy_num ALREADY in the feature spec is a different thing: it is that
measurement, folded back by the curator after apply_copy_num.py. It is passed
through. Dropping it was a silent data loss -- apply_copy_num wrote the flags,
the next routine regeneration of the JSON wiped them, and the module shipped
marking nothing essential. Every Caliciviridae module did exactly that.

  python3 gen_module_json.py --workdir Crinivirus --module Crinivirus \
      --features Crinivirus/features.json --out Crinivirus/Crinivirus_Viral_PSSM.json
"""
import argparse, collections, glob, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from json_canon import write_json   # one canonical format, or diffs are noise

#  bit_cutoff convention, read off the shipped closterovirid and Bromoviridae
#  modules. Role first, then length for anything unnamed.
ROLE = {"POLY": 400, "REP": 400, "REP1": 400, "REP2": 300, "RDRP": 300,
        "HSP70": 150, "HSP90": 150, "CP": 80, "CPM": 80, "MP": 80}


def cutoff(key, medlen):
    if key in ROLE:
        return ROLE[key]
    if medlen >= 250:
        return 60
    if medlen >= 80:
        return 30
    return 25


def lens(path):
    out, cur = [], 0
    for l in open(path):
        if l.startswith(">"):
            if cur: out.append(cur)
            cur = 0
        else:
            cur += len(l.strip())
    if cur: out.append(cur)
    return sorted(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--module", required=True)
    ap.add_argument("--features", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--segments", default=None,
                    help='JSON: {"RNA1": {"min_len":.., "max_len":..}, ...}')
    ap.add_argument("--pad", type=float, default=0.06,
                    help="fractional slack on the collection length range")
    a = ap.parse_args()

    W, M = a.workdir, a.module
    feats = json.load(open(a.features))
    cdir = os.path.join(W, "collections", M)

    cg_path = os.path.join(W, "Rep-Contigs", "close_genomes.json")
    if not os.path.exists(cg_path):
        sys.exit("no %s -- run build_rep_contigs.py first" % cg_path)
    close = json.load(open(cg_path))
    #  Never declare a reference whose .dna is gone. close_genomes.json is
    #  written once by build_rep_contigs.py and then edited by hand when a
    #  reference turns out to be wrong -- deleting the file is the obvious
    #  move and updating the manifest is the easy thing to forget. Vesivirus.19
    #  was "Bovine calicivirus Mukti/India", a NEBOVIRUS misfiled in the
    #  vesivirus reference set; it was removed from Rep-Contigs but left in the
    #  manifest, so the next regeneration silently re-declared it and the
    #  module shipped pointing at a BLASTn subject that does not exist.
    orphaned = sorted(k for k in close
                      if not os.path.exists(os.path.join(W, "Rep-Contigs", k)))
    for k in orphaned:
        del close[k]
    if orphaned:
        print("  dropped %d close_genome(s) with no .dna in Rep-Contigs: %s"
              % (len(orphaned), ", ".join(orphaned)))

    if a.segments:
        segments = json.load(open(a.segments)) if os.path.exists(a.segments) \
                   else json.loads(a.segments)
    else:
        L = []
        for f in glob.glob(os.path.join(W, "Contigs", "*.fna")):
            n = sum(len(x.strip()) for x in open(f) if not x.startswith(">"))
            if n: L.append(n)
        L.sort()
        segments = {"Single RNA Segment": {
            "min_len": int(L[0] * 0.95), "max_len": int(L[-1] * 1.05),
            "replicon_geometry": "linear"}}

    out = {}
    missing = []
    for key in sorted(feats):
        spec = feats[key]
        fa = os.path.join(cdir, key + ".fasta")
        if not os.path.exists(fa):
            missing.append(key); continue
        L = lens(fa)
        if not L:
            missing.append(key); continue
        lo = int(L[0] * (1 - a.pad)); hi = int(L[-1] * (1 + a.pad))
        med = L[len(L) // 2]
        seg = spec.get("segment") or (list(segments)[0] if len(segments) == 1 else None)
        if seg is None:
            sys.exit("feature %s has no 'segment' and the module has %d segments"
                     % (key, len(segments)))
        #  upstream_ext=1 makes the annotator run scan_to_met_start, which walks
        #  LEFT in 3-nt steps until it meets a Met or a stop. An internal mature
        #  peptide never begins with Met, so the scan fires on every call, and
        #  inside a polyprotein ORF there is no stop to halt it -- it runs into
        #  the upstream neighbour until it reaches an arbitrary Met. Norovirus
        #  shipped that way and every one of the five junctions overlapped on
        #  40/40 panel genomes (NS5/NS6 by a median 496 nt); the six peptides
        #  summed to a median 371 aa MORE than the polyprotein containing them.
        #  So a mat_peptide defaults to no N-terminal scan. It is still only a
        #  default: a product whose N-terminus really is the parent's start Met
        #  may set upstream_ext back to 1, and for that product the scan is a
        #  no-op anyway because the hit already begins with M.
        ftype = spec.get("feature_type", "CDS")
        d = {
            "anno": spec["anno"],
            "bit_cutoff": spec.get("bit_cutoff", cutoff(key, med)),
            "coverage_cutoff": 0.65,
            "downstream_ext": spec.get("downstream_ext", 0),
            "feature_type": ftype,
            #  The key is filename-safe and the symbol is what the annotator
            #  and the genome-quality tool report, so they are not always the
            #  same string -- Enterovirus VPG/VPg, Cardiovirus LSTAR/L*, and
            #  any polymerase, which this project names RDRP whatever its
            #  positional label in the taxon. Honour an explicit gene_symbol;
            #  it used to be silently overwritten with the key, so Orbivirus
            #  VP1 shipped as symbol "VP1" while carrying the SAME annotation
            #  string as Rotavirus VP1, which ships as RDRP.
            "gene_symbol": spec.get("gene_symbol", key),
            "kmers": spec.get("kmers", 1),
            "max_len": hi,
            "min_len": lo,
            "segment": seg,
            "type": ftype,
            "upstream_ext": spec.get("upstream_ext",
                                     0 if ftype == "mat_peptide" else 1),
        }
        if spec.get("comment"): d["comment"] = spec["comment"]
        #  Citations. PMID is the curator's: a reference a person read and
        #  chose. PMID_claude_generated is this project's own proposal and is
        #  kept in a separate field precisely so the two are never confused.
        #  Both were being dropped on the floor here, the same way the
        #  extension flags were, so a feature declared with a literature
        #  reference shipped without one.
        for field in ("PMID", "PMID_claude_generated"):
            if spec.get(field):
                v = spec[field]
                d[field] = v if isinstance(v, list) else [str(v)]
        #  copy_num is still never INVENTED here -- see the note at the top --
        #  but one already in the spec is a MEASURED value the curator folded
        #  back after apply_copy_num.py, and silently dropping it destroyed the
        #  measurement. Every Caliciviridae module lost its essentiality flags
        #  this way: apply_copy_num wrote them, the next routine regeneration
        #  wiped them, and the module shipped flagging nothing. Pass through
        #  what is explicitly there; invent nothing.
        if spec.get("copy_num") is not None:
            d["copy_num"] = spec["copy_num"]
        out[key] = d

    block = {M: {"close_genomes": close, "features": out, "segments": segments}}
    #  Canonical on the way out. json.dump(indent=1) is not the project
    #  format, and a module JSON that is not canonical makes its first real
    #  diff unreadable -- Perl rewrites it in hash order the moment any tool
    #  touches it.
    write_json(block, a.out)
    print("  %s: %d features, %d references, %d segment(s) -> %s"
          % (M, len(out), len(close), len(segments), a.out))
    if missing:
        print("  NO COLLECTION (not declared): %s" % ", ".join(missing))
    ncn = sum(1 for v in out.values() if "copy_num" in v)
    if ncn:
        print("  copy_num carried through on %d feature(s) from the spec "
              "(measured, not invented)" % ncn)
    else:
        print("  copy_num set on 0 features by design -- run run_gto_eval.py, "
              "then add copy_num where called/routed >= 0.85")


main()
