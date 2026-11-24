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
    ap.add_argument("--contigs", default=None,
                    help="directory of <genome>.fna to measure the segment "
                         "length range from (default: <workdir>/Contigs)")
    ap.add_argument("--segments", default=None,
                    help='JSON: {"RNA1": {"min_len":.., "max_len":..}, ...}')
    ap.add_argument("--allow-dropped", action="store_true",
                    help="permit a feature declared `special` in the existing "
                         "--out JSON to disappear from it. Use only when that "
                         "feature is genuinely being retired.")
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
        #  The derived path below sets replicon_geometry; this one did not, so
        #  every module built with --segments shipped a segments block that
        #  install_module.py reports as incomplete. Measured 29 Oct 2026: 187
        #  segments in the runtime repo carry "linear" and the only 9 missing
        #  it were the nine Parvoviridae modules, i.e. every module that took
        #  this branch. Default it here rather than making each caller
        #  remember, and let a caller that knows better override it.
        for _s in segments.values():
            if isinstance(_s, dict):
                _s.setdefault("replicon_geometry", "linear")
    else:
        #  --contigs, because the contigs are not always under the workdir.
        #  Every Retroviridae module keeps them in a sibling
        #  Contigs.<Genus>/ instead, and the hard-wired <workdir>/Contigs
        #  then matched nothing: L stayed empty and the next line raised
        #  IndexError on L[0], which reads as a corrupt file rather than a
        #  path that is simply somewhere else.
        cd = a.contigs or os.path.join(W, "Contigs")
        L = []
        for f in glob.glob(os.path.join(cd, "*.fna")):
            n = sum(len(x.strip()) for x in open(f) if not x.startswith(">"))
            if n: L.append(n)
        if not L:
            sys.exit("no *.fna with sequence under %s -- pass --contigs <dir> "
                     "or --segments '<json>'" % cd)
        L.sort()
        segments = {"Single RNA Segment": {
            "min_len": int(L[0] * 0.95), "max_len": int(L[-1] * 1.05),
            "replicon_geometry": "linear"}}

    #  THE FEATURES AND THE SEGMENTS BLOCK MUST AGREE, and nothing used to
    #  check that they did. The derived branch above always writes exactly one
    #  segment called "Single RNA Segment", because it measures contig files
    #  and cannot know a genome is segmented. A segmented module whose features
    #  each declare their own segment therefore got a block naming a segment no
    #  feature belongs to, and no feature's segment declared at all.
    #
    #  That is the malformed-segments failure this file's own docstring opens
    #  with: annotation succeeds, every feature is called, and quality scoring
    #  is dead on every genome, because the scorer looks each feature's segment
    #  up in this block and finds nothing.
    #
    #  Measured on Birnaviridae, 19 Nov 2026: 11 features declared `Segment A`
    #  and 1 declared `Segment B`, and the generated block declared
    #  `Single RNA Segment` with a 4692-7203 nt range taken from the merged
    #  A+B genomes. It printed `1 segment(s)` and exited 0.
    #
    #  Refuse instead, and say which segments the features actually want, so
    #  the caller passes --segments with their measured ranges.
    declared = {e.get("segment") for e in feats.values() if e.get("segment")}
    unknown = sorted(d for d in declared if d not in segments)
    if unknown:
        sys.exit(
            "REFUSED: %d feature segment name(s) are not in the segments "
            "block: %s\n"
            "  the block declares: %s\n"
            "  Annotation would still succeed and quality scoring would be "
            "dead on every genome -- the scorer resolves a feature's segment "
            "through this block.\n"
            "  Pass --segments with one entry per segment, e.g.\n"
            "    --segments '%s'"
            % (len(unknown), ", ".join(unknown), ", ".join(sorted(segments)),
               json.dumps({d: {"min_len": 0, "max_len": 0} for d in unknown})))

    out = {}
    missing = []
    for key in sorted(feats):
        spec = feats[key]
        fa = os.path.join(cdir, key + ".fasta")
        #  A `special` FEATURE HAS NO COLLECTION, BY DESIGN. It is called by
        #  get_splice_variant_features.pl or the transcript-edit pass from the
        #  references in Splice-Variants/<Module>/<FEAT>.fasta, and never gets a
        #  PSSM, so there is nothing in collections/ for it. Requiring a
        #  collection here silently DROPPED it: Betaretrovirus REM was declared
        #  in features_all.json, had 7 verified splice references on disk, and
        #  came out of this script as a 4-feature module with one line of
        #  "NO COLLECTION (not declared)" among the normal output.
        #
        #  Its length range cannot come from a collection either, so the spec
        #  must carry min_len/max_len and it is an error if it does not.
        if spec.get("special") and not os.path.exists(fa):
            if "min_len" not in spec or "max_len" not in spec:
                sys.exit("feature %s is special:%s with no collection, so it must "
                         "declare min_len and max_len in the features file"
                         % (key, spec["special"]))
            lo, hi = int(spec["min_len"]), int(spec["max_len"])
            med = (lo + hi) // 2
        else:
            #  A MERGED FEATURE HAS NO COLLECTION OF ITS OWN EITHER, for a
            #  different reason than `special`: it was built from two or more
            #  collections by merge_feature_profiles.py, because its endpoint
            #  is decided by something the clustering cannot see.
            #  Orthohepadnavirus P17 is built from P17_NONA and P17_A, the
            #  furin site being genotype-dependent.
            #
            #  Those sources are deliberately NOT concatenated into a
            #  P17.fasta. That would add a collection holding a second copy of
            #  every sequence in two others -- which the bin audit would
            #  rightly call a duplicate bin -- and it would invalidate both
            #  audits for a file existing only to carry a length range.
            #
            #  So take the range from the union of the sources named in
            #  MERGED_FROM. Without this the feature lands in `missing` and is
            #  dropped SILENTLY, which is how Betaretrovirus REM shipped as a
            #  4-feature module (see the comment above).
            if not os.path.exists(fa):
                mf = os.path.join(W, "Alignments", M, key, "MERGED_FROM")
                srcs = []
                if os.path.exists(mf):
                    txt = "".join(l for l in open(mf) if not l.startswith("#"))
                    try:
                        srcs = sorted({r["from"] for r in json.loads(txt)})
                    except Exception:
                        srcs = []
                if not srcs:
                    missing.append(key); continue
                L = sorted(x for s in srcs
                           for x in lens(os.path.join(cdir, s + ".fasta")))
                print("  %-16s no collection; length range from the merged "
                      "sources %s" % (key, ", ".join(srcs)))
            else:
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
            #  0.65 is the default, not the only value. Since the 26 Oct
            #  caller, coverage is measured against the PROFILE's length
            #  rather than the HSP's own query span, so this threshold
            #  finally does something: it is what rejects a profile that
            #  matches part of an open reading frame and is then extended
            #  over the whole of it. The value was hardcoded here, so a
            #  features file that set one was silently ignored.
            "coverage_cutoff": spec.get("coverage_cutoff", 0.65),
            #  Scan to the stop unless the curator says otherwise. Defaulting
            #  this to 0 made every adenovirus CDS end where its tBLASTn HSP
            #  ended (15 Nov 2026: protease short in 2,447/3,777 genomes). A
            #  mat_peptide still states both flags from its own termini.
            "downstream_ext": spec.get("downstream_ext",
                                       0 if ftype == "mat_peptide" else 1),
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
        #  Fields that change HOW a feature is called, not merely how it is
        #  described. Each was being dropped here, which is silent and total:
        #  a feature declared `special: transcript_edit` in the spec arrived in
        #  the JSON as an ordinary PSSM feature, so the annotator looked for a
        #  profile that by definition does not exist and the protein was never
        #  called at all. Astrovirus ORF1ab is reached by a -1 ribosomal
        #  frameshift and is the first feature in this project to need it.
        #
        #    special            transcript_edit | splice -- called by an
        #                       external program from a nucleotide reference
        #    internal_stop      the match may legitimately contain a stop
        #    translate          emit a protein translation for this feature
        #    non_pssm_partner   this feature's coordinates come from a
        #                       neighbour's PSSM match
        #    begin / end        the anchor pair that does it
        #    require_met_start  veto a hit whose first codon is not ATG.
        #                       This one is NOT cosmetic and was being dropped
        #                       with the rest: it is the ONLY thing separating
        #                       a nested in-frame isoform from the shorter
        #                       form it contains, where the sole difference is
        #                       an upstream ATG and a PSSM cannot veto on one
        #                       column. Orthohepadnavirus WHOLEX carried it in
        #                       the shipped JSON, set by hand; the next
        #                       regeneration would have removed it and the
        #                       HBwx profile would have called plain X in
        #                       every X-only genome. WHOLEL has exactly the
        #                       same relationship to L.
        for field in ("special", "internal_stop", "translate",
                      "non_pssm_partner", "begin", "end",
                      "require_met_start"):
            if spec.get(field) is not None:
                d[field] = spec[field]
        #  A `special` feature has no PSSM, so the fields that only describe
        #  how a PSSM match is scored and extended describe nothing. Writing
        #  them is not harmless: a bit_cutoff on a feature with no profile
        #  reads as a threshold someone chose and invites a future curator to
        #  tune it. 51 of the 54 special features already in the repo omit
        #  them; this makes that the rule rather than the habit.
        if d.get("special"):
            for field in ("bit_cutoff", "coverage_cutoff", "kmers",
                          "upstream_ext", "downstream_ext"):
                d.pop(field, None)
        #  Citations. PMID is the curator's: a reference a person read and
        #  chose. PMID_claude_generated is this project's own proposal and is
        #  kept in a separate field precisely so the two are never confused.
        #  Both were being dropped on the floor here, the same way the
        #  extension flags were, so a feature declared with a literature
        #  reference shipped without one.
        #
        #  THE SHAPE MATTERS, AND THIS LINE IS WHERE IT WAS DESTROYED. Two
        #  value forms are in use and both are correct: a bare list of ids,
        #  and a {"function": [...], "coordinates": [...]} object that says
        #  which claim each citation supports. `[str(v)]` turned the object
        #  form into a PYTHON REPR INSIDE A JSON STRING --
        #  ["{'function': ['399327']}"] -- which nothing reads as a citation.
        #  audit_pmids.py found twelve of these in the two hepadnavirus
        #  modules on 24 Nov 2026 and repaired them; they came straight back
        #  on the next regeneration, because the repair fixed the file and
        #  this line is what writes the file. Only a scalar needs wrapping.
        for field in ("PMID", "PMID_claude_generated"):
            if spec.get(field):
                v = spec[field]
                d[field] = v if isinstance(v, (list, dict)) else [str(v)]
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

    #  A FEATURE THIS SCRIPT ALREADY SHIPPED MUST NOT VANISH WITHOUT SAYING SO.
    #
    #  This JSON is regenerated from the features file, so anything the
    #  features file does not list is deleted from the module -- silently,
    #  because a feature that is gone produces no output of any kind. The way
    #  it happens is specific and it has happened: `special` is a curation
    #  decision that lived ONLY in this JSON, the per-taxon setup script read
    #  it from here to EXCLUDE the feature from the features file (correctly --
    #  a special must not get profiles), and the next regeneration then had no
    #  record of it. All four Adenoviridae modules lost K33 (33K, spliced,
    #  PMID 25636034) that way on 5 Nov 2026, with 578/180/38/40 validated
    #  splice references still on disk. Nothing reported a gap: there is no
    #  profile to be stale, no collection to be empty and no call rate to
    #  fall, because the feature simply stopped being asked for.
    #
    #  So diff against what is already there. A retirement is a real decision
    #  and is allowed -- it just has to be stated, not discovered later.
    if os.path.exists(a.out):
        try:
            prev = (json.load(open(a.out)).get(M) or {}).get("features") or {}
        except Exception:
            prev = {}
        gone = sorted(set(prev) - set(out))
        if gone:
            lost = [k for k in gone if prev[k].get("special")]
            print("  DROPPED from the previous JSON: %s" % ", ".join(gone))
            if lost and not a.allow_dropped:
                sys.exit(
                    "REFUSED: %d feature(s) declared `special` in %s are not in\n"
                    "  %s, so this would delete them:\n    %s\n"
                    "  A special has no collection and no profile, so nothing\n"
                    "  downstream will notice it is gone. Declare it in the\n"
                    "  features file (carrying its `special` flag), or pass\n"
                    "  --allow-dropped if the feature is genuinely retired."
                    % (len(lost), a.out, a.features, ", ".join(lost)))

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
