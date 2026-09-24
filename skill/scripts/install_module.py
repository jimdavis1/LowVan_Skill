#!/usr/bin/env python3
"""Install a built LowVan module into a Viral_Annotation checkout.

Reads a module working directory and writes the four things the annotator
actually loads:

    Viral_PSSM.json                       <- module block merged in
    Viral-PSSMs/<Module>.pssms/<FEAT>/    <- the PSSMs
    Viral-Rep-Contigs/<Module>.<n>.dna    <- BLASTn routing subjects
    PSSM-Alignments/<Module>/<FEAT>/      <- alignments (housekeeping; unused at runtime)

Expected working-directory layout:

    <workdir>/
      *_Viral_PSSM.json                  one or more module blocks
      Alignments/<Module>/<FEAT>/pssms/<Module>.<FEAT>.<n>.pssm
      Alignments/<Module>/<FEAT>/corrected_alis/<n>.fa
      Rep-Contigs/<Module>.<n>.dna

Validation runs first and refuses to install a module whose JSON and PSSMs
disagree, because a feature declared in the JSON with no PSSM behind it is
silently never called, and a PSSM with no JSON entry is silently never loaded.

    python3 install_module.py --workdir . --repo ~/Viral_Annotation
    python3 install_module.py --workdir . --repo ~/Viral_Annotation --check
"""

import argparse
import glob
import json
import os
import re
import shutil
import sys
from collections import OrderedDict

from json_canon import dumps as canon_dumps


def find_module_json(workdir):
    cands = [p for p in glob.glob(os.path.join(workdir, "*_Viral_PSSM.json"))
             if os.path.basename(p) != "Viral_PSSM.json"]
    if not cands:
        sys.exit("no *_Viral_PSSM.json found in %s" % workdir)
    if len(cands) > 1:
        sys.exit("several module JSONs found; pass --json: %s" % ", ".join(cands))
    return cands[0]


def scan_pssms(workdir, module):
    """{feature_key: [pssm paths]} for one module."""
    out = {}
    base = os.path.join(workdir, "Alignments", module)
    if not os.path.isdir(base):
        return out
    for feat in sorted(os.listdir(base)):
        pdir = os.path.join(base, feat, "pssms")
        if os.path.isdir(pdir):
            hits = sorted(glob.glob(os.path.join(pdir, "*.pssm")))
            if hits:
                out[feat] = hits
    return out


def validate(workdir, mod_json):
    """Return {module: pssms} for installable modules, printing every problem."""
    ok, problems = {}, []
    for module, block in mod_json.items():
        feats = set(block.get("features", {}))
        pssms = scan_pssms(workdir, module)
        built = set(pssms)

        # A feature declared with "special" is called by an external program
        # (transcript_edit / splice), so it legitimately has no PSSM.
        special = {f for f in feats if block["features"][f].get("special")}
        #  A feature that takes its coordinates from a NEIGHBOUR'S PSSM match
        #  has no profiles of its own and is not missing. annotate_by_viral_pssm.pl
        #  calls these through non_pssm_partner: the partner declares the
        #  partnership and the derived feature declares `begin`/`end` with a
        #  pssm name, a START/STOP anchor and an offset. Pestiviridae NS2 is the
        #  first of them -- it is cut from p7's C-terminus to the residue before
        #  NS3's N-terminal Gly, the Gly1590 of the conserved Arg1589-Gly1590
        #  site -- and without this the installer reported it as "DECLARED BUT
        #  UNBUILT -> would never be called", which is the opposite of true.
        #  The segments block must be {name: {min_len, max_len, replicon_geometry}}.
        #  viral_genome_quality.pl dereferences it as a hash, so a bare scalar
        #  kills quality scoring on EVERY genome with
        #      Can't use string ("1") as a HASH ref ... line 302
        #  while annotation itself succeeds -- so the module looks built and
        #  scores nothing. --check passed this shape until now.
        segs = block.get("segments")
        if not isinstance(segs, dict) or not segs:
            print("       segments block is missing or not an object"); problems += 1
        else:
            for sname, sv in segs.items():
                if not isinstance(sv, dict):
                    print("       segment %r is %r, expected an object with "
                          "min_len/max_len/replicon_geometry" % (sname, sv)); problems += 1
                    continue
                for need in ("min_len", "max_len", "replicon_geometry"):
                    if need not in sv:
                        print("       segment %r has no %s" % (sname, need)); problems += 1
            declared = {v.get("segment") for v in block["features"].values() if v.get("segment")}
            for d0 in sorted(declared - set(segs)):
                print("       feature(s) declare segment %r, which the segments "
                      "block does not define" % d0); problems += 1

        derived = {f for f in feats
                   if block["features"][f].get("begin") or block["features"][f].get("end")}
        missing = feats - built - special - derived
        orphan = built - feats

        n = sum(len(v) for v in pssms.values())
        print("  %-24s %3d features declared, %3d with PSSMs (%d files)"
              % (module, len(feats), len(built), n))
        if special:
            print("       %2d special (no PSSM expected): %s"
                  % (len(special), ", ".join(sorted(special))))
        if derived:
            print("       %2d derived from a neighbour's PSSM (no PSSM expected): %s"
                  % (len(derived), ", ".join(sorted(derived))))
            for f in sorted(derived):
                v = block["features"][f]
                for side in ("begin", "end"):
                    d = v.get(side) or {}
                    anc = d.get("%s_pssm" % side)
                    if anc and anc not in feats:
                        print("          %s: %s anchors on %r, which this module "
                              "does not declare" % (f, side, anc)); problems += 1
                    elif anc and f not in (block["features"][anc].get("non_pssm_partner") or []):
                        print("          %s: %s anchors on %s, but %s does not list it "
                              "in non_pssm_partner" % (f, side, anc, anc)); problems += 1
        if missing:
            print("       %2d DECLARED BUT UNBUILT -> would never be called:" % len(missing))
            print("          %s" % ", ".join(sorted(missing)))
            problems.append((module, "unbuilt", sorted(missing)))
        if orphan:
            print("       %2d PSSMs WITH NO JSON ENTRY -> would never be loaded:" % len(orphan))
            print("          %s" % ", ".join(sorted(orphan)))
            problems.append((module, "orphan", sorted(orphan)))

        # rep contigs
        rc = sorted(glob.glob(os.path.join(workdir, "Rep-Contigs", module + ".*.dna")))
        declared = set(block.get("close_genomes", {}))
        have = set(os.path.basename(p) for p in rc)
        if declared - have:
            print("       DECLARED BUT ABSENT rep contigs: %s" % ", ".join(sorted(declared - have)))
            problems.append((module, "repcontig", sorted(declared - have)))
        # The other direction matters just as much: a .dna file that nothing
        # declares still gets BLASTed (the annotator globs the directory), so it
        # routes genomes to this module while carrying no genome_name -- and a
        # regeneration of the JSON will not know it exists. This is how the
        # Rice yellow stunt rep contig went missing from close_genomes while
        # staying installed.
        if have - declared:
            print("       PRESENT BUT UNDECLARED rep contigs: %s" % ", ".join(sorted(have - declared)))
            problems.append((module, "repcontig-undeclared", sorted(have - declared)))
        thin = [k for k, v in block.get("close_genomes", {}).items()
                if not (v.get("genome_ids") or "").strip()
                or not (v.get("genome_name") or "").strip()]
        if thin:
            print("       rep contigs with an empty genome_ids or genome_name: %s"
                  % ", ".join(sorted(thin)))
            problems.append((module, "repcontig-blank", sorted(thin)))
        if not built:
            print("       SKIPPING: no PSSMs at all")
            continue
        ok[module] = pssms
    return ok, problems



def _retire(parent, keep, what, module, dry):
    """Remove feature subdirectories under parent that are not in keep."""
    if not os.path.isdir(parent):
        return
    gone = sorted(d for d in os.listdir(parent)
                  if os.path.isdir(os.path.join(parent, d)) and d not in keep)
    for d in gone:
        print("  %-24s retiring %s/%s (%s no longer built)"
              % (module, os.path.basename(parent), d, what))
        if not dry:
            shutil.rmtree(os.path.join(parent, d), ignore_errors=True)


def install(workdir, repo, modules, mod_json, dry=False):
    vp = os.path.join(repo, "Viral_PSSM.json")
    if not os.path.exists(vp):
        sys.exit("not a Viral_Annotation checkout (no Viral_PSSM.json): %s" % repo)

    for module, pssms in modules.items():
        # 1. PSSMs -> Viral-PSSMs/<Module>.pssms/<FEAT>/
        dest = os.path.join(repo, "Viral-PSSMs", module + ".pssms")
        #  Prune features this build no longer produces, BEFORE copying the
        #  ones it does. Installing without pruning leaves a retired feature's
        #  profiles live in the repo, where the annotator will happily go on
        #  calling them: an Alphagymnorhavirus nucleocapsid kept being emitted
        #  from a profile whose feature had no sequences left anywhere in the
        #  working tree, and the UNCHAR directories outlived the UNC groups
        #  that replaced them. A profile the current build cannot reproduce
        #  must not survive an install of that build.
        _retire(dest, set(pssms), "profiles", module, dry)
        _retire(os.path.join(repo, "PSSM-Alignments", module), set(pssms),
                "alignments", module, dry)
        for feat, files in pssms.items():
            fd = os.path.join(dest, feat)
            if not dry:
                os.makedirs(fd, exist_ok=True)
            #  Prune stale PROFILES inside a feature that is still built.
            #  _retire above only removes whole retired features. A cluster that
            #  was renumbered or dropped leaves its .pssm behind in a directory
            #  that still exists, and the annotator globs that directory, so the
            #  orphan goes on being used. A stale Togaviridae.P123.6.pssm from an
            #  abandoned first build survived exactly this way -- its cluster 1
            #  began 190 residues in and placed P123 122 residues late on
            #  Chikungunya. Nothing reported it; the installed count simply did
            #  not match the built count.
            want = {os.path.basename(f) for f in files}
            if os.path.isdir(fd):
                for old_file in sorted(os.listdir(fd)):
                    if old_file.endswith(".pssm") and old_file not in want:
                        print("  %-24s retiring %s/%s (profile no longer built)"
                              % (module, feat, old_file))
                        if not dry:
                            os.remove(os.path.join(fd, old_file))
            if not dry:
                for f in files:
                    shutil.copy2(f, os.path.join(fd, os.path.basename(f)))
        n = sum(len(v) for v in pssms.values())
        print("  %-24s %3d PSSMs -> Viral-PSSMs/%s.pssms/" % (module, n, module))

        # 2. alignments -> PSSM-Alignments/<Module>/<FEAT>/
        na = 0
        for feat in pssms:
            #  reclustered_alis holds the N-terminal splits, which are live
            #  profiles exactly like corrected_alis. Copying only the latter
            #  shipped Trivirinae with 6 of its 153 PSSMs having no
            #  alignment behind them -- unauditable, and unrebuildable, since
            #  rebuild_pssms.py refreshes an existing profile from its
            #  alignment and cannot recreate one that has none.
            fd = os.path.join(repo, "PSSM-Alignments", module, feat)
            srcs = []
            for sub in ("corrected_alis", "reclustered_alis"):
                src = os.path.join(workdir, "Alignments", module, feat, sub)
                if os.path.isdir(src):
                    srcs.extend(sorted(glob.glob(os.path.join(src, "*.fa"))))
            if srcs and not dry:
                os.makedirs(fd, exist_ok=True)
            #  Prune stale ALIGNMENTS inside a feature that is still built --
            #  the mirror of the profile prune above, which was missing. The
            #  copy below never clears the destination, so every alignment a
            #  previous build wrote survives a rebuild that renumbers or drops
            #  its cluster. Alphaflexiviridae shipped 276 profiles against 289
            #  alignments this way: 13 files from two earlier builds, including
            #  an NABP cluster 14 whose members had since fallen to the leftover
            #  pool. Alignments are provenance, not runtime, so nothing was
            #  miscalled -- but the counts no longer described the module, and
            #  the same gap is what made Hepeviridae read 15 profiles / 13
            #  alignments. Provenance that outlives what it documents is worse
            #  than absent: it is read as current.
            want_ali = {os.path.basename(f) for f in srcs}
            if os.path.isdir(fd):
                for old_file in sorted(os.listdir(fd)):
                    if old_file.endswith(".fa") and old_file not in want_ali:
                        print("  %-24s retiring %s/%s (alignment no longer built)"
                              % (module, feat, old_file))
                        if not dry:
                            os.remove(os.path.join(fd, old_file))
            for f in srcs:
                if not dry:
                    shutil.copy2(f, os.path.join(fd, os.path.basename(f)))
                na += 1
        print("  %-24s %3d alignments -> PSSM-Alignments/%s/" % ("", na, module))

        # 2b. special-feature references
        #
        #  A feature declaring `special: splice` or `special: transcript_edit`
        #  is not called by a PSSM. It is called by
        #  get_splice_variant_features.pl or get_transcript_edited_features.pl
        #  from hand-curated nucleotide references under
        #  <Splice-Variants|Transcript-Editing>/<Module>/<FEAT>.fasta.
        #
        #  Installing without them leaves the feature declared and unable to
        #  fire, and nothing reports it: the programs warn on stderr and exit
        #  0. Merhavirus L_SPLICED sat in exactly that state in
        #  Alpharhabdovirinae for the life of that module.
        for kind in ("Splice-Variants", "Transcript-Editing"):
            src = os.path.join(workdir, kind)
            #  accept both the module-scoped layout used under modules/ and
            #  the runtime layout <kind>/<Module>/
            for cand in (os.path.join(src, module), src):
                if os.path.isdir(cand) and glob.glob(os.path.join(cand, "*.fasta")):
                    dd = os.path.join(repo, kind, module)
                    n = 0
                    for f in sorted(glob.glob(os.path.join(cand, "*.fasta"))):
                        if not dry:
                            os.makedirs(dd, exist_ok=True)
                            shutil.copy2(f, os.path.join(dd, os.path.basename(f)))
                        n += 1
                    print("  %-24s %3d %s reference file(s) -> %s/%s/"
                          % ("", n, kind.split("-")[0].lower(), kind, module))
                    break

        # 3. rep contigs
        nr = 0
        for f in sorted(glob.glob(os.path.join(workdir, "Rep-Contigs", module + ".*.dna"))):
            if not dry:
                os.makedirs(os.path.join(repo, "Viral-Rep-Contigs"), exist_ok=True)
                shutil.copy2(f, os.path.join(repo, "Viral-Rep-Contigs", os.path.basename(f)))
            nr += 1
        print("  %-24s %3d rep contigs -> Viral-Rep-Contigs/" % ("", nr))

    # 4. merge the JSON and write it canonically
    with open(vp) as fh:
        before = fh.read()
    master = json.loads(before, object_pairs_hook=OrderedDict)
    added, replaced = [], []
    for module in modules:
        (replaced if module in master else added).append(module)
        master[module] = mod_json[module]
    after = canon_dumps(master)
    if not dry:
        shutil.copy2(vp, vp + ".bak")
        with open(vp, "w", encoding="utf8") as fh:
            fh.write(after)
    print("  Viral_PSSM.json: %d module(s) added %s, %d replaced %s (backup .bak)"
          % (len(added), added or "", len(replaced), replaced or ""))
    # The master is shared. Say plainly when a write reformats the whole file,
    # because that happens exactly once and looks alarming in a diff.
    if before != canon_dumps(json.loads(before, object_pairs_hook=OrderedDict)):
        print("  NOTE: the master was not in canonical form, so this write reformats")
        print("        every line (JSON::XS pretty+canonical). One-time; after this a")
        print("        diff shows only real changes, and Perl tools stop reordering it.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".", help="module working directory")
    ap.add_argument("--repo", required=True, help="Viral_Annotation checkout")
    ap.add_argument("--json", help="module JSON (default: <workdir>/*_Viral_PSSM.json)")
    ap.add_argument("--check", action="store_true", help="validate only, write nothing")
    ap.add_argument("--only", help="comma-separated module names to install")
    args = ap.parse_args()

    workdir = os.path.abspath(args.workdir)
    repo = os.path.abspath(os.path.expanduser(args.repo))
    jf = args.json or find_module_json(workdir)
    with open(jf) as fh:
        mod_json = json.load(fh, object_pairs_hook=OrderedDict)
    print("module JSON: %s (%d modules)\n" % (jf, len(mod_json)))

    print("validating:")
    modules, problems = validate(workdir, mod_json)
    if args.only:
        keep = set(args.only.split(","))
        modules = {k: v for k, v in modules.items() if k in keep}

    if args.check:
        print("\n--check: nothing written. %d installable, %d problem(s)."
              % (len(modules), len(problems)))
        return 0
    if not modules:
        sys.exit("\nnothing installable")

    print("\ninstalling into %s:" % repo)
    install(workdir, repo, modules, mod_json)
    print("\nrun with:  export LOWVAN_DATA_DIR=%s" % repo)
    return 0


if __name__ == "__main__":
    sys.exit(main())
