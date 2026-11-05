#!/usr/bin/env python3
"""Audit nested alternative-start collections for mislabelled members.

VP1/VP2/VP3 and Rep78/68/52/40 are ~100% identical where they overlap, so a
deposit filed under the wrong form cannot be separated again by homology, by
self-recall or by call rate. The only signal is LENGTH, read per species and
against what the literature says each form measures.

Why this exists. Parvoviridae, 30 Oct 2026: 40 unique deposits named "VP2" in
Erythroparvovirus were 781 aa, which is B19 VP1 (VP2 is 554 aa -- VP1 without
the 227-residue VP1u, PMID 11799199). fasta-cluster-pssm then dropped the two
genuine VP2 clusters because a consensus that is a SUBSTRING of another
cluster's reads as truncated, and VP2 shipped as a VP1 profile. The module
called VP2 at VP1's exact coordinates in 281 of 295 genomes while every other
check passed. TetraCopiparvovirus carried the same defect in both directions.

The test that settles which form a sequence is, is structural rather than
nominal: VP2 is an exact C-terminal substring of VP1 FROM THE SAME GENOME.
Where a genome deposits both, the longer is VP1 and the shorter is VP2
whatever the labels say.

    python3 capsid_length_audit.py --workdir . --module Erythroparvovirus
    python3 capsid_length_audit.py --workdir . --module TetraCopiparvovirus \\
            --forms VP,VP2,VP3

Exit status is 1 if anything is reported, so it can gate a build.
"""
import argparse, collections, glob, os, re, sys


def read_fasta(path):
    out, h, s = [], None, []
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if h is not None:
                    out.append((h, "".join(s)))
                h, s = line[1:].split()[0], []
            else:
                s.append(line.strip())
    if h is not None:
        out.append((h, "".join(s)))
    return out


def genome_of(fid):
    """'172297.3|VP_1' -> '172297.3'."""
    return fid.split("|")[0]


def species_of(name, gid=None):
    """Species key: the taxon id when there is one, else the genome name.

    The point is to group deposits that share a VP1 length, not to do
    taxonomy. A per-GENUS comparison is not good enough: Tetraparvovirus
    alone holds VP1 at 914, 925, 931, 933 and 1032 aa, so its genus median
    declares every real VP1 below 1032 a VP2.

    PREFER THE TAXON ID. The genome id's prefix IS the NCBI taxon id --
    `10524.1311` -> `10524` -- so it is exact and free, where the name is a
    string that submitters vary. The first-three-words heuristic below holds
    for "Human adenovirus 41 Hannover_2022" but gives a DIFFERENT key for
    every isolate of "unidentified adenovirus MCAdV/XJ11C", and a per-isolate
    key means a median computed over one member, which clears every finding
    it should raise. The name remains the fallback for callers that have a
    name and no id.
    """
    if gid:
        tax = str(gid).split("|")[0].split(".")[0]
        if tax:
            return tax
    n = re.sub(r"\s+(isolate|strain|clone)\b.*", "", name or "")
    bits = n.split()
    return " ".join(bits[:3]) if bits else "(unknown)"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=".")
    ap.add_argument("--module", required=True)
    ap.add_argument("--forms", default="",
                    help="comma-separated nested forms to audit, longest "
                         "first (e.g. VP,VP2,VP3). Default: every feature "
                         "that nests with another inside one genome.")
    ap.add_argument("--names",
                    help="<id>\\t<name>\\t<family>\\t<genus> table. Default: "
                         "<workdir>/*.id_name_gs")
    ap.add_argument("--tolerance", type=float, default=0.97,
                    help="a member within this fraction of ANOTHER form's "
                         "length class for its own species is suspect")
    a = ap.parse_args()

    W, M = a.workdir, a.module
    cdir = os.path.join(W, "collections", M)
    if not os.path.isdir(cdir):
        sys.exit("no collections at %s" % cdir)

    names = {}
    npath = a.names
    if not npath:
        hits = sorted(glob.glob(os.path.join(W, "*.id_name_gs")))
        npath = hits[0] if hits else None
    if npath and os.path.exists(npath):
        for line in open(npath, errors="replace"):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                names[f[0]] = f[1]

    #  feature -> [(genome, seq)]
    feats = {}
    for p in sorted(glob.glob(os.path.join(cdir, "*.fasta"))):
        k = os.path.basename(p)[:-6]
        feats[k] = [(genome_of(h), s) for h, s in read_fasta(p)]
    if not feats:
        sys.exit("no collections in %s" % cdir)

    forms = [x.strip() for x in a.forms.split(",") if x.strip()]
    problems = 0

    #  ---- 1. the structural test: nesting inside ONE genome --------------
    #
    #  No names involved. If one genome's deposit is an exact C-terminal
    #  substring of another of its own deposits, the longer is the parent
    #  form. Where that contradicts the feature each one sits in, the labels
    #  are inverted and the collections are cross-contaminated.
    bygen = collections.defaultdict(list)
    for k, rows in feats.items():
        if forms and k not in forms:
            continue
        for g, s in rows:
            bygen[g].append((k, s))

    inverted = collections.Counter()
    nested_pairs = collections.Counter()
    for g, rows in bygen.items():
        for i, (k1, s1) in enumerate(rows):
            for k2, s2 in rows[i + 1:]:
                if k1 == k2 or not s1 or not s2:
                    continue
                lo, hi = (s1, s2) if len(s1) < len(s2) else (s2, s1)
                klo = k1 if lo is s1 else k2
                khi = k2 if lo is s1 else k1
                if len(lo) == len(hi) or not hi.endswith(lo):
                    continue
                nested_pairs[(khi, klo)] += 1
                if forms and khi in forms and klo in forms \
                        and forms.index(khi) > forms.index(klo):
                    inverted[(khi, klo, len(hi), len(lo))] += 1

    if nested_pairs:
        print("  nested pairs within one genome (longer form -> shorter):")
        for (khi, klo), n in nested_pairs.most_common():
            print("    %-6s contains %-6s  x%d" % (khi, klo, n))
    if inverted:
        problems += sum(inverted.values())
        print("\n  LABELS INVERTED -- the longer sequence sits in the form")
        print("  declared shorter. The deposited names disagree with the")
        print("  structure, and the structure is the evidence:")
        for (khi, klo, lhi, llo), n in inverted.most_common():
            print("    %s=%d aa is the parent of %s=%d aa  x%d"
                  % (khi, lhi, klo, llo, n))

    #  ---- 2. length classes per species ----------------------------------
    #
    #  Per SPECIES. A length class is claimed by whichever feature holds most
    #  of its members; the rest are the suspects. This is what catches the
    #  case where a whole species deposits only the wrong name, so there is
    #  no same-genome pair to compare against.
    #
    #  Only forms that ACTUALLY nest. Without --forms this test used to run
    #  over every feature in the module, and two unrelated proteins that
    #  happen to measure the same reported each other: Bocaparvovirus NS1 and
    #  VP1 are both ~685 aa in porcine bocavirus 3, which is a coincidence
    #  and not a mislabelled alternative start.
    #  Forms are compared only against the ones they actually nest WITH, as
    #  connected components of the nesting graph. A flat set is not enough:
    #  Bocaparvovirus nests VP/VP2 and, in two deposits, NS1/NP1, and with
    #  one combined set NS1 and VP1 -- both ~685 aa in porcine bocavirus 3
    #  and unrelated -- reported each other.
    if forms:
        comps = [set(forms)]
    else:
        comps = []
        for khi, klo in nested_pairs:
            hit = [c for c in comps if khi in c or klo in c]
            new = set([khi, klo])
            for c in hit:
                new |= c
                comps.remove(c)
            comps.append(new)

    owner = {}
    cls = collections.defaultdict(collections.Counter)
    for k, rows in feats.items():
        grp = [i for i, c in enumerate(comps) if k in c]
        if not grp:
            continue
        for g, s in rows:
            cls[(grp[0], species_of(names.get(g, ""), g), len(s))][k] += 1
    for key, q in cls.items():
        owner[key] = q.most_common(1)[0][0]

    minority = collections.Counter()
    for key, q in sorted(cls.items()):
        win = owner[key]
        _, sp, L = key
        for k, n in q.items():
            if k != win:
                minority[(sp, L, k, win)] += n

    if minority:
        problems += sum(minority.values())
        print("\n  MINORITY LABEL ON A LENGTH CLASS -- for this species, this")
        print("  length is mostly filed as the other form:")
        for (sp, L, k, win), n in minority.most_common(20):
            print("    %-34s %4d aa  %-6s x%-4d (majority says %s)"
                  % (sp[:34], L, k, n, win))

    #  ---- 3. a form carrying its parent's length class --------------------
    if forms:
        med = {}
        for k in forms:
            per = collections.defaultdict(list)
            for g, s in feats.get(k, []):
                per[species_of(names.get(g, ""), g)].append(len(s))
            for sp, L in per.items():
                L.sort()
                med[(k, sp)] = L[len(L) // 2]
        carry = collections.Counter()
        for idx, k in enumerate(forms):
            if idx == 0:
                continue
            for g, s in feats.get(k, []):
                sp = species_of(names.get(g, ""), g)
                parent = med.get((forms[0], sp))
                if parent and len(s) >= a.tolerance * parent:
                    carry[(sp, k, len(s), parent)] += 1
        if carry:
            problems += sum(carry.values())
            print("\n  NESTED FORM AT THE PARENT'S LENGTH -- these cannot be")
            print("  separated by any bit_cutoff and will poison the profile:")
            for (sp, k, L, par), n in carry.most_common(20):
                print("    %-34s %-6s %4d aa vs %s median %d aa  x%d"
                      % (sp[:34], k, L, forms[0], par, n))

    #  ---- 4. which END does it share with the parent? --------------------
    #
    #  A nested form is an alternative START: it carries the parent's
    #  C-terminus and is missing the parent's N-terminus. A deposit that
    #  begins at the PARENT's residue 1 and stops short is the opposite of
    #  that -- a truncated parent under the wrong name -- and no length rule
    #  can tell the two apart, because both are simply "shorter".
    #
    #  Every "VP2" deposit of porcine parvovirus 2 (606, 654, 713 aa) and of
    #  porcine parvovirus 4 (394, 427, 485 aa) starts at VP1 residue 1. They
    #  built two VP2 profiles out of VP1's N-terminus, which then called VP2
    #  at VP1's own coordinates in 132 genomes.
    for comp in comps:
        order = [f for f in (forms or sorted(comp)) if f in comp]
        if len(order) < 2:
            continue
        par = order[0]
        #  Against EVERY parent member of the species, not just its longest:
        #  porcine parvovirus 4's longest VP deposit is a 914 aa outlier
        #  while the species sits at 728, and its three N-anchored fragments
        #  start at the 728 aa form.
        ref = collections.defaultdict(list)
        for g, s in feats.get(par, []):
            ref[species_of(names.get(g, ""), g)].append(s)
        nanch = collections.Counter()
        for k in order[1:]:
            for g, s in feats.get(k, []):
                sp = species_of(names.get(g, ""), g)
                #  By IDENTITY over the first 25 residues, not an exact
                #  prefix. Two B19 deposits of 602 and 603 aa carry VP1's own
                #  MSKESGKWWES... start but differ from the 781 aa VP1 at
                #  residue 18, so an exact match missed them and they shipped
                #  a VP2 profile that called VP2 at VP1's coordinates in 269
                #  of 300 genomes. A real VP2 head is VP1's middle: 2-3 of 25.
                hit = [r for r in ref.get(sp, ())
                       if len(s) < 0.95 * len(r)
                       and sum(a == b for a, b in zip(s[:25], r[:25])) >= 20]
                if hit:
                    nanch[(sp, k, len(s), len(max(hit, key=len)))] += 1
        if nanch:
            problems += sum(nanch.values())
            print("\n  N-ANCHORED TO THE PARENT -- these share %s's START, not"
                  % par)
            print("  its end, so they are truncated %s, not an alternative"
                  % par)
            print("  start. They belong in %s, where truncation QC drops"
                  % par)
            print("  them, never in a nested sub-form:")
            for (sp, k, L, R), n in nanch.most_common(20):
                print("    %-34s %-6s %4d aa begins at %s residue 1 (%d aa)"
                      " x%d" % (sp[:34], k, L, par, R, n))

    if not problems:
        print("  %s: no mislabelled nested member found" % M)
        return 0
    print("\n  %d suspect member(s). Re-bin on the literature's size, then"
          % problems)
    print("  REBUILD the alignment and PSSM -- a contaminated collection does")
    print("  not become clean by renaming rows in a report.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
