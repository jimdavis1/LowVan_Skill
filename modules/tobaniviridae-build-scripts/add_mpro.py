#!/usr/bin/env python3
"""Declare MPRO in both modules that have a published cleavage site, and give
it a literature length range for the audit."""
import json, os
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
TORO = ("The torovirus main proteinase, cut from pp1a at its own two published "
        "cleavage sites rather than modelled from deposits -- nothing in BV-BRC "
        "is annotated as a torovirus mature peptide. Both sites were mapped by "
        "Edman degradation of the cleavage products: equine torovirus "
        "Gln3251|Ser3252 in SNFSFQ|SVFSKA and Gln3543|Ser3544 in QPFKKQ|SVSNNV, "
        "mature protein 292 aa (PMID 16571831, GenBank DQ310701); porcine "
        "torovirus Gln3095|Ser3096 and Gln3382|Ser3383, mature protein 287 aa "
        "(PMID 32727876, GenBank MH603532). Measured over all 32 pp1a sequences "
        "in this module's collection, P1 is Gln and P1' is Ser at 100% at BOTH "
        "ends, and the derived products are 287-292 aa. It is a chymotrypsin-"
        "like SERINE protease (nucleophile Ser, mature numbering Ser165) that "
        "nevertheless cleaves after Gln, consensus FxxQ|(S,A).")
PISC = ("The bafinivirus main protease, cut from pp1a by anchoring on the "
        "protein whose ends were sequenced -- white bream virus pp1a 3443-3709, "
        "the N1-C1 marker of PMID 21068254 -- and snapping each boundary to the "
        "nearest Gln|(Ala,Thr,Ser), that protease's published consensus. All "
        "four of its autoprocessing sites were mapped by Edman degradation: N1 "
        "Gln3442|Ala3443 (ITRQ|ALQR, predominant), N2 Gln3449|Ala3450, C1 "
        "Gln3709|Thr3710 (VTQQ|TSVT) and C2 Gln3716|Thr3717; catalytic triad "
        "Ser3589-His3492-Asp3518. Measured here over the collection: 5 of 6 "
        "pp1a sequences yield a product, P1 is Gln at 100% at both ends, "
        "products 267-275 aa. The sixth is refused rather than guessed -- no "
        "Gln|(Ala,Thr,Ser) within 12 residues of the projected boundary.")
for M, note, lo, hi, pmids in (
        ('Torovirinae', TORO, 280, 300, ["16571831", "32727876"]),
        ('Piscanivirinae', PISC, 255, 285, ["21068254"])):
    f = json.load(open('%s/features.json' % M))
    f['MPRO'] = {
        "anno": "Main protease involved in polyprotein cleavage",
        "gene_symbol": "Mpro",
        "segment": "Single RNA Segment",
        "feature_type": "mat_peptide",
        #  Both termini are protease cuts, so neither may scan outward: there
        #  is no codon at either end for the annotator to find.
        "upstream_ext": 0,
        "downstream_ext": 0,
        "params": {"m": 2},
        "note": note,
        "PMID_claude_generated": {"function": pmids, "coordinates": pmids},
    }
    json.dump(f, open('%s/features.json' % M, 'w'), indent=1, sort_keys=True)
    r = json.load(open('%s/length_refs.json' % M))
    r['MPRO'] = {"min": lo, "max": hi, "pmid": pmids[0],
                 "note": ("mature main protease: 292 aa in equine torovirus and 287 aa in "
                          "porcine torovirus" if M == 'Torovirinae' else
                          "mature main protease: white bream virus pp1a 3443-3709 = 267 aa")}
    json.dump(r, open('%s/length_refs.json' % M, 'w'), indent=1, sort_keys=True)
    print(M, 'MPRO declared')
