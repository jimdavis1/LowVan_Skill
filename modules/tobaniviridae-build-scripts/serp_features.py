#!/usr/bin/env python3
"""Write Serpentovirinae/features.json."""
import json, os
os.chdir(os.path.expanduser('~/lowvan_scratch/tobani'))
f = {
 "ORF1A": {"anno": "ORF1a polyprotein", "gene_symbol": "ORF1a", "segment": "Single RNA Segment",
  "downstream_ext": 1,
  "note": "pp1a, 3,840-6,570 aa; ball python nidovirus ORF1a is 5,803 aa and every domain in it -- protein kinase, transmembrane regions, main protease, and in ORF1b the RdRp, Zn-helicase, exoribonuclease, NendoU and O-methyltransferase -- is an HMM prediction. The describing paper says functional experiments are still needed to validate the cleavage pattern, so this module declares no mature peptide.",
  "PMID_claude_generated": {"function": ["25205093"]}},
 "ORF1AB": {"anno": "ORF1ab polyprotein", "gene_symbol": "ORF1ab", "segment": "Single RNA Segment",
  "special": "transcript_edit", "feature_type": "CDS",
  "note": "Expressed by -1 ribosomal frameshifting. ORF1a and ORF1b overlap by 52 nt and the overlap carries the A-rich slippery heptamer AAAAAAC, not the UUUAAAC of the rest of this family (PMID 25205093; NC_024709 misc_feature 18405..18411, CDS join(1018..18411,18411..25343)); the signal was independently described as the astrovirus/retrovirus type (PMID 25063552). It is putative: no serpentovirus frameshift has been measured. No PSSM.",
  "PMID_claude_generated": {"function": ["25205093", "25063552"]}},
 "S": {"anno": "Spike protein", "gene_symbol": "S", "segment": "Single RNA Segment", "downstream_ext": 1,
  "note": "ORF2, 821-1,141 aa. The S2 fusion domain is detectable (18% identity over 333 aa to a coronavirus spike) while S1 has no detectable homology (PMID 25205093). No cleavage site is published, so no subunits are declared.",
  "PMID_claude_generated": {"function": ["25205093"]}},
 "M": {"anno": "Membrane glycoprotein", "gene_symbol": "M", "segment": "Single RNA Segment", "downstream_ext": 1,
  "note": "ORF4, 208-264 aa, three N-terminal transmembrane segments, 21% identity over 116 aa to the white bream virus M (PMID 25205093). The subfamily holds two orthogroups of it, 28% identical to each other over 163 aa; both are bound here by homology to the protein that paper names, not by position.",
  "PMID_claude_generated": {"function": ["25205093"]}},
 "N": {"anno": "Nucleocapsid protein", "gene_symbol": "N", "segment": "Single RNA Segment", "downstream_ext": 1,
  "note": "ORF5a, 129-246 aa; the ball python nidovirus protein is 152 aa with a Pfam Arteri_nucleo hit and a predicted pI of 10.6 (PMID 25205093). Antisera against it detect infection in green tree pythons (PMID 28794044).",
  "PMID_claude_generated": {"function": ["25205093", "28794044"]}},
}
UNC = {
 "UNC1": "the 3'-terminal accessory of Pregotovirus (499-510 aa) together with the ORF reported to resemble the paramyxovirus haemagglutinin-neuraminidase (479-506 aa, PMID 25063552). They are ONE feature by measurement: as separate features their profiles called each other's loci at identical coordinates in Xinzhou nematode virus 6 and Bellinger River virus, although blastp finds no HSP between the two collections even at e=10. The HN resemblance is a sequence similarity, not a demonstrated function",
 "UNC2": "ORF5b-position accessory of Pregotovirus, 273-417 aa",
 "UNC3": "ORF3-position accessory of Pregotovirus, 232-256 aa, deposited both as 'minor M protein' and as 'TM-glycoprotein'",
 "UNC4": "Infratovirus accessory between S and M, 460-479 aa",
 "UNC5": "accessory immediately after S in Sectovirus and the snake viruses BV-BRC files under Torovirus, 269-321 aa",
 "UNC6": "second accessory after S in Sectovirus and those snake viruses, 338-381 aa",
 "UNC7": "3'-terminal accessory of Sectovirus and those snake viruses, 328 aa",
 "UNC8": "accessory of Lyctovirus and those snake viruses, 215-227 aa",
 "UNC9": "3'-terminal accessory of Lyctovirus and those snake viruses, 675-681 aa",
}
for k, n in UNC.items():
    f[k] = {"anno": "Uncharacterized lineage-specific protein", "gene_symbol": k,
            "segment": "Single RNA Segment", "downstream_ext": 1,
            "note": "%s. No function is published for any serpentovirus accessory ORF (PMID 25205093). One profile set per orthogroup; genome position was used only to tell one orthogroup from another." % n,
            "PMID_claude_generated": {"function": ["25205093", "25063552"] if k == "UNC1" else ["25205093"]}}
json.dump(f, open('Serpentovirinae/features.json', 'w'), indent=1, sort_keys=True)
print(len(f), 'features')
