#!/bin/bash
#  Literature length references, one file per module. Every range needs a PMID.
set -eu
cd "$HOME/lowvan_scratch/tobani"
cat > Torovirinae/length_refs.json <<'EOF'
{
 "ORF1A": {"min": 4300, "max": 4700, "pmid": "16571831",
   "note": "EToV Berne pp1a 4,569 aa (Smits 2006; GenBank DQ310701 ORF1a 822..14531); BToV/PToV pp1a 4,417-4,418 aa"},
 "S":     {"min": 1500, "max": 1650, "pmid": "2219698",
   "note": "Berne virus peplomer apoprotein 1,581 aa (Snijder 1990); RefSeq NC_007447 S 20975..25729 = 1,584 aa"},
 "M":     {"min": 220, "max": 245, "pmid": "2024492",
   "note": "torovirus triple-spanning envelope protein, now called M, 233 aa in Berne virus (den Boon 1991)"},
 "HE":    {"min": 400, "max": 440, "pmid": "9188596",
   "note": "BToV hemagglutinin-esterase, 65 kDa glycoprotein with sialate-O-acetylesterase activity (Cornelissen 1997); RefSeq NC_007447 HE 26477..27736 = 419 aa"},
 "N":     {"min": 155, "max": 175, "pmid": "2607340",
   "note": "Berne virus nucleocapsid protein 160 aa (Snijder 1989)"}
}
EOF
cat > Piscanivirinae/length_refs.json <<'EOF'
{
 "ORF1A": {"min": 4400, "max": 5000, "pmid": "16987966",
   "note": "white bream virus pp1a 4,555 aa (Schuetze 2006; NC_008516 ORF1a 906..14573); fathead minnow nidovirus 4,860 aa"},
 "S":     {"min": 1150, "max": 1250, "pmid": "16987966",
   "note": "WBV S 1,220 aa (NC_008516 21525..25187); FHMNV S 1,190 aa (Batts 2012, PMID 22422065); CSBV 1,205 aa"},
 "M":     {"min": 215, "max": 250, "pmid": "16987966",
   "note": "WBV M 227 aa; FHMNV M 225 aa; CSBV M 245 aa"},
 "N":     {"min": 155, "max": 200, "pmid": "16987966",
   "note": "WBV N 161 aa; FHMNV N 165 aa (PMID 22422065); CSBV N 177-178 aa"}
}
EOF
cat > Remotovirinae/length_refs.json <<'EOF'
{
 "ORF1A": {"min": 3000, "max": 3200, "pmid": "25918239",
   "note": "bovine nidovirus TCH5; KM589359 CDS join(487..9825,9825..15818), pp1a 3,113 aa (Tokarz 2015)"},
 "S":     {"min": 540, "max": 580, "pmid": "25918239",
   "note": "ORF2 spike, 1,689 nt / 562 aa, the shortest nidovirus spike reported"},
 "M":     {"min": 225, "max": 240, "pmid": "25918239", "note": "ORF3 membrane protein, 231 aa, multiple TMs"},
 "N":     {"min": 170, "max": 185, "pmid": "25918239", "note": "ORF4 nucleocapsid, 178 aa, predicted pI 12.0"},
 "G2":    {"min": 420, "max": 470, "pmid": "25918239",
   "note": "ORF5 'G2', 455 aa putative accessory glycoprotein, TM 429-448, several glycosylation sites"},
 "ORF6":  {"min": 80, "max": 95, "pmid": "25918239",
   "note": "the unnamed 267 nt / 88 aa ORF near the 3' end; no similarity to any known protein"}
}
EOF
cat > Serpentovirinae/length_refs.json <<'EOF'
{
 "ORF1A": {"min": 3800, "max": 6700, "pmid": "25205093",
   "note": "ball python nidovirus ORF1a 17,412 nt / 5,803 aa (Stenglein 2014); the subfamily spans 3,840 aa (shingleback nidovirus) to 6,570 aa (veiled chameleon serpentovirus B)"},
 "S":     {"min": 800, "max": 1200, "pmid": "25205093",
   "note": "BPNV ORF2 S, 959 aa; subfamily range 821-1,141 aa"},
 "M":     {"min": 200, "max": 270, "pmid": "25205093",
   "note": "BPNV ORF4 M, 215 aa, three N-terminal TMs, 21% identity over 116 aa to white bream virus M"},
 "N":     {"min": 125, "max": 250, "pmid": "25205093",
   "note": "BPNV ORF5a N, 152 aa, Arteri_nucleo Pfam hit, pI 10.6; divergent genera reach 129-246 aa and are bound to this protein by homology, never by length"}
}
EOF
echo "wrote 4 length_refs.json"
