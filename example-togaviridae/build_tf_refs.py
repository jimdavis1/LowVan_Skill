#!/usr/bin/env python3
"""Construct the transcript-edited TF reference set for Togaviridae.

TF arises from a -1 ribosomal frameshift at the conserved UUUUUUA heptamer inside
6K (Firth 2008, PMID 18822126).  get_transcript_edited_features.pl reconstructs it
by blastn-ing a curated reference that carries one extra templated base, then
filling the resulting subject gap from the query.  So the reference is the 6K CDS
read in frame up to the slip, with the slipped base duplicated, continuing in the
-1 frame to the stop, WITH the stop codon retained.

The star is the convention: a protein that really terminates at a stop carries a
trailing "*", one that ends at a cleavage site does not. gjoseqlib::translate_seq
already implements exactly that -- it emits "*" only for an unambiguous stop
codon, and "x" for an ambiguous or partial one -- so the reference decides. TF
ends at a genuine TAA in the -1 frame, so its reference keeps it. The precedent
is ssGP (74/74 references end in a stop), not NSP12 (0/3006), which is a protease
cleavage product with no stop of its own.

The heptamer sits at phase 2 of the 6K reading frame: X_XXY_YYZ puts the lone X
at the third position of a codon.  An earlier survey tested for phase 0 and so
rejected every genome; that was the bug, not a finding.
"""
import collections, glob, os, sys

T={}; _b="TCAG"; _a="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"; i=0
for x in _b:
    for y in _b:
        for z in _b: T[x+y+z]=_a[i]; i+=1
tr=lambda s:"".join(T.get(s[j:j+3],"X") for j in range(0,len(s)-2,3))
rc=lambda s:s[::-1].translate(str.maketrans("ACGTN","TGCAN"))

MOTIF="TTTTTTA"
FLANK=900          # nt downstream of the 6K CDS to search for the -1 stop
WINDOW=(100,170)   # offset range of the slip site from the 6K start

def load_fna(p):
    s=[]
    for L in open(p,errors="replace"):
        if not L.startswith(">"): s.append(L.strip())
    return "".join(s).upper()

stat=collections.Counter(); recs=[]; fails=[]
for tbl in sorted(glob.glob("coverage_eval_h/ann/*.feature.tbl")):
    if not os.path.getsize(tbl): continue
    gid=os.path.basename(tbl).replace(".feature.tbl","")
    rows=[]
    for line in open(tbl,errors="replace"):
        c=line.rstrip("\n").split("\t")
        if len(c)>=16 and c[6]=="6K":
            try: rows.append((int(c[7]),int(c[8]),c[9],c[15]))
            except ValueError: pass
    if len(rows)!=1:
        stat["6K not called exactly once"]+=1; continue
    start,end,strand,species=rows[0]
    fna=f"coverage_eval/ex/{gid}.fna"
    if not os.path.exists(fna): stat["no sequence"]+=1; continue
    g=load_fna(fna)

    # extended CDS: 6K plus downstream, in the coding orientation
    if strand=="+": ext=g[start-1:end+FLANK]
    else:           ext=rc(g[max(0,start-1-FLANK):end])
    if len(ext)<200: stat["sequence too short"]+=1; continue

    cds_len=end-start+1
    # The slip site sits at codon 39-46 of 6K in every genome that has one
    # (190/279 at exactly codon 46).  Anchoring the search to that window rather
    # than to the called 6K interval keeps two spurious upstream heptamers out,
    # and still finds the site when the 6K call is truncated by an N-run.
    lo,hi=WINDOW
    hits=[j for j in range(lo,min(hi,len(ext)-len(MOTIF))+1)
          if ext[j:j+len(MOTIF)]==MOTIF and j%3==2]
    if not hits:
        any_motif=MOTIF in ext[:cds_len]
        stat["motif at phase 2 absent"+(" (motif present, other phase)" if any_motif else "")]+=1
        fails.append((gid,species,"no phase-2 heptamer in 6K")); continue
    if len(hits)>1: stat["multiple phase-2 heptamers (took first)"]+=1
    h=hits[0]

    RES=3*(h//3+2)          # first nt past the last in-frame codon
    DUP=RES-1               # the nucleotide decoded twice
    tail=tr(ext[DUP:])
    k=tail.find("*")
    if k<0: stat["no -1 stop within flank"]+=1; fails.append((gid,species,"no -1 stop")); continue
    END=DUP+3*(k+1)         # exclusive, includes the stop codon

    ref=ext[:RES]+ext[DUP:END]          # stop codon RETAINED -- see below
    prot=tr(ref)
    if prot.count("*")!=1 or not prot.endswith("*"):
        stat["internal stop in construct"]+=1; fails.append((gid,species,"internal stop")); continue
    if "X" in prot:  stat["ambiguous base in construct"]+=1; fails.append((gid,species,"ambiguous base")); continue
    if len(ref)%3:   stat["construct not in frame"]+=1; continue
    # the pre-slip portion must agree with the 6K reading frame
    if prot[:RES//3]!=tr(ext[:RES]): stat["prefix disagrees with 6K frame"]+=1; continue

    stat["OK"]+=1
    recs.append((gid,species,ref,prot,h,k))

print("  " + "\n  ".join(f"{v:>5}  {k}" for k,v in stat.most_common()))
print(f"\n  constructed {len(recs)} TF references")
if recs:
    import statistics as st
    L=[len(r[3].rstrip("*")) for r in recs]
    print(f"  protein length  min {min(L)}  median {int(st.median(L))}  max {max(L)}")
    U=len({r[2] for r in recs})
    print(f"  unique nucleotide sequences: {U}")

os.makedirs("TF_refs",exist_ok=True)
with open("TF_refs/TF.fasta","w") as f:
    for gid,species,ref,prot,h,k in recs:
        f.write(f">TF|{gid} Togaviridae TF protein [{species} | {gid}]\n")
        for j in range(0,len(ref),60): f.write(ref[j:j+60].lower()+"\n")
with open("TF_refs/TF.faa","w") as f:
    for gid,species,ref,prot,h,k in recs:
        f.write(f">TF|{gid} {species}\n")
        for j in range(0,len(prot),60): f.write(prot[j:j+60]+"\n")
with open("TF_refs/failures.tsv","w") as f:
    f.write("genome\tspecies\treason\n")
    for r in fails: f.write("\t".join(r)+"\n")
print(f"\n  wrote TF_refs/TF.fasta, TF_refs/TF.faa, TF_refs/failures.tsv")
