"""Solve the single 33K intron per genome, exactly where a deposit exists.

build_splice_refs.py places the exon boundary at the best canonical GT..AG
inside a nudge window around the tblastn HSP ends.  Barthadenovirus 33K
defeats that: its intron INTERRUPTS A CODON, so the true boundary is not at
a codon boundary, and the window then finds a different, downstream AG that
also looks canonical.  On 2169706.5 that gave a 150 nt intron where the real
one is 99, and a 120 aa product where the deposit is 137.

Here the deposit is the target, not a guide: a junction is accepted only when
the spliced translation reproduces a known 33K protein.  Genomes with their
own deposit are solved to byte-identity; the rest are solved against the
validated products of their nearest relatives.
"""
import subprocess, sys, os, json

B1='T'*16+'C'*16+'A'*16+'G'*16; B2=('T'*4+'C'*4+'A'*4+'G'*4)*4; B3='TCAG'*16
AA="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
COD={B1[i]+B2[i]+B3[i]:AA[i] for i in range(64)}
tr=lambda s:''.join(COD.get(s.upper()[i:i+3],'X') for i in range(0,len(s)-2,3))
rc=lambda s: s.translate(str.maketrans("acgtACGT","tgcaTGCA"))[::-1]

def rdfa(p):
    d=[];n=None;s=[]
    for l in open(p):
        if l[0]=='>':
            if n is not None: d.append((n,''.join(s)))
            n=l[1:].strip();s=[]
        else: s.append(l.strip())
    if n is not None: d.append((n,''.join(s)))
    return d

def idf(a,b):
    return sum(x==y for x,y in zip(a,b))/max(len(a),len(b),1)

def hsps(qfile, contig):
    r=subprocess.run(["tblastn","-query",qfile,"-subject",contig,
        "-outfmt","6 qstart qend sstart send bitscore","-seg","no","-evalue","1e-3"],
        capture_output=True,text=True)
    return [x.split('\t') for x in r.stdout.strip().split('\n') if x]

def solve(contig, targets, tmpq, intron_lo=50, intron_hi=200, want=0.80):
    """Return (best_product, start, donor, acceptor, strand, identity, target)."""
    G=''.join(l.strip() for l in open(contig) if l[0]!='>')
    L=len(G)
    best=None
    for tname,P in targets:
        open(tmpq,"w").write(">q\n%s\n"%P)
        rows=hsps(tmpq, contig)
        if not rows: continue
        for strand in (True,False):
            sel=[r for r in rows if (int(r[2])<int(r[3]))==strand]
            if not sel: continue
            S = G if strand else rc(G)
            conv=(lambda p:p) if strand else (lambda p: L-p+1)
            lo=min(min(conv(int(r[2])),conv(int(r[3]))) for r in sel)
            hi=max(max(conv(int(r[2])),conv(int(r[3]))) for r in sel)
            #  candidate ATG starts just upstream of where the match begins
            starts=[s for s in range(max(1,lo-30), lo+31) if S[s-1:s+2].upper()=="ATG"]
            region=S[lo-1:min(len(S),hi+400)]
            gts=[lo-1+i for i in range(len(region)-1) if region[i:i+2].upper()=="GT"]
            for s in starts:
                for d in gts:
                    if d < s+2: continue
                    for a in range(d+intron_lo, min(d+intron_hi, len(S))):
                        if S[a-2:a].upper()!="AG": continue
                        #  Translate a WINDOW, not the rest of the contig. The
                        #  product cannot be much longer than the target, and
                        #  translating ~30 kb on every one of the millions of
                        #  (start, donor, acceptor) candidates was the whole
                        #  cost: Aviadenovirus ran at 50 s/genome, which put
                        #  Mastadenovirus's 343 deposits near five hours.
                        win=3*(len(P)+60)
                        pre=S[s-1:d]+S[a:a+win]
                        p=tr(pre)
                        stop=p.find('*')
                        if stop<0: continue
                        p=p[:stop]
                        if 'X' in p or not p.startswith('M'): continue
                        #  The product must actually USE the junction. If the
                        #  ORF from the start codon hits its stop before the
                        #  donor, exon B is never reached and (d, a) are
                        #  arbitrary -- the "solution" is an unspliced ORF that
                        #  happens to match the target. That inflated the exact
                        #  counts badly before it was caught: 52 of 307 on
                        #  Mastadenovirus, 42 of 124 on Aviadenovirus.
                        if (d-(s-1)) >= 3*(len(p)+1): continue
                        sc=idf(p,P)
                        if best is None or sc>best[5]:
                            best=(p,s,d,a,strand,sc,tname)
                        #  Nothing beats reproducing the deposit exactly.
                        if sc==1.0: return best
    if best and best[5]>=want: return best
    return None
