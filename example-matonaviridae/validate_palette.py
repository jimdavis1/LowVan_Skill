#!/usr/bin/env python3
"""Python port of the dataviz skill's validate_palette.js (no node on this box).

Same thresholds and the same Machado-Oliveira-Fernandes (2009) severity-1.0
CVD transforms, so the numbers are comparable to the JS tool's.
"""
import math,sys
BAND={'light':(0.43,0.77),'dark':(0.48,0.67)}
CHROMA_FLOOR=0.10; CVD_TARGET=8.0; CVD_FLOOR=6.0; NORMAL_FLOOR=15.0; CONTRAST_MIN=3.0
MACHADO={'protan':[[0.152286,1.052583,-0.204868],[0.114503,0.786281,0.099216],[-0.003882,-0.048116,1.051998]],
         'deutan':[[0.367322,0.860646,-0.227968],[0.280085,0.672501,0.047413],[-0.011820,0.042940,0.968881]],
         'tritan':[[1.255528,-0.076749,-0.178779],[-0.078411,0.930809,0.147602],[0.004733,0.691367,0.303900]]}
def hex2rgb(h):
    h=h.lstrip('#'); return [int(h[i:i+2],16)/255 for i in (0,2,4)]
def s2l(c): return c/12.92 if c<=0.04045 else ((c+0.055)/1.055)**2.4
def l2s(c):
    c=max(0.0,min(1.0,c)); return c*12.92 if c<=0.0031308 else 1.055*c**(1/2.4)-0.055
def lin(h): return [s2l(c) for c in hex2rgb(h)]
def lrgb2oklab(r,g,b):
    l=0.4122214708*r+0.5363325363*g+0.0514459929*b
    m=0.2119034982*r+0.6806995451*g+0.1073969566*b
    s=0.0883024619*r+0.2817188376*g+0.6299787005*b
    l_,m_,s_=[x**(1/3) if x>0 else -((-x)**(1/3)) for x in (l,m,s)]
    return (0.2104542553*l_+0.7936177850*m_-0.0040720468*s_,
            1.9779984951*l_-2.4285922050*m_+0.4505937099*s_,
            0.0259040371*l_+0.7827717662*m_-0.8086757660*s_)
def oklch(h):
    L,a,b=lrgb2oklab(*lin(h)); return L,math.hypot(a,b),math.degrees(math.atan2(b,a))%360
def sim(h,kind):
    M=MACHADO[kind]; r,g,b=lin(h)
    return [M[i][0]*r+M[i][1]*g+M[i][2]*b for i in range(3)]
def dE(c1,c2):
    a=lrgb2oklab(*c1); b=lrgb2oklab(*c2)
    return 100*math.dist(a,b)
def relL(h):
    r,g,b=lin(h); return 0.2126*r+0.7152*g+0.0722*b
def contrast(a,b):
    la,lb=relL(a),relL(b); hi,lo=max(la,lb),min(la,lb); return (hi+0.05)/(lo+0.05)
def run(pal,mode,surface):
    print(f"\n=== palette {pal}  mode={mode}  surface={surface} ===")
    lo,hi=BAND[mode]; fails=[]; warns=[]
    print(f"{'hex':<10}{'L':>7}{'C':>7}{'band':>8}{'chroma':>8}{'contrast':>10}")
    for h in pal:
        L,C,_=oklch(h); cr=contrast(h,surface)
        bok='ok' if lo<=L<=hi else 'FAIL'; cok='ok' if C>=CHROMA_FLOOR else 'FAIL'
        crs='ok' if cr>=CONTRAST_MIN else 'WARN'
        if bok=='FAIL': fails.append(f"{h} lightness {L:.3f} outside {lo}-{hi}")
        if cok=='FAIL': fails.append(f"{h} chroma {C:.3f} below {CHROMA_FLOOR}")
        if crs=='WARN': warns.append(f"{h} contrast {cr:.2f}:1 below {CONTRAST_MIN}")
        print(f"{h:<10}{L:>7.3f}{C:>7.3f}{bok:>8}{cok:>8}{cr:>9.2f}:1")
    print(f"\n{'pair':<22}{'normal':>9}{'protan':>9}{'deutan':>9}{'tritan':>9}{'verdict':>10}")
    for i in range(len(pal)-1):
        a,b=pal[i],pal[i+1]
        n=dE(lin(a),lin(b))
        p=dE(sim(a,'protan'),sim(b,'protan'))
        d=dE(sim(a,'deutan'),sim(b,'deutan'))
        t=dE(sim(a,'tritan'),sim(b,'tritan'))
        worst=min(p,d)
        v='ok' if worst>=CVD_TARGET else ('WARN' if worst>=CVD_FLOOR else 'FAIL')
        if v=='FAIL': fails.append(f"{a}/{b} CVD dE {worst:.1f} below floor {CVD_FLOOR}")
        if v=='WARN': warns.append(f"{a}/{b} CVD dE {worst:.1f} in the {CVD_FLOOR}-{CVD_TARGET} band")
        if n<NORMAL_FLOOR: fails.append(f"{a}/{b} normal-vision dE {n:.1f} below hard floor {NORMAL_FLOOR}")
        print(f"{a+' / '+b:<22}{n:>9.1f}{p:>9.1f}{d:>9.1f}{t:>9.1f}{v:>10}")
    print()
    for f in fails: print(f"  FAIL  {f}")
    for w in warns: print(f"  WARN  {w}")
    if not fails: print(f"  PASS  {len(pal)} slots, no hard failures" + (f", {len(warns)} warn(s)" if warns else ""))
    return len(fails)
if __name__=='__main__':
    # default: the report series' categorical pair, light then dark
    if len(sys.argv)>1:
        pal=sys.argv[1].split(','); mode=sys.argv[2] if len(sys.argv)>2 else 'light'
        surf=sys.argv[3] if len(sys.argv)>3 else ('#fcfcfb' if mode=='light' else '#1a1a19')
        sys.exit(1 if run(pal,mode,surf) else 0)
    n=run(['#2a78d6','#eb6834'],'light','#fcfcfb')
    n+=run(['#3987e5','#d95926'],'dark','#1a1a19')
    sys.exit(1 if n else 0)
