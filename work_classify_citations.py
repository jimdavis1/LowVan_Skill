import json, glob, collections
from collections import OrderedDict
#  coordinates: the paper fixes a BOUNDARY -- a cleavage site, a start codon,
#  a frameshift or readthrough position, a splice junction, a domain span, or
#  the genome sequence/organisation that places the ORF.
#  function:    the paper says what the protein DOES.
COORD = {
 "26296665","9191923","1413521","8337841","11447125","11350035","15755749",
 "1518855","23255617","31399875","2214022","2273395","10799588","10823845",
 "6897030","21507977","4038520","8445732","19196463","19906906","33494395",
 "2174669","12827468","6964389","16453503","16453524","18822126","7241658",
 "1985194","3553612","8627739","6322438","12466484","12021349","2141206",
 "12917405","9125055",
}
FUNC = {
 "1853569","32525472","24672528","11722009","15098118","12771402","24257609",
 "8607260","38819305","26700068","36852909","19622744","28096411","34523963",
 "11160697","39207105","9557742","35723908","8683214","28276468","7753193",
 "2769760","28786782","17794341","7968923","9094737","11301009","2200886",
 "23035234","10357827","36792007","30825649","19036396","21430054","23864626",
 "21124458","23720714","36399124","16847135",
}
#  Verified against the fetched title: these do not support the feature at all.
WRONG = {("Allexivirus","ALLEXI_40K","11752591"),
         ("Orthoflavivirus","NS3","8107205"),
         ("Orthoflavivirus","NS5","2170676"),
         ("Orthoflavivirus","PRM","18184701")}

stat=collections.Counter(); removed=[]
for p in sorted(glob.glob("modules/*/*_Viral_PSSM.json")):
    doc=json.load(open(p),object_pairs_hook=OrderedDict); touched=False
    for mod,blk in doc.items():
        if not isinstance(blk,dict): continue
        f=blk.get("features",blk)
        for k,v in (f.items() if isinstance(f,dict) else []):
            if not isinstance(v,dict): continue
            cg=v.get("PMID_claude_generated")
            ids = cg.get("unclassified",[]) if isinstance(cg,dict) else (cg or [])
            if not ids: continue
            out=OrderedDict(); leftover=[]
            for i in ids:
                if (mod,k,i) in WRONG:
                    removed.append((mod,k,i)); stat["removed"]+=1; continue
                cat = "coordinates" if i in COORD else "function" if i in FUNC else None
                if cat is None: leftover.append(i); stat["unclassified"]+=1; continue
                out.setdefault(cat,[]).append(i); stat[cat]+=1
            if leftover: out["unclassified"]=leftover
            if out: v["PMID_claude_generated"]=out
            else: v.pop("PMID_claude_generated",None)
            touched=True
    if touched: json.dump(doc,open(p,"w"),indent=1)
print("  coordinates %d   function %d   still unclassified %d   removed %d"
      %(stat["coordinates"],stat["function"],stat["unclassified"],stat["removed"]))
for r in removed: print("     removed %s %s -> PMID %s"%r)
