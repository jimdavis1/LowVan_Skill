import pymupdf, sys
def cols(page):
    d = page.get_text("dict")
    lines = []
    for b in d["blocks"]:
        for l in b.get("lines", []):
            txt = "".join(s["text"] for s in l["spans"]).rstrip()
            if txt.strip():
                lines.append((l["bbox"][0], l["bbox"][1], txt))
    if not lines:
        return ""
    mid = page.rect.width / 2
    L = sorted([x for x in lines if x[0] < mid - 10], key=lambda x: x[1])
    R = sorted([x for x in lines if x[0] >= mid - 10], key=lambda x: x[1])
    return "\n".join(t for _x, _y, t in L) + "\n" + "\n".join(t for _x, _y, t in R)
f = sys.argv[1]
d = pymupdf.open(f)
out = []
for i, p in enumerate(d):
    out.append("\n===== PAGE %d =====\n" % (i + 1) + cols(p))
open(sys.argv[2], "w").write("\n".join(out))
print("wrote", sys.argv[2], d.page_count, "pages")
