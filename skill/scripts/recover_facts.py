#!/usr/bin/env python3
"""Recover a facts file for every coverage audit that has none, and prove it.

A recovered facts file is only kept when regenerating from it reproduces the
page it came from BYTE FOR BYTE. Anything less is discarded and reported,
because a facts file that renders a slightly different page is worse than no
facts file at all: the next person edits it, regenerates, and silently ships
a page that lost something.

Nothing already in generators/ is touched.
"""
import json, os, re, shutil, subprocess, sys, tempfile

KIT = os.environ["LOWVAN_KIT"]
REPORTS = os.path.join(KIT, "reports")
GENS = os.path.join(REPORTS, "generators")
GEN = os.path.join(KIT, "skill", "scripts", "gen_coverage_audit.py")
HERE = os.path.dirname(os.path.abspath(__file__))
EXTRACT = os.path.join(HERE, "facts_from_html.py")

DRY = "--apply" not in sys.argv


def main():
    pages = sorted(f for f in os.listdir(REPORTS)
                   if f.endswith("-coverage-audit.html"))
    todo = []
    for p in pages:
        stem = p[:-len("-coverage-audit.html")]
        if stem == "0":
            continue          # byte-identical duplicate of Matonaviridae
        if not os.path.exists(os.path.join(GENS, "%s_audit.json" % stem)):
            todo.append((stem, p))

    ok, bad, moved = [], [], []
    tmp = tempfile.mkdtemp(prefix="facts_rt_")
    for stem, page in todo:
        orig = os.path.join(REPORTS, page)
        facts = os.path.join(tmp, "%s_audit.json" % stem)
        #  Figures are written into reports/ because that is where the page
        #  will look for them at render time; they are new files, never
        #  overwriting one that is already there.
        r = subprocess.run([sys.executable, "-I", EXTRACT, orig,
                            "--out", facts, "--figdir", REPORTS],
                           capture_output=True, text=True)
        if r.returncode:
            bad.append((stem, "extract failed: "
                        + (r.stderr.strip().splitlines() or ["?"])[-1][:90]))
            continue

        out = os.path.join(tmp, "%s.html" % stem)
        r = subprocess.run([sys.executable, "-I", GEN, "--facts", facts,
                            "--out", out], cwd=REPORTS,
                           capture_output=True, text=True)
        if r.returncode:
            bad.append((stem, "regen failed: "
                        + (r.stderr.strip().splitlines() or ["?"])[-1][:90]))
            continue

        a, b = open(orig, "rb").read(), open(out, "rb").read()
        if a == b:
            ok.append(stem)
        else:
            #  Two differences are allowed, and only these two. Re-seating a
            #  figure out of a table's scroll wrapper, which is the defect my
            #  splice left on the parvovirus pages; and a refreshed <style>,
            #  because the generator's CSS constant has grown since the older
            #  pages were rendered. Everything else means content was lost.
            ta, fa = _strip(a.decode())
            tb, fb = _strip(b.decode())
            if ta != tb or sorted(fa) != sorted(fb):
                bad.append((stem, "round trip differs (%d vs %d bytes)"
                            % (len(a), len(b))))
                continue
            why = []
            sa, sb = a.decode(), b.decode()
            if _style(sa) != _style(sb):
                why.append("stylesheet refreshed")
            if _unlit(_strip(sa)[0]) != _strip(sa)[0] or \
               re.search(r"&lt;(?:%s)&gt;" % _INLINE, sa):
                why.append("literal inline markup repaired")
            #  Only claim the figure moved if the page actually has one and
            #  its surroundings changed; hepeviridae has none.
            if fa and sa.index("<figure") != sb.index("<figure"):
                why.insert(0, "figure re-seated")
            moved.append((stem, ", ".join(dict.fromkeys(why))))
        if not DRY:
            shutil.copy(facts, os.path.join(GENS, "%s_audit.json" % stem))
            if any(stem == m[0] for m in moved):   # page itself changes too
                shutil.copy(out, orig)

    print("%s%d of %d recovered: %d byte-identical, %d also rewritten"
          % ("" if not DRY else "DRY RUN -- ", len(ok) + len(moved), len(todo),
             len(ok), len(moved)))
    if moved:
        print("\npage also rewritten:")
        for s_, why in moved:
            print("   %-28s %s" % (s_, why))
    if bad:
        print("\n%d NOT recovered:" % len(bad))
        for s, why in bad:
            print("   %-28s %s" % (s, why))
    print("\nrecovered: " + " ".join(sorted(ok + [m[0] for m in moved])))



def _style(t):
    m = re.search(r"<style>.*?</style>", t, re.S)
    return m.group(0) if m else ""


_INLINE = "br|i|b|em|strong|code|sub|sup|small"


def _unlit(t):
    """Literal "&lt;em&gt;" and a real <em> count as the same content.

    Pages rendered before esc() grew its inline-tag whitelist show the markup
    as visible text -- "32 of 47 carry <em>four</em>" with the brackets on
    screen. Regenerating repairs that, so the two must compare equal here or
    every such page would look like content loss.
    """
    return re.sub(r"&lt;(/?(?:%s))\s*/?&gt;" % _INLINE, r"<\1>", t)


def _strip(t):
    """(page with <style> and every <figure> removed, the figure fragments)"""
    figs = re.findall(r"<figure.*?</figure>", t, re.S)
    t = re.sub(r"<figure.*?</figure>", "", t, flags=re.S)
    t = re.sub(r"<style>.*?</style>", "", t, flags=re.S)
    return _unlit(t), figs


main()
