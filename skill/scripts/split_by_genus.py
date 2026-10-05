#!/usr/bin/env python3
"""Move one genus's sequences out of a shared feature into a key of its own.

A slot name is a homology group inside a genus and not across one. Where the
per-genus medians inside a grouped module differ by more than about twofold,
the slot is holding two different proteins and a single profile set would
blur them -- the same trap that keeps Cardiovirus and Aphthovirus apart
despite both declaring an "L".

    split_by_genus.py <GROUP> <KEY> <NEWKEY|DROP> <Genus>[,<Genus>...]
"""
import collections, os, sys

S = os.path.expanduser("~/lowvan_scratch/picorna")
group, key, newkey, genera = sys.argv[1], sys.argv[2], sys.argv[3], set(sys.argv[4].split(","))
sp = {}
for l in open(f"{S}/meta/{group}.tsv"):
    f = l.rstrip("\n").split("\t"); sp[f[0]] = f[3]
src = f"{S}/collections/{group}/{key}.fasta"
keep, move = [], []
h = None
for l in open(src):
    if l.startswith(">"): h = l.rstrip()
    else:
        gid = h[1:].split()[0].split("|")[1].rsplit(".", 2)[0]
        (move if sp.get(gid) in genera else keep).append((h, l.strip()))
with open(src, "w") as f:
    for h, s in keep: f.write("%s\n%s\n" % (h, s))
if newkey != "DROP":
    with open(f"{S}/collections/{group}/{newkey}.fasta", "w") as f:
        for h, s in move: f.write("%s\n%s\n" % (h, s))
def med(v): 
    L = sorted(len(s) for _, s in v); return (len(L), L[len(L)//2] if L else 0)
print("  %-16s %-5s kept n=%d median %d   %s n=%d median %d"
      % (group, key, *med(keep), newkey, *med(move)))
