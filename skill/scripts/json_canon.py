#!/usr/bin/env python3
"""One canonical JSON serialisation, so a diff shows only real changes.

`Viral_PSSM.json` is written by tools in two languages. Without a single agreed
format every write reformats the whole file and buries the actual change. Worse,
a plain Perl round-trip is not even stable *with itself*: Perl randomises hash
order, so decoding and re-encoding the master reorders all 7,373 lines
differently on every run. That is the bug `JSON::XS ->canonical` exists to fix.

The canonical form here is exactly what Perl produces with

    JSON::XS->new->pretty->canonical->encode($data)

which is: 3-space indent, `" : "` between key and value, keys sorted at every
level, and a trailing newline. Indent width is fixed at 3 in JSON::XS and is not
configurable, so that is the number, not a preference.

Python matches it byte for byte given one normalisation: an integral float is
written `1`, not `1.0`, because Perl cannot tell the two apart and would write
`1`. That is the only place the two languages disagree on this schema.

    python3 json_canon.py --check FILE...     is it canonical? (exit 1 if not)
    python3 json_canon.py --write FILE...     rewrite in place

From another script:

    from json_canon import dumps, write_json
    write_json(obj, path)

The equivalent Perl, for anything on that side of the fence:

    perl -MJSON::XS -e 'local $/; my $d = JSON::XS->new->decode(<>);
        print JSON::XS->new->pretty->canonical->encode($d);' in.json > out.json
"""

import argparse
import json
import os
import sys

# JSON::XS ->pretty is indent(3) + space_before + space_after
INDENT = 3
SEPARATORS = (",", " : ")


def normalize(obj):
    """Make the structure encode identically in Perl and Python.

    Perl has no separate float type for an integral value, so JSON::XS writes
    1.0 as 1. Fold those down before encoding or the two languages produce
    different bytes for the same data.
    """
    if isinstance(obj, bool):
        return obj
    if isinstance(obj, float) and obj.is_integer():
        return int(obj)
    if isinstance(obj, dict):
        return {k: normalize(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [normalize(v) for v in obj]
    return obj


def dumps(obj):
    """Canonical text, byte-identical to JSON::XS ->pretty->canonical."""
    return json.dumps(normalize(obj), indent=INDENT, separators=SEPARATORS,
                      sort_keys=True, ensure_ascii=False) + "\n"


def write_json(obj, path):
    """Write canonical JSON to path. Returns True if the bytes changed."""
    new = dumps(obj)
    old = None
    if os.path.exists(path):
        with open(path, encoding="utf8") as fh:
            old = fh.read()
    if old == new:
        return False
    with open(path, "w", encoding="utf8") as fh:
        fh.write(new)
    return True


def canonical_text(path):
    with open(path, encoding="utf8") as fh:
        return dumps(json.load(fh))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--check", action="store_true", help="report, change nothing")
    g.add_argument("--write", action="store_true", help="rewrite in place")
    args = ap.parse_args()

    bad = 0
    for path in args.files:
        try:
            want = canonical_text(path)
        except Exception as exc:
            print("  %-46s UNREADABLE: %s" % (path, exc))
            bad += 1
            continue
        with open(path, encoding="utf8") as fh:
            have = fh.read()
        if have == want:
            print("  %-46s canonical" % path)
            continue
        bad += 1
        if args.write:
            with open(path, "w", encoding="utf8") as fh:
                fh.write(want)
            print("  %-46s REWRITTEN (%d -> %d bytes)"
                  % (path, len(have.encode()), len(want.encode())))
        else:
            # say how much of the difference is only formatting
            same_data = json.loads(have) == json.loads(want)
            print("  %-46s NOT canonical (%s)"
                  % (path, "formatting only" if same_data else "CONTENT DIFFERS"))

    if bad and args.check:
        print("\n  %d file(s) not canonical. Rerun with --write." % bad)
    return 1 if bad and args.check else 0


if __name__ == "__main__":
    sys.exit(main())
