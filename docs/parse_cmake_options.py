#!/usr/bin/env python3
"""
This script accepts a cmake lists file as an argument extracts all
`option` and `set(... CACHE ...)` variables. It then writes a
markdown table to stdout
"""

import argparse
from pathlib import Path
import re
from tabulate import tabulate
import textwrap

p = argparse.ArgumentParser(description=__doc__)

p.add_argument("cmakefile", help="Input cmake lists file to parse")
p.add_argument(
    "--prefix", default="ACTS_", help="Prefix to identify relevant variables to extract"
)
p.add_argument(
    "--width",
    type=int,
    default=40,
    help="Width of second column generated from cmake doc strings",
)

args = p.parse_args()

cmakefile = Path(args.cmakefile)

with cmakefile.open() as fh:
    rows = []
    for line in fh:
        if m := re.match(rf"option\( *({args.prefix}\w*) \"(.*)\" (ON|OFF)\ *\)", line):
            name, doc, default = m.groups()
            type = "bool"
        elif m := re.match(
            rf"set\( *({args.prefix}\w*) \"(.*)\" CACHE (\w+) \"(.*)\"( FORCE)? *\)",
            line,
        ):
            name, default, type, doc, _ = m.groups()
            type = type.lower()
            if default == "":
                default = '""'
        else:
            continue
        doc = "<br>".join(textwrap.wrap(doc, width=args.width))
        rows.append((name, f"{doc}<br> type: `{type}`, default: `{default}`"))

print(tabulate(rows, headers=("Option", "Description"), tablefmt="github"))
