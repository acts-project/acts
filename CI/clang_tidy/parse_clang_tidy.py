#!/usr/bin/env python3

"""
This produces a structured normalized report json file based on warnings generated by another tool.
Currently implemented is clang-tidy warnings.
"""

import argparse
import re
from collections import namedtuple
from itertools import groupby
import os
import html
from fnmatch import fnmatch
import json
import sys
from dataclasses import dataclass
from pathlib import Path

from item import Item, ItemCollection


def parse_clang_tidy_item(itemstr):

    try:
        m = re.match(
            r"(?P<file>[/.\-+\w]+):(?P<line>\d+):(?P<col>\d+): (?P<sev>.*?):(?P<msg>[\s\S]*)\[(?P<code>.*)\]\n(?P<info>[\s\S]*)",
            itemstr,
        )

        lines = itemstr.split("\n")

        item = Item(
            path=Path(m.group("file")),
            line=int(m.group("line")),
            col=int(m.group("col")),
            #  message=m.group("msg").strip(),
            message=m.group("msg").strip() + "\n" + "\n".join(lines[1:]),
            code=m.group("code"),
            severity=m.group("sev"),
        )

        #  print(repr(item))

        return item
    except:
        print("Failed parsing clang-tidy item:")
        print("-" * 20)
        print(itemstr)
        print("-" * 20)
        raise


def parse_clang_tidy_output(output):

    # cleanup
    itemstr = output
    itemstr = re.sub(r"Enabled checks:\n[\S\s]+?\n\n", "", itemstr)
    itemstr = re.sub(r"clang-tidy-\d\.\d.*\n?", "", itemstr)
    itemstr = re.sub(r"clang-apply-.*", "", itemstr)
    itemstr = re.sub(r".*-header-filter.*", "", itemstr)

    items = []
    prevstart = 0

    matches = list(
        re.finditer(r"([\w/.\-+]+):(\d+):(\d+): (?:(?:warning)|(?:error)):", itemstr)
    )
    for idx, m in enumerate(matches):
        # print(m)
        start, end = m.span()
        if idx > 0:
            item = itemstr[prevstart:start]
            items.append(item)
        if idx + 1 == len(matches):
            item = itemstr[start:]
            items.append(item)
        prevstart = start

    items = set(map(parse_clang_tidy_item, sorted(items)))

    return items


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("inputfile", help="The input file containing the warnings")
    p.add_argument(
        "output", default="codereport_clang_tidy.json", help="The resulting JSON file"
    )
    p.add_argument(
        "--exclude",
        "-e",
        action="append",
        default=[],
        help="Exclude files that match any of these patterns",
    )
    p.add_argument(
        "--filter",
        action="append",
        default=[],
        help="Only include files that match any of these patterns",
    )
    p.add_argument(
        "--ignore",
        action="append",
        default=[],
        help="Ignore items with codes matching any of these patterns",
    )
    p.add_argument("--cwd", type=Path)
    p.add_argument("--strip-common", action="store_true")

    args = p.parse_args()

    with open(args.inputfile, "r", encoding="utf-8") as f:
        inputstr = f.read()
    items = parse_clang_tidy_output(inputstr)

    def select(item):
        accept = True
        if len(args.filter) > 0:
            accept = accept and all(fnmatch(item.path, e) for e in args.filter)

        accept = accept and not any(fnmatch(item.path, e) for e in args.exclude)

        accept = accept and not any(fnmatch(item.code, i) for i in args.ignore)

        return accept

    items = list(filter(select, items))

    if args.cwd:
        for item in items:
            item.path = (args.cwd / item.path).resolve()

    if args.strip_common:
        prefix = Path(os.path.commonprefix([i.path for i in items]))

        def subpath(m):
            path, line, col = m.groups()
            path = Path(path).resolve().relative_to(prefix)
            return f"{path}:{line}:{col}:"

        for item in items:
            item.path = item.path.relative_to(prefix)

            item.message = re.sub(r"([\w/.\-+]+):(\d+):(\d+):", subpath, item.message)

    print("Write to", args.output)
    with open(args.output, "w+") as jf:
        jf.write(ItemCollection(pydantic.RootModel=items).json(indent=2))


if "__main__" == __name__:
    main()
