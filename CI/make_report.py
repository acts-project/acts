#!/usr/bin/env python3

import argparse
import re
from collections import namedtuple
from itertools import groupby
import os
import html
from fnmatch import fnmatch
import json


from codereport import CodeReport, ReportItem

def parse_clang_tidy_item(itemstr):

    m = re.match(r"([/.\-\w]+):(\d+):(\d+): (.*?):(.*?)\n((?:.|\n)*)", itemstr)
    # keyword / category
    mkey = re.match(r"(.*)\[(.*)\]$", m.group(5))
    if mkey is None:
        code = "unknown"
    else:
        code = mkey.group(2)
    msg = mkey.group(1).strip() +"\n"+ m.group(6)

    if m is None:
        print(itemstr)

    item = ReportItem(
        path=m.group(1),
        line=int(m.group(2)),
        col=int(m.group(3)),
        message=msg,
        code=code,
        severity=m.group(4)
    )

    return item

def parse_clang_tidy_output(output):

    # cleanup
    itemstr = re.sub(r"Enabled checks:\n(?:.|\n)+\n\n", "", output)
    itemstr = re.sub(r"clang-tidy-\d\.\d.*\n?", "", itemstr)
    itemstr = re.sub(r".*-header-filter.*", "", itemstr)

    items = []
    prevstart = 0

    matches = list(re.finditer(r"([\w/.\-]+):(\d+):(\d+): warning:", itemstr))
    for idx, m in enumerate(matches):
        # print(m)
        start, end = m.span()
        if idx > 0:
            item = itemstr[prevstart:start]
            items.append(item)
        if idx+1 == len(matches):
            item = itemstr[start:]
            items.append(item)
        prevstart = start

    items = set(map(parse_clang_tidy_item, sorted(items)))

    return items


def main():
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=("clang-tidy"))
    p.add_argument("inputfile")
    p.add_argument("output", default="codereport_clang_tidy.json")
    p.add_argument("--exclude", "-e", action="append", default=[])
    p.add_argument("--filter", action="append", default=[])

    args = p.parse_args()

    if args.mode == "clang-tidy":
        with open(args.inputfile, "r", encoding="utf-8") as f:
            inputstr = f.read()
        items = parse_clang_tidy_output(inputstr)

        def select(item):
            accept = True
            if len(args.filter) > 0:
                accept = accept and all(fnmatch(item.path, e) for e in args.filter)

            accept = accept and not any(fnmatch(item.path, e) for e in args.exclude)
            return accept

        items = filter(select, items)


        data = [i.dict() for i in items]
        with open(args.output, "w+") as jf:
            json.dump(data, jf, indent=2)


if "__main__" == __name__:
    main()


