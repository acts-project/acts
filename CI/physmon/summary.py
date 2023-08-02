#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import functools
import os


parser = argparse.ArgumentParser()
parser.add_argument("inputs", nargs="+")
parser.add_argument("--base", required=True)
parser.add_argument("--html")
parser.add_argument("--md")
args = parser.parse_args()

re_title = re.compile(r'<p class="title">\s*(.*)\s*<\/p>', re.RegexFlag.MULTILINE)
re_check = re.compile(r'<a.*title="(.*)">\s*(.)\s*<\/a>', re.RegexFlag.MULTILINE)

summary = {}

for h in args.inputs:
    with open(h, mode="r", encoding="utf-8") as f:
        try:
            content = f.read()
            print(h, re_title.findall(content))
            title = re_title.findall(content)[0]
            checks = re_check.findall(content)
            parsed_checks = list(map(lambda c: c[1] == "✅", checks))
            summary[h] = {
                "title": title,
                "checks": checks,
                "parsed_checks": parsed_checks,
                "total": functools.reduce(lambda a, b: a and b, parsed_checks),
            }
        except Exception as e:
            print(r"could not parse {h}", e)

output = "# physmon summary\n"
for h, s in summary.items():
    path = os.path.relpath(h, args.base)
    output += f"  - {'✅' if s['total'] else '🔴'} [{s['title']}]({path})\n"

if args.html:
    with open(args.html, mode="w", encoding="utf-8") as f:
        f.write(
            """
<!DOCTYPE html>
<script>window.texme = { style: 'plain' }</script>
<script src="https://cdn.jsdelivr.net/npm/texme@1.2.2"></script><textarea>

                """
        )
        f.write(output)

if args.md:
    with open(args.md, mode="w", encoding="utf-8") as f:
        f.write(output)
