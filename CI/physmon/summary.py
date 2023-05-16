#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import functools
import os


parser = argparse.ArgumentParser()
parser.add_argument("html", nargs="+")
parser.add_argument("output")
args = parser.parse_args()

re_title = re.compile(r'<p class="title">\s*(.*)\s*<\/p>', re.RegexFlag.MULTILINE)
re_check = re.compile(r'<a.*title="(.*)">\s*(.)\s*<\/a>', re.RegexFlag.MULTILINE)

summary = {}

for h in args.html:
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

with open(args.output, mode="w", encoding="utf-8") as f:
    f.write(
        """<!DOCTYPE html>
<html>
<head>
  <title>physmon summary</title>
  <meta charset="UTF-8">
</head>
<body>
  <h1>physmon summary</h1>
  <ul>"""
    )

    for h, s in summary.items():
        path = os.path.relpath(h, os.path.dirname(args.output))
        f.write(
            f"""
    <li>{"✅" if s["total"] else "🔴"} <a href="{path}">{s["title"]}</a></li>"""
        )

    f.write(
        """
  </ul>
</body>
</html>"""
    )
