#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import functools
import os

HERALD_URL = "https://herald.dokku.paulgessinger.com/view/{repo}/runs/{run_id}/artifacts/{artifact_name}/{path}"
IS_CI = "GITHUB_ACTIONS" in os.environ


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

if args.html:
    with open(args.html, mode="w", encoding="utf-8") as f:
        f.write(
            """<!DOCTYPE html>
<html>
<head>
  <title>physmon summary</title>
  <meta charset="UTF-8">
</head>
<body>
  <h1>physmon summary</h1>
  <ul>
            """
        )

        for h, s in summary.items():
            path = os.path.relpath(h, args.base)
            f.write(
                f"""
        <li>{"✅" if s["total"] else "🔴"} <a href="{path}">{s["title"]}</a></li>"""
            )

        f.write(
            """
      </ul>
    </body>
    </html>
            """
        )

if args.md:
    with open(args.md, mode="w", encoding="utf-8") as f:
        f.write("# physmon summary\n")
        for h, s in summary.items():
            path = os.path.relpath(h, args.base)
            if IS_CI:
                url = HERALD_URL.format(
                    repo=os.environ["GITHUB_REPOSITORY"],
                    run_id=os.environ["GITHUB_RUN_ID"],
                    artifact_name="physmon",
                    path=path,
                )
            else:
                url = path
            f.write(f"  - {'✅' if s['total'] else '🔴'} [{s['title']}]({url})\n")
