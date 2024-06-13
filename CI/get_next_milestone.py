#!/usr/bin/env python3

import os
import urllib.request
import json

GITHUB_TOKEN = os.environ["GH_TOKEN"]
GITHUB_REPO = os.environ["GH_REPO"]

with urllib.request.urlopen(
    urllib.request.Request(
        f"https://api.github.com/repos/{GITHUB_REPO}/milestones?state=open",
        headers={"Authorization": f"Bearer {GITHUB_TOKEN}"},
    )
) as response:
    milestones = json.loads(response.read())


titles = [m["title"] for m in milestones if m["title"] != "next"]
titles.sort()

print(titles[0])
