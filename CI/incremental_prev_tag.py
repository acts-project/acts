#!/usr/bin/env python3

import os
import sys
from subprocess import check_output, run


def split_version(s):
    assert s.startswith("v")
    return tuple(int(c) for c in s[1:].split(".", 3))


tags_raw = check_output(["git", "tag", "-l"]).decode("utf8").strip().split("\n")

desc = run(
    ["git", "describe", "--exact-match", "HEAD"], capture_output=True, encoding="utf8"
)
on_tagged_commit = desc.returncode == 0

sys.stderr.write(f"On tagged commit? {on_tagged_commit}\n")

tags = []

for t in tags_raw:
    assert t.startswith("v")
    version = split_version(t)
    sys.stderr.write(f"{t} -> {version}\n")
    tags.append(version)

tags = list(sorted(tags))


if not on_tagged_commit:
    sys.stderr.write("Not on tagged commit, testing against highest tag\n")
    prev_tag = tags[-1]
else:
    current_tag = split_version(desc.stdout.strip())
    sys.stderr.write(f"Currently on {current_tag}")
    idx = tags.index(current_tag)
    assert idx != 0, "This doesn't work on the first ever tag"
    prev_tag = tags[idx - 1]

sys.stderr.write(f"Testing previous tag: v{'.'.join(map(str, prev_tag))}\n")

# this is the only thing that goes into stdout so we can catch it
sys.stdout.write("v%d.%02d.%02d" % prev_tag)
