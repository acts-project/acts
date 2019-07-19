#!/usr/bin/env python
from __future__ import print_function

import os
import sys
from subprocess import check_output


def split_version(s):
    assert s.startswith("v")
    return tuple(int(c) for c in s[1:].split(".", 3))


tags_raw = check_output(["git", "tag", "-l"]).decode("utf8").strip().split("\n")

try:
    desc = check_output(["git", "describe", "--exact-match", "HEAD"]).strip()
    on_tagged_commit = True
except Exception as e:
    on_tagged_commit = False

sys.stderr.write("On tagged commit? %s\n" % on_tagged_commit)

tags = []

for t in tags_raw:
    assert t.startswith("v")
    version = split_version(t)
    sys.stderr.write("%s -> %s\n" % (t, version))
    tags.append(version)

tags = list(sorted(tags))


if not on_tagged_commit:
    sys.stderr.write("Not on tagged commit, testing against highest tag\n")
    prev_tag = tags[-1]
else:
    current_tag = split_version(desc)
    sys.stderr.write("Currently on %s\n" % str(current_tag))
    idx = tags.index(current_tag)
    assert idx != 0, "This doesn't work on the first ever tag"
    prev_tag = tags[idx - 1]

sys.stderr.write("Testing previous tag: v%s\n" % ".".join(map(str, prev_tag)))

# this is the only thing that goes into stdout so we can catch it
sys.stdout.write("v%d.%02d.%02d" % prev_tag)
