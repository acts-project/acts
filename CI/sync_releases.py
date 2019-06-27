#!/usr/bin/env python3

import os
import json
import re
import argparse
import contextlib
import sys
import math

from util import def_arguments, Spinner, gitlab
from release_notes import (
    collect_milestone,
    make_release_notes,
    get_label_groups,
    parse_version,
)


def main():
    p = argparse.ArgumentParser()
    p = def_arguments(p, gl=True)

    p.add_argument("--dry-run", "-s", action="store_true")
    p.add_argument("--verbose", "-v", action="store_true")

    print("Label groups:", ", ".join(get_label_groups()))

    args = p.parse_args()

    gl = gitlab(args)

    project = gl.projects.get("acts/acts-core")

    with Spinner(text="Loading tags"):
        tags = project.tags.list(all=True)

    with Spinner(text="Loading milestones"):
        milestones = project.milestones.list(all=True)
        ms_map = {}
        for ms in milestones:
            ms_map[ms.title] = ms

    for tag in tags:
        version = parse_version(tag.name)
        if not version in ms_map:
            print(f"No milestone found for tag f{tag.name} => skipping")
        milestone = ms_map[version]
        print(tag.name, milestone.title)

        mrs_grouped, issues_grouped = collect_milestone(milestone)

        if args.verbose:
            print("Issues:", ", ".join([str(i.iid) for i in issues]))

            for g, issues in issues_grouped.items():
                print(g, ", ".join([str(i.iid) for i in issues]))

            print("MRs:", ", ".join([str(mr.iid) for mr in mrs]))
            for g, mrs in mrs_grouped.items():
                print(g, ", ".join([str(mr.iid) for mr in mrs]))

        with Spinner(text="Assembling release notes"):
            md = make_release_notes(milestone, mrs_grouped, issues_grouped)

        # print(md)
        if not args.dry_run:
            with Spinner(text=f"Saving release notes on {tag.name}"):
                tag.set_release_description(md)
        if args.verbose:
            print("---")


if "__main__" == __name__:
    main()
