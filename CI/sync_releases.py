#!/usr/bin/env python3

import os
import json
import re
import argparse
import gitlab
import contextlib
import sys
import math

try:
    from halo import Halo
except:
    Halo = None


@contextlib.contextmanager
def spinner(text, *args, **kwargs):
    if sys.stdout.isatty() and Halo is not None:
        with Halo(text, *args, **kwargs):
            yield
    else:
        sys.stdout.write(text + "\n")
        yield


def parse_version(tag):
    return re.match(r"v(\d\.\d\d\.\d\d)", tag.name).group(1)


def group_items(labels, items):
    groups = {l: [] for l in labels}
    groups["Uncategorized"] = []

    for item in items:
        assigned = False
        for label in item.labels:
            if label in labels:
                groups[label].append(item)
                assigned = True
                break
        # is we get here, we didn't group this
        if not assigned:
            groups["Uncategorized"].append(item)
    return groups


def main():
    p = argparse.ArgumentParser()

    p.add_argument(
        "--access-token",
        help="Gitlab access token to update the releases",
        default=os.getenv("ATSJENKINS_ACCESS_TOKEN", None),
    )
    p.add_argument("--dry-run", "-s", action="store_true")
    p.add_argument("--verbose", "-v", action="store_true")

    label_groups = os.getenv(
        "RELEASE_NOTES_LABEL_GROUPS", "New Feature;Bug;Improvement"
    ).split(";")

    print("Label groups:", ", ".join(label_groups))

    args = p.parse_args()

    gl = gitlab.Gitlab("https://gitlab.cern.ch", private_token=args.access_token)
    if not args.dry_run:
        assert args.access_token is not None
        gl.auth()
    project = gl.projects.get("acts/acts-core")

    with spinner(text="Loading tags"):
        tags = project.tags.list(all=True)

    with spinner(text="Loading milestones"):
        milestones = project.milestones.list(all=True)
        ms_map = {}
        for ms in milestones:
            ms_map[ms.title] = ms

    for tag in tags:
        version = parse_version(tag)
        if not version in ms_map:
            print(f"No milestone found for tag f{tag.name} => skipping")
        milestone = ms_map[version]
        print(tag.name, milestone.title)
        # issues = list(milestone.issues())
        with spinner(text=f"Loading merge requests associated with %{milestone.iid}"):
            mrs = list(milestone.merge_requests())

        # need to get issues from merged MRs
        with spinner(text=f"Collecting issues from {len(mrs)} merged MRs"):
            issue_ids = []
            issues = []
            for mr in mrs:
                if mr.state != "merged":
                    continue
                for issue in mr.closes_issues():
                    if issue.id not in issue_ids:
                        issue_ids.append(issue.id)
                        issues.append(issue)

        issues_grouped = group_items(label_groups, issues)
        mrs_grouped = group_items(label_groups, mrs)

        if args.verbose:
            print("Issues:", ", ".join([str(i.iid) for i in issues]))

            for g, issues in issues_grouped.items():
                print(g, ", ".join([str(i.iid) for i in issues]))

            print("MRs:", ", ".join([str(mr.iid) for mr in mrs]))
            for g, mrs in mrs_grouped.items():
                print(g, ", ".join([str(mr.iid) for mr in mrs]))

        with spinner(text="Assembling release notes"):
            # make the Markdown
            md = ""
            # md = f"# Release {tag.name}\n"
            md += f"Milestone: [%{milestone.title}]({milestone.web_url})\n"

            if len(mrs) > 0:
                md += f"### {len(mrs)} Merge Requests in this release:\n"
                for g in label_groups + ["Uncategorized"]:
                    if len(mrs_grouped[g]) == 0:
                        continue
                    md += f"#### {g}\n"
                    for mr in mrs_grouped[g]:
                        md += f"- [!{mr.iid} - {mr.title}]({mr.web_url})\n"
                md += "\n"

            if len(issues) > 0:
                md += f"### {len(issues)} issues addressed in this release:\n"
                for g in label_groups + ["Uncategorized"]:
                    if len(issues_grouped[g]) == 0:
                        continue
                    md += f"#### {g}\n"
                    for issue in issues_grouped[g]:
                        md += f"- [#{issue.iid} - {issue.title}]({issue.web_url})\n"

        # print(md)
        if not args.dry_run:
            with spinner(text=f"Saving release notes on {tag.name}"):
                tag.set_release_description(md)
        if args.verbose:
            print("---")


if "__main__" == __name__:
    main()
