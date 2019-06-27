#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess as sp
import requests
import tempfile
from urllib.parse import urljoin

from fs.osfs import OSFS
import fs.copy
import gitlab.exceptions
from datetime import datetime
from dateutil.parser import parse

from concurrent.futures import ThreadPoolExecutor, wait

from util import get_lxplus_fs, def_arguments, Spinner, gitlab

import logging

# logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)


def main():
    p = argparse.ArgumentParser()
    p = def_arguments(p, acc=True, gl=True)
    p.add_argument("--coverage-source", required=True)
    p.add_argument(
        "--ref", default=os.getenv("CI_COMMIT_TAG", os.getenv("CI_COMMIT_SHA", None))
    )
    p.add_argument(
        "--coverage-commit-limit",
        default=int(os.getenv("COVERAGE_COMMIT_LIMIT", 10)),
        type=int,
    )
    p.add_argument(
        "--coverage-root",
        default=os.getenv(
            "COVERAGE_WEBSITE_ROOT", "/eos/user/a/atsjenkins/www/ACTS/coverage"
        ),
    )
    p.add_argument(
        "--website-public-url",
        default=os.getenv(
            "COVERAGE_WEBSITE_URL", "https://acts.web.cern.ch/ACTS/coverage/"
        ),
    )
    p.add_argument("--project-id", default=3031, type=int)
    p.add_argument("--dry-run", "-s", action="store_true")

    args = p.parse_args()

    try:
        www_fs = get_lxplus_fs(args).opendir(args.coverage_root)
        # www_fs = OSFS("www")
        listdir = www_fs.listdir(".")
    except:
        print("Unable to establish SSH connection to lxplus")
        print("This might indicate a problem with the credentials")
        print("or a temporary connection / configuration problem")

        raise
        sys.exit(1)

    gl = gitlab(args)
    project = gl.projects.get(args.project_id)

    if len(args.ref) == 40:
        # is commit hash
        deploy_name = args.ref[:8]
    else:
        # probably tag
        deploy_name = args.ref

    coverage_dest = os.path.join(args.coverage_root, deploy_name)
    print("Going to deploy coverage for", deploy_name, "to", coverage_dest)
    print(
        "Will be publicly available under",
        urljoin(args.website_public_url, deploy_name),
    )

    src_fs = OSFS(args.coverage_source)

    with Spinner(f"Publishing ref {deploy_name}"):
        if not args.dry_run:
            fs.copy.copy_dir(src_fs, ".", www_fs, deploy_name)

    # cleanup
    # get all deployed commits
    with Spinner(text="Getting deployed commits"):
        deployed_commits = set()
        for item in www_fs.listdir("."):
            if not www_fs.isdir(item):
                continue
            if item.startswith("v"):  # skip versions
                continue
            deployed_commits.add(item)

    with Spinner(text="Getting info for deployed commits"):
        with ThreadPoolExecutor(max_workers=20) as tp:
            # deployed_commit_info = p.map(project.commits.get, deployed_commits)
            futures = [tp.submit(project.commits.get, c) for c in deployed_commits]
            wait(futures)

    deployed_commits_with_time = []
    for commit, future in zip(deployed_commits, futures):
        try:
            info = future.result()
            date = parse(info.committed_date)
            deployed_commits_with_time.append((commit, date))
        except gitlab.exceptions.GitlabGetError as e:
            print("Commit", commit, "not found, will remove")

    deployed_commits_with_time = list(
        reversed(sorted(deployed_commits_with_time, key=lambda i: i[1]))
    )

    # take the n newest commits
    commits_to_keep = set(
        h for h, _ in deployed_commits_with_time[: args.coverage_commit_limit]
    )

    print("Currently deployed commits:")
    for idx, (h, t) in enumerate(deployed_commits_with_time):
        if idx < args.coverage_commit_limit:
            print(" o", h, "-", t)
        else:
            print(" x", h, "-", t)

    print("Keeping commits:", ", ".join(commits_to_keep))

    commits_to_delete = deployed_commits - commits_to_keep

    if len(commits_to_delete) > 0:
        with Spinner("Removing: %s" % ", ".join(commits_to_delete)):
            if not args.dry_run:
                for commit in commits_to_delete:
                    www_fs.removetree(commit)

    # install / update indexfile
    latest_commit = deployed_commits_with_time[0][0]
    latest_coverage_url = urljoin(args.website_public_url, latest_commit)
    index_content = """
<!DOCTYPE html>
<html>
<head>
<meta http-equiv="refresh" content="0; url={0}" />
</head>
<body>
Redirecting to <a href"{0}">{0}</a>
</body>
</html>
    """.format(
        latest_coverage_url
    )

    with Spinner("Writing index file redirecting to %s" % latest_coverage_url):
        if not args.dry_run:
            with www_fs.open("index.html", "w") as f:
                f.write(index_content)


if "__main__" == __name__:
    main()
