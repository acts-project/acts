#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess as sp
import requests
import tempfile
from urllib.parse import urljoin

# from fs.sshfs import SSHFS
from sshfs import SSHFS
from fs.osfs import OSFS
import fs.copy
import gitlab
import gitlab.exceptions
from datetime import datetime
from dateutil.parser import parse

from concurrent.futures import ThreadPoolExecutor, wait

import logging
# logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.DEBUG)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--deploy-user", default="atsjenkins")
    p.add_argument("--deploy-pwd", default=os.getenv("ATSJENKINS_PASSWORD"))
    p.add_argument("--coverage-source", required=True)
    p.add_argument("--commit-hash", default=os.getenv("CI_COMMIT_SHA"))
    p.add_argument("--coverage-commit-limit", default=int(os.getenv("COVERAGE_COMMIT_LIMIT", 10)), type=int)
    p.add_argument("--coverage-root", default=os.getenv("COVERAGE_WEBSITE_ROOT", "/eos/user/a/atsjenkins/www/ACTS/coverage"))
    p.add_argument("--website-public-url", default=os.getenv("COVERAGE_WEBSITE_URL", "https://acts.web.cern.ch/ACTS/coverage/"))
    p.add_argument("--project-id", default=3031, type=int)
    p.add_argument("--dry-run", "-s", action="store_true")

    args = p.parse_args()

    try:
        www_fs = SSHFS(host="lxplus.cern.ch",
                       user=args.deploy_user,
                       passwd=args.deploy_pwd).opendir(args.coverage_root)
        # www_fs = OSFS("www")
        listdir = www_fs.listdir(".")
    except:
        print("Unable to establish SSH connection to lxplus")
        print("This might indicate a problem with the credentials")
        print("or a temporary connection / configuration problem")

        raise
        sys.exit(1)

    gl = gitlab.Gitlab("https://gitlab.cern.ch/")
    project = gl.projects.get(args.project_id)

    commit_slug = args.commit_hash[:7]
    coverage_dest = os.path.join(args.coverage_root, commit_slug)
    print("Going to deploy coverage for", commit_slug, "to", coverage_dest)
    print("Will be publicly available under", urljoin(args.website_public_url, commit_slug))

    src_fs = OSFS(args.coverage_source)

    if not args.dry_run:
        fs.copy.copy_dir(src_fs, ".", www_fs, commit_slug)

    # cleanup
    # get all deployed commits
    deployed_commits = set(filter(www_fs.isdir, www_fs.listdir(".")))

    with ThreadPoolExecutor(max_workers=8) as tp:
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

    deployed_commits_with_time = list(reversed(sorted(deployed_commits_with_time, key=lambda i: i[1])))

    # take the n newest commits
    commits_to_keep = set(h for h,_ in deployed_commits_with_time[:args.coverage_commit_limit])

    print("Currently deployed commits:")
    for idx, (h, t) in enumerate(deployed_commits_with_time):
        if idx < args.coverage_commit_limit:
            print(" o", h, "-", t)
        else:
            print(" x", h, "-", t)


    print("Keeping commits:", ", ".join(commits_to_keep))

    commits_to_delete = deployed_commits - commits_to_keep

    if len(commits_to_delete) > 0:
        print("Removing:", ", ".join(commits_to_delete))

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
    """.format(latest_coverage_url)

    with www_fs.open("index.html", "w") as f:
        print("Writing index file redirecting to", latest_coverage_url)
        if not args.dry_run:
            f.write(index_content)

if "__main__" == __name__:
    main()
