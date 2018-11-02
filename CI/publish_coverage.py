#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess as sp
import requests
import tempfile
from urllib.parse import urljoin

from fs.sshfs import SSHFS
from fs.osfs import OSFS
import fs.copy
import gitlab
from datetime import datetime
from dateutil.parser import parse

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--deploy-user", default="atsjenkins")
    p.add_argument("--deploy-pwd", default=os.getenv("ATSJENKINS_PASSWORD"))
    p.add_argument("--coverage-source", required=True)
    p.add_argument("--commit-hash", default=os.getenv("CI_COMMIT_SHA"))
    p.add_argument("--coverage-commit-limit", default=int(os.getenv("COVERAGE_COMMIT_LIMIT", 10)), type=int)
    p.add_argument("--coverage-root", default=os.getenv("COVERAGE_WEBSITE_ROOT", "/eos/user/a/atsjenkins/www/ACTS/coverage"))
    p.add_argument("--website-public-url", default=os.getenv("COVERAGE_WEBSITE_URL", "https://acts.web.cern.ch/ACTS/coverage/"))

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
        sys.exit(1)

    gl = gitlab.Gitlab("https://gitlab.cern.ch/")
    project = gl.projects.get(3031) # acts/acts-core

    commit_slug = args.commit_hash[:7]
    coverage_dest = os.path.join(args.coverage_root, commit_slug)
    print("Going to deploy coverage for", commit_slug, "to", coverage_dest)
    print("Will be publicly available under", urljoin(args.website_public_url, commit_slug))

    src_fs = OSFS(args.coverage_source)

    fs.copy.copy_dir(src_fs, ".", www_fs, commit_slug)

    # cleanup
    # get all deployed commits
    deployed_commits = set(filter(www_fs.isdir, www_fs.listdir(".")))
    deployed_time = [parse(project.commits.get(c).committed_date) for c in deployed_commits]

    deployed_commits_with_time = list(reversed(sorted(zip(deployed_commits, deployed_time), key=lambda i: i[1])))

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
        f.write(index_content)

if "__main__" == __name__:
    main()
