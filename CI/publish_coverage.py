#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import argparse
import subprocess as sp
import requests
from urllib import quote_plus

p = argparse.ArgumentParser()
p.add_argument("--commit-hash", default=os.environ["CI_COMMIT_SHA"])

args = p.parse_args()

def ssh_cmd(cmd):
    return "ssh atsjenkins@lxplus.cern.ch -F {} \"{}\"".format(ssh_config_file, cmd)

COMMIT_HASH = args.commit_hash
coverage_commit_limit = 10
WEBSITE_ROOT = "/eos/user/a/atsjenkins/www/ACTS"
ssh_config_file = os.path.join(os.path.dirname(__file__), "ssh_config")
commit_slug = COMMIT_HASH[:7]
coverage_base = os.path.join(WEBSITE_ROOT, "coverage")
coverage_dest = os.path.join(coverage_base, commit_slug)
coverage_src = os.path.join(os.getcwd(), "build/coverage/")

print("Going to deploy coverage for", COMMIT_HASH, "to", coverage_dest)

if not os.path.exists(coverage_src):
    print("Coverage path does not exist")
    sys.exit(1)


mkdir_cmd = "ssh atsjenkins@lxplus.cern.ch -F {} \"mkdir -p {}\"".format(ssh_config_file, coverage_dest)
print(mkdir_cmd)
print(sp.check_output(mkdir_cmd, shell=True))

copy_cmd = "rsync -e \"ssh -F {}\" -ruv {} atsjenkins@lxplus.cern.ch:{}".format(ssh_config_file, coverage_src, coverage_dest)
print(copy_cmd)
print(sp.check_output(copy_cmd, shell=True))


# figure out what's deployed right now
ls_cmd = ssh_cmd("ls {}".format(coverage_base))
print(ls_cmd)
output = sp.check_output(ls_cmd, shell=True).decode("utf-8")
commits = output.strip().split("\n")

if len(commits) > coverage_commit_limit:
    base_url = "https://gitlab.cern.ch/api/v4/"
    def get_commit(project, commit):
        slug = quote_plus(project, safe="")
        url = base_url+"projects/"+slug+"/repository/commits/{}".format(commit)
        # print(url)
        res = requests.get(url, params={})
        data = res.json()
        return data


    commit_info = map(lambda c: get_commit("acts/acts-core", c), commits)

    commit_sorted = list(reversed(sorted(commit_info, key=lambda c: c["committed_date"])))

    print()

    for c in commit_sorted:
        print(c["id"], c["committed_date"])


    commits_to_delete = map(lambda c: c["id"], commit_sorted[coverage_commit_limit:])
    for c in commits_to_delete:
        commit_dir = c[:7]
        print("Deleting", commit_dir)
        assert commit_dir != ""
        rm_cmd = ssh_cmd("rm -rf {}".format(os.path.join(coverage_base, commit_dir)))
        print(rm_cmd)
        print(sp.check_outputr(rm_cmd, shell=True))
