#!/usr/bin/env python3

import argparse
import yaml
import os
import json
import sys
from fnmatch import fnmatch
import codecs
from tabulate import tabulate
from operator import itemgetter
from gitlab import Gitlab
import yaml
from sshfs import SSHFS
from fs.osfs import OSFS
import fs.copy
from urllib.parse import urljoin

from codereport import CodeReport, ReportItem, report_from_json
from static_analysis_results import analysis

def publish(srcfs, destfs, prefix):
    destfs.makedirs(prefix, recreate=True)
    targetdir = destfs.opendir(prefix)
    # clear it
    print("Publishing to", targetdir.geturl("/"))
    targetdir.removetree("/")
    fs.copy.copy_dir(src_fs=srcfs, src_path="/", 
                     dst_fs=targetdir, dst_path="/")
    

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--limitfile", required=True)
    p.add_argument("--results", required=True)
    p.add_argument("--report-dir", required=True)
    p.add_argument("--gitlab-token", required=True)
    p.add_argument("--ref", required=True)
    p.add_argument("--eos-user", required=True)
    p.add_argument("--eos-pwd", required=True)

    args = p.parse_args()

    assert os.path.exists(args.limitfile)
    assert os.path.exists(args.results)

    exit, string = analysis(args.limitfile, args.results, False, True)
    print(string)
    print()

    with open(os.path.join(os.path.dirname(__file__), "config.yml"), "r") as f:
        config = yaml.load(f) 

    gl = Gitlab(config["gitlab"]["gitlab_url"], private_token=args.gitlab_token)
    gl.auth()

    # we need to figure out if this ref belongs to an MR
    project = gl.projects.get(config["gitlab"]["gitlab_project"])
    # print(project)

    target_mr = None
    for mr in project.mergerequests.list(as_list=False, state="opened"):
        if mr.sha == args.ref:
            target_mr = mr
            break

    target_branch = None
    for branch in project.branches.list(as_list=False):
        if branch.commit["id"] == args.ref:
            target_branch = branch
            break

    sshfs = SSHFS(host=config["www"]["ssh_host"], user=args.eos_user, passwd=args.eos_pwd)
    destdir = sshfs.opendir(os.path.join(config["www"]["eos_www_root"], 
                                         config["www"]["static_analysis_dir"]))
    # destdir = OSFS("/tmp/eos_www_acts")

    if target_mr is not None:
        # this commit is the tip of an MR that is active
        prefix = os.path.join(project.name, "mr", str(target_mr.iid), "clang-tidy")
        publish(OSFS(args.report_dir), 
                destdir,
                prefix=prefix)

        public_url = urljoin(config["www"]["public_root"], config["www"]["static_analysis_dir"]+"/"+prefix)
        # comment and publish 
        note_body = "For commit: %s\n\n%s" % (args.ref, string)
        note_body += "\n\nAnalysis results at: [{0}]({0})".format(public_url)

        print("Commenting on", project.name, "!%d"%target_mr.iid)
        target_mr.notes.create({"body": note_body})


    if target_branch is not None:
        prefix = os.path.join(project.name, "branch", target_branch.name, "clang-tidy")

        publish(OSFS(args.report_dir),
                destdir,
                prefix=prefix)

    sys.exit(exit)

if "__main__" == __name__:
    main()

