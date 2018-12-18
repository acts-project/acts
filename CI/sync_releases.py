#!/usr/bin/env python3

import cern_sso
import os
from urllib.parse import quote as quote_plus
import requests
import json
import re
import argparse
import gitlab
from concurrent.futures import ThreadPoolExecutor, as_completed
import dateutil.parser
import contextlib
import sys
import math


# these are optional, but nice
try:
    from tqdm import tqdm
except:
    tqdm = None
try:
    from halo import Halo
except:
    Halo = None

class JIRAException(Exception):
    def __init__(self, messages, *args, **kwargs):
        self.messages = messages

class JIRA:
    def __init__(self, cookies, url):
        self.cookies = cookies
        self.url = url

    def jql(self, q):
        return requests.get(self.url + "/rest/api/2/search?jql={}&maxResults=500".format(q), cookies=self.cookies)

    def get_version_issues(self, version):
        res = self.jql("project= ACTS AND fixVersion = {} AND status = Closed".format(version)).json()
        try:
            assert res["maxResults"] > res["total"]
            assert res["startAt"] == 0

            return res["issues"]
        except:
            raise JIRAException(res["errorMessages"])


@contextlib.contextmanager
def spinner(text, *args, **kwargs):
    if sys.stdout.isatty() and Halo is not None:
        with Halo(text, *args, **kwargs):
            yield
    else:
        sys.stdout.write(text+"\n")
        yield

if sys.stdout.isatty() and tqdm is not None:
    Progress = tqdm
    prog_iter = tqdm
    prog_write = tqdm.write
else:
    class Progress():
        def __init__(self, total, desc, *args, **kwargs):
            self.total = total
            self.current = 0
            sys.stdout.write(desc+"\n")
        def update(self, n=1, *args, **kwargs):
            self.current += n
            perc = self.current / self.total * 100
            inc = math.ceil(self.total / 10)
            if self.current % inc == 0:
                sys.stdout.write("%.2f"%perc + "%\n")
        def close(self):
            pass

    def prog_iter(values, *args, **kwargs):
        p = Progress(len(values), *args, **kwargs)
        for value in values:
            yield value
            p.update()
        p.close()

    def prog_write(*args, **kwargs):
        sys.stdout.write(*args, **kwargs)
        sys.stdout.write("\n")

def mtmap(tp, func, values, desc=None):
    prog = Progress(total=len(values), leave=False, desc=desc)
    futures = []
    for v in values:
        futures.append(tp.submit(func, v))
    
    for _ in as_completed(futures):
        prog.update()
    prog.close()
    return [f.result() for f in futures]


def make_release_notes(version, issues, mrs):

    issues_by_type = {}

    for issue in issues:
        issue_type = issue["fields"]["issuetype"]["name"]
        summary = issue["fields"]["summary"]
        key = issue["key"]
        url = "https://its.cern.ch/jira/browse/{}".format(key)

        if not issue_type in issues_by_type:
            issues_by_type[issue_type] = []

        issues_by_type[issue_type].append((key, url, summary))


    markdown = ""
    markdown += "# Release v{}\n\n".format(version)

    markdown += "\n"*2 + "Merge requests for this release:\n\n"
    for mr in mrs:
        markdown += " - [!%d - %s](%s)" % (mr.iid, mr.title, mr.web_url) + "\n"
    markdown += "\n"*2

    for issue_type, issues in issues_by_type.items():

        markdown += "## {}\n".format(issue_type)

        for key, url, summary in sorted(issues, key=lambda i: i[0]):
            markdown += " - [[{key}] {summary}]({url})\n".format(key=key, url=url, summary=summary)

        markdown += "\n"*2

    return markdown

def parse_version(tag):
    return re.match(r"v(\d\.\d\d\.\d\d)", tag.name).group(1)

def main():
    p = argparse.ArgumentParser()

    p.add_argument("--access-token",
                   help="Gitlab access token to update the releases",
                   default=os.getenv("ATSJENKINS_ACCESS_TOKEN", None))
    p.add_argument("--dry-run", "-s", action="store_true")

    args = p.parse_args()

    jira_url = "https://its.cern.ch/jira"
    cookies = cern_sso.krb_sign_on(jira_url)
    jira = JIRA(cookies=cookies, url=jira_url)

    gl = gitlab.Gitlab("https://gitlab.cern.ch", private_token=args.access_token)
    if not args.dry_run:
        assert args.access_token is not None
        gl.auth()
    project = gl.projects.get("acts/acts-core")

    with spinner(text="Loading tags"):
        tags = project.tags.list(all=True)

    with spinner(text="Loading merge requests"):
        mrlist = project.mergerequests.list(state="merged", target_branch="master", all=True)

    with ThreadPoolExecutor(max_workers=15) as tp:
        mrs = mrlist

        for tag in tags:
            date = dateutil.parser.parse(tag.commit["created_at"])
            tag.created_at = date
        tags = list(sorted(tags, key=lambda t: t.created_at))

        def augment_with_commit(mr):
            commit = project.commits.get(mr.sha)
            date = dateutil.parser.parse(commit.created_at)
            mr.commit_date = date
            return mr

        mrs = mtmap(tp, augment_with_commit, mrs, desc="Loading MR commit info")

        def load_issues(tag):
            version = parse_version(tag)
            try:
                return tag, jira.get_version_issues(version)
            except JIRAException:
                return tag, []
        version_issues = dict(mtmap(tp, load_issues, tags, desc="Loading issues from JIRA"))


    tag_mrs = {}


    tag_mrs[tags[0]] = []
    for mr in mrs:
        if tags[0].created_at > mr.commit_date:
            tag_mrs[tags[0]].append(mr)

    for tag, tagn in zip(tags, tags[1:]):
        tag_mrs[tagn] = []
        for mr in mrs:
            if tag.created_at < mr.commit_date < tagn.created_at:
                tag_mrs[tagn].append(mr)

    print("Found", len(tags), "tags")


    for tag in prog_iter(tags, desc="Updating tag release notes"):
        name = tag.name
        version = parse_version(tag)
        has_release = tag.release is not None

        prog_write(name)
        relnotes = make_release_notes(version, version_issues[tag], tag_mrs[tag])

        if not has_release:
            prog_write("Creating release for tag %s"%name)
        else:
            prog_write("Updating release notes for tag %s"%name)

        if not args.dry_run:
            tag.set_release_description(relnotes)

            prog_write("Release notes for %s set" % name)

    print("Release note synchronization complete")


if "__main__" == __name__:
    main()
