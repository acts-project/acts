#!/usr/bin/env python
from __future__ import print_function

import cern_sso

import os
from urllib import quote_plus
import requests
import json
import re
import argparse

p = argparse.ArgumentParser()

token_default = os.environ["ATSJENKINS_ACCESS_TOKEN"] if "ATSJENKINS_ACCESS_TOKEN" in os.environ else ""
p.add_argument("--access-token", help="Gitlab access token to update the releases", default=token_default)

args = p.parse_args()

BASE_URL = "https://its.cern.ch/jira"
GITLAB_BASE_URI = "https://gitlab.cern.ch/api/v4/"
GITLAB_TOKEN = args.access_token

cookies = cern_sso.krb_sign_on(BASE_URL)

class JIRAException(Exception):
    def __init__(self, messages, *args, **kwargs):
        self.messages = messages
        # super().__init__(*args, **kwargs)

def jql(q):
    return requests.get(BASE_URL + "/rest/api/2/search?jql={}&maxResults=100".format(q), cookies=cookies)

def get_version_issues(version):
    res = jql("project= ACTS AND fixVersion = {} AND status = Closed".format(version)).json()
    try:
        assert res["maxResults"] > res["total"]
        assert res["startAt"] == 0

        return res["issues"]
    except:
        raise JIRAException(res["errorMessages"])


def get_tags(project):
    slug = quote_plus(project, safe="")
    url = GITLAB_BASE_URI+"projects/"+slug+"/repository/tags"
    res = requests.get(url, params={})
    assert res.ok, "Failure getting tags for {} (HTTP {})".format(project, res.status_code)
    data = res.json()
    return data

def update_tag(project, tag, desc, create=False):
    slug = quote_plus(project, safe="")
    url = GITLAB_BASE_URI + "projects/"+slug+"/repository/tags/" + tag + "/release"
    # print(url)

    if create:
        func = requests.post
    else:
        func = requests.put

    res = func(url, json={"tag_name": tag, "description": desc}, headers={"Private-Token": GITLAB_TOKEN})
    
    assert res.ok, "Failure updating tag {} for {} (HTTP {})".format(tag, project, res.status_code)

def make_release_notes(version):

    print("Making release notes for version", version)

    issues = get_version_issues(version)

    issues_by_type = {}

    for issue in issues:
        # print(issue)
        issue_type = issue["fields"]["issuetype"]["name"]
        summary = issue["fields"]["summary"]
        key = issue["key"]
        url = "https://its.cern.ch/jira/browse/{}".format(key)

        # print(summary, "\n", key, url, issue_type)

        if not issue_type in issues_by_type:
            issues_by_type[issue_type] = []

        issues_by_type[issue_type].append((key, url, summary))

    print("Found", len(issues), "issues with fix version ", version)

    markdown = ""
    markdown += "# Release v{}\n\n".format(version)

    for issue_type, issues in issues_by_type.iteritems():
            

        markdown += "## {}\n".format(issue_type)

        for key, url, summary in sorted(issues, key=lambda i: i[0]):
            markdown += " - [[{key}] {summary}]({url})\n".format(key=key, url=url, summary=summary)

        markdown += "\n"*2


    print("Release notes assembled")
    return markdown

tags = get_tags("acts/acts-core")

print("Found", len(tags), "tags")

for tag in tags:
    try:
        name = tag["name"]
        version = re.match(r"v(\d\.\d\d\.\d\d)", name).group(1)
        has_release = tag["release"] is not None

        print(name)
        relnotes = make_release_notes(version)

        if not has_release:
            print("Creating release for tag", name)
        else:
            print("Updating release notes for tag", name)

        update_tag("acts/acts-core", name, relnotes, not has_release)



    except JIRAException as e:
        print("Skipping", tag["name"], ", reason:")
        for m in e.messages:
            print(" ->", m)

    print()


print("Release note synchronization complete")
