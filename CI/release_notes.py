from util import Spinner
import os
import re
import urllib.parse
import requests


def parse_version(name):
    return re.match(r"v(\d\.\d\d\.\d\d)", name).group(1)


def get_label_groups():
    return os.getenv("RELEASE_NOTES_LABEL_GROUPS", "New Feature;Bug;Improvement").split(
        ";"
    )


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


def collect_milestone(milestone):
    label_groups = get_label_groups()
    mrs = []
    with Spinner(text=f"Loading merge requests associated with %{milestone.iid}"):
        for mr in milestone.merge_requests():
            if mr.state == "merged":
                mrs.append(mr)

    # need to get issues from merged MRs
    with Spinner(text=f"Collecting issues from {len(mrs)} merged MRs"):
        issue_ids = []
        issues = []
        for mr in mrs:
            for issue in mr.closes_issues():
                if issue.id not in issue_ids:
                    issue_ids.append(issue.id)
                    issues.append(issue)

    issues_grouped = group_items(label_groups, issues)
    mrs_grouped = group_items(label_groups, mrs)

    return mrs_grouped, issues_grouped


def make_release_notes(milestone, mrs_grouped, issues_grouped):
    label_groups = get_label_groups()
    # make the Markdown
    md = ""
    # md = f"# Release {tag.name}\n"
    # md += f"Milestone: [%{milestone.title}]({milestone.web_url})\n"

    ms_badge = "https://badgen.net/badge/milestone/%s/green" % urllib.parse.quote(
        f"%{milestone.title}"
    )

    md += f"[![]({ms_badge})]({milestone.web_url})\n"

    badge = f"![](https://gitlab.cern.ch/acts/acts-core/badges/v{milestone.title}/coverage.svg)"
    cov_url = f"https://acts.web.cern.ch/ACTS/coverage/v{milestone.title}/"

    r = requests.get(cov_url)
    if r.status_code == 200:
        md += f"[{badge}]({cov_url})\n\n"
    else:
        md += f"[{badge}](#)\n\n"

    nmrs = sum([len(v) for v in mrs_grouped.values()])
    nissues = sum([len(v) for v in issues_grouped.values()])

    if nmrs > 0:
        md += f"### {nmrs} Merge Requests in this release:\n"
        for g in label_groups + ["Uncategorized"]:
            if len(mrs_grouped[g]) == 0:
                continue
            md += f"#### {g}\n"
            for mr in mrs_grouped[g]:
                md += f"- [!{mr.iid} - {mr.title}]({mr.web_url})\n"
        md += "\n"

    if nissues > 0:
        md += f"### {nissues} issues addressed in this release:\n"
        for g in label_groups + ["Uncategorized"]:
            if len(issues_grouped[g]) == 0:
                continue
            md += f"#### {g}\n"
            for issue in issues_grouped[g]:
                md += f"- [#{issue.iid} - {issue.title}]({issue.web_url})\n"
    return md
