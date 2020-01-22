#!/usr/bin/env python3

import re
import click
from sh import git
from pathlib import Path
import gitlab as gitlab_api
import sys
from datetime import datetime, date
from dateutil.relativedelta import relativedelta, FR
import dateutil.parser
import jinja2
import json
import humanize
import click_config_file
import requests
import os
from pprint import pprint
import tempfile

from util import Spinner
from release_notes import (
    collect_milestone,
    make_release_notes,
    get_label_groups,
)

version_ex = re.compile(r"^v?(\d+)\.(\d{1,2})\.(\d{1,2})$")


def split_version(version):
    m = version_ex.match(version)
    assert m is not None, f"Version {version} is not in valid format"
    return [int(m.group(i)) for i in range(1, 4)]


def format_version(version):
    return "v{:d}.{:>2d}.{:>02d}".format(*version)


def gitlab_instance(ctx, param, token):
    gl = gitlab_api.Gitlab("https://gitlab.cern.ch", private_token=token)
    gl.auth()
    return gl


def gitlab_option(f):
    f = click.option(
        "--gitlab-token", "-t", "gitlab", required=True, callback=gitlab_instance
    )(f)
    return f


def get_milestones(project, **kwargs):
    milestones = project.milestones.list(**kwargs)
    ms_dict = {}
    for ms in milestones:
        try:
          ms_dict[tuple(split_version(ms.title))] = ms
        except:
          pass
    return ms_dict


def find_milestone(version, milestones):
    vstr = "{:d}.{:>02d}.{:>02d}".format(*version)
    ms_titles = (vstr, "v" + vstr)

    milestone = None
    for ms in milestones:
        #  print(ms.title, ms_titles)
        if ms.title in ms_titles:
            milestone = ms
            break

    return milestone


def get_branches():
    branches = (
        git("for-each-ref", "refs/heads", format="%(refname:short)").strip().split("\n")
    )
    return branches


def current_branch():
    branch = git.branch(show_current=True).strip()
    print(branch)
    return branch


@click.group()
def main():
    pass


@main.command()
@gitlab_option
@click.option("--dry-run", is_flag=True)
@click.argument("version")
def minor(version, dry_run, gitlab):
    project = gitlab.projects.get("acts/acts-core")
    version = split_version(version)
    milestone = find_milestone(version, project.milestones.list(state="active"))
    assert (
        milestone is not None
    ), f"Didn't find milestone for {version}. Is it closed already?"

    branches = get_branches()

    release_branch = "release/v{:d}.{:>02d}.X".format(*version)
    source_branch = "master"  # always master for minor version
    version_file = Path() / "version_number"
    tag_name = format_version(version)

    print(
        "Will make new release with version %s from milestone %s and branch %s"
        % (format_version(version), milestone.title, source_branch)
    )

    if click.confirm("Do you want to run local preparation?"):
        if source_branch not in branches:
            print("Source branch", source_branch, "not found.")
            sys.exit(1)

        with Spinner(text=f"Checkout and update source branch {source_branch}"):
            if not dry_run:
                git.checkout(source_branch)
                assert current_branch() == source_branch
                git.pull()

        if release_branch in branches and not dry_run:
            print("Release branch", release_branch, "exists. I'm bailing")
            sys.exit(1)

        with Spinner(text=f"Creating {release_branch} from {source_branch}"):
            if not dry_run:
                git.checkout("-b", release_branch)

        with Spinner(text=f"Bumping version to {format_version(version)}"):
            if not dry_run:
                assert current_branch() == release_branch
                with version_file.open("w") as fh:
                    fh.write(".".join(map(str, version)))

        with Spinner(
            text=f"Committing bumped version on release branch {release_branch}"
        ):
            if not dry_run:
                git.add(str(version_file))
                git.commit(message="Bump version to %s" % ".".join(map(str, version)))

        with Spinner(text=f"Creating local tag {tag_name} on {release_branch}"):
            if not dry_run:
                git.tag(tag_name)
        print(f"You might want to run 'git push REMOTE {tag_name}'")

    if click.confirm(f"Do you want me to try to push {release_branch}?"):
        with Spinner(text=f"Pushing {release_branch}"):
            if not dry_run:
                git.push()

    if click.confirm(f"Do you want me to try to push tag {tag_name}?"):
        with Spinner(text=f"Pushing {tag_name}"):
            if not dry_run:
                git.push("REMOTE", tag_name)

    if click.confirm(f"Do you want me to close %{milestone.title}?"):
        with Spinner(text=f"Closing milestone %{milestone.title}"):
            if not dry_run:
                milestone.state_event = "close"
                milestone.save()

    print("Done!")
    if dry_run:
        print("THIS WAS A DRY RUN!")

@main.command()
@gitlab_option
@click.option("--dry-run", is_flag=True)
@click.argument("version")
def patch(version, dry_run, gitlab):
    project = gitlab.projects.get("acts/acts-core")

    version = split_version(version)
    milestone = find_milestone(version, project.milestones.list(state="active"))
    assert (
        milestone is not None
    ), f"Didn't find milestone for {version}. Is it closed already?"

    branches = get_branches()

    release_branch = "release/v{:d}.{:>02d}.X".format(*version)
    version_file = Path() / "version_number"
    tag_name = format_version(version)

    if release_branch not in branches:
      print("Release branch", release_branch, "does not exist. I'm bailing")

    print(
        "Will make new patch version tag %s from milestone %s on branch %s"
        % (format_version(version), milestone.title, release_branch)
    )

    if click.confirm("Do you want to run local preparation?"):

      with Spinner(text=f"Checkout and update release branch {release_branch}"):
          if not dry_run:
              git.checkout(release_branch)
              assert current_branch() == release_branch
              git.pull()

      with Spinner(text=f"Bumping version to {format_version(version)}"):
          if not dry_run:
              assert current_branch() == release_branch
              with version_file.open("w") as fh:
                  fh.write(".".join(map(str, version)))

      with Spinner(
          text=f"Committing bumped version on release branch {release_branch}"
      ):
          if not dry_run:
              git.add(str(version_file))
              git.commit(message="Bump version to %s" % ".".join(map(str, version)))

      with Spinner(text=f"Creating local tag {tag_name} on {release_branch}"):
          if not dry_run:
              git.tag(tag_name)
      print(f"You might want to run 'git push REMOTE {tag_name}'")

    if click.confirm(f"Do you want me to try to push {release_branch}?"):
        with Spinner(text=f"Pushing {release_branch}"):
            if not dry_run:
                git.push()

    if click.confirm(f"Do you want me to close %{milestone.title}?"):
        with Spinner(text=f"Closing milestone %{milestone.title}"):
            if not dry_run:
                milestone.state_event = "close"
                milestone.save()



@main.command()
@gitlab_option
@click.argument("version")
def message(version, gitlab):
    dfmt = "%Y-%m-%d"

    current_version = split_version(version)
    next_version = current_version[:]
    next_version[1] += 1

    print(current_version, next_version, file=sys.stderr)

    project = gitlab.projects.get("acts/acts-core")

    milestones = project.milestones.list()

    current_milestone = find_milestone(current_version, milestones)
    assert current_milestone is not None

    next_milestone = find_milestone(next_version, milestones)

    if next_milestone is None:
        print("Milestone for", format_version(next_version), "does not exist")
        if click.confirm("Want me to create it?"):
            title = click.prompt("What title?", format_version(next_version))
            next_milestone = project.milestones.create(title=title)
        else:
            sys.exit(1)

    if current_milestone.due_date != date.today().strftime(dfmt):
        if sys.stdout.isatty():
          if click.confirm(
              f"Do you want me to set due date of %{current_milestone.title} to {date.today()}? (is {current_milestone.due_date})"
          ):
              current_milestone.due_date = date.today().strftime(dfmt)
              current_milestone.save()

    if next_milestone.due_date is None:
      dt = date.today()
      delta = relativedelta(weekday=FR(1))
      next_due = dt + delta
    else:
      next_due = datetime.strptime(next_milestone.due_date, dfmt)
    if sys.stdout.isatty():
        next_due = datetime.strptime(
            click.prompt(
                f"Due date for milestone %{next_milestone.title}",
                next_due.strftime(dfmt),
            ),
            dfmt,
        )

    start_date = datetime.strptime(next_milestone.start_date, dfmt) or date.today()
    start_date_str = start_date.strftime(dfmt)
    next_due_str = next_due.strftime(dfmt)

    if (
        next_milestone.start_date != start_date_str
        or next_milestone.due_date != next_due_str
    ):
        if click.confirm(f"Update milestone %{next_milestone.title}?"):
            with Spinner(text=f"Updating milestone %{next_milestone.title}"):
                next_milestone.start_date = start_date_str
                next_milestone.due_date = next_due_str
                next_milestone.save()

    release_branch = "release/v{:d}.{:>02d}.X".format(*current_version)

    tpl = jinja2.Template(
        """
I've just tagged [`{{cv}}`](https://gitlab.cern.ch/acts/acts-core/-/tags/{{cv}}) from milestone [`%{{cm.title}}`](https://gitlab.cern.ch/acts/acts-core/-/milestones/{{cm.iid}}). 
Bugfixes should be targeted at [`{{ release_branch }}`](https://gitlab.cern.ch/acts/acts-core/tree/{{ release_branch  }}).

We will tag the next release `{{nv}}` on {{humanize.naturaldate(next_due)}} from [`%{{nm.title}}`](https://gitlab.cern.ch/acts/acts-core/-/milestones/{{nm.iid}}). 
This release can be cancelled if a sufficient number of merges does not  happen before that date.
""".strip()
    )

    tpl.globals["humanize"] = humanize

    text = tpl.render(
        next_due=datetime.strptime(next_milestone.due_date, dfmt).date(),
        release_branch=release_branch,
        cm=current_milestone,
        cv=format_version(current_version),
        nm=next_milestone,
        nv=format_version(next_version),
    )

    print(text)

@main.command()
@gitlab_option
@click.argument("start", callback=lambda c, p, s: split_version(s))
@click.argument("end", callback=lambda c, p, s: split_version(s))
def relnotes(start, end, gitlab):
    start = tuple(start)
    end = tuple(end)
    print(start, end, file=sys.stderr)
    project = gitlab.projects.get("acts/acts-core")

    all_milestones = get_milestones(project)
    milestones = []
    for ms in all_milestones.values():
      try:
        ver = split_version(ms.title)
        milestones.append(ms)
      except:
        pass

    sorted_milestones = list(sorted(all_milestones.keys()))

    start_ms = all_milestones[start]
    end_ms = all_milestones[end]

    ms_range = (sorted_milestones.index(start), sorted_milestones.index(end))

    md = ""

    for mst in sorted_milestones[ms_range[0]+1:ms_range[1]+1]:
      ms = all_milestones[mst]
      print(ms.title, file=sys.stderr)
      

      mrs_grouped, issues_grouped = collect_milestone(ms)
      with Spinner(text="Assembling release notes", stream=sys.stderr):
          md += f"## {format_version(mst)}\n\n"
          md += make_release_notes(ms, mrs_grouped, issues_grouped, badges=False, links=False)


    print(md)

class Zenodo:
  def __init__(self, token, base_url = "https://zenodo.org/api/"):
    self.token = token
    self.base_url = base_url

  def get(self, url, params = {}, **kwargs):
    _params = {"access_token": self.token}
    _params.update(params)
    r =  requests.get(os.path.join(self.base_url, url), params=_params, **kwargs)
    return r.json()

  def post(self, url, params = {}, headers = {}, **kwargs):
    _headers = {"Content-Type": "application/json"}
    _headers.update(headers)
    _params = {"access_token": self.token}
    _params.update(params)
    r = requests.post(os.path.join(self.base_url, url), params=_params, 
                      headers=_headers, 
                      **kwargs)
    assert r.status_code == 201, r.json()
    return r

  def put(self, url, data, params = {}, headers = {}, **kwargs):
    _headers = {"Content-Type": "application/json"}
    _headers.update(headers)
    #  _params = {"access_token": self.token}
    #  _params.update(params)
    _url = os.path.join(self.base_url, url)+f"?access_token={self.token}"
    print(_url)
    r = requests.put(_url,
                      data=json.dumps(data),
                      headers=_headers, 
                      **kwargs)
    assert r.status_code == 200, f"Status {r.status_code}, {r.json()}"
    return r

  def upload(self, deposition, name, fh):
    _params = {"access_token": self.token}
    data={"name": name}
    files = {"file": fh}
    r = requests.post(os.path.join(self.base_url, f"deposit/depositions/{deposition}/files"), 
                      params=_params, data=data, files=files)
    assert r.status_code == 201, r.status_code
    return r

  def delete(self, url, params = {}, **kwargs):
    _params = {"access_token": self.token}
    return requests.delete(os.path.join(self.base_url, url), params=_params)

@main.command()
@gitlab_option
@click.argument("version")
@click.option("--zenodo-token", "-z", required=True)
@click.option("--deposition", "-d", required=True)
@click_config_file.configuration_option()
def zenodo(version, gitlab, zenodo_token, deposition):
  version = split_version(version)
  print(version, gitlab, zenodo_token)
  zenodo = Zenodo(zenodo_token)

  create_res = zenodo.post(f"deposit/depositions/{deposition}/actions/newversion")
  print(create_res)
  create_res = create_res.json()

  draft_id = create_res["links"]["latest_draft"].split("/")[-1]
  pprint(create_res)

  print("Created new version with id", draft_id)

  print("Delete all files for draft")
  draft = zenodo.get(f"deposit/depositions/{draft_id}")
  pprint(draft)

  for file in draft["files"]:
    file_id = file["id"]
    r = zenodo.delete(f"deposit/depositions/{draft_id}/files/{file_id}")
    assert r.status_code == 204

  creator_file = os.path.join(os.path.dirname(__file__), "../AUTHORS.md")
  with open(creator_file) as fh:
    md = fh.read().strip().split("\n")
  md = [l.strip() for l in md if not l.strip().startswith("#") and not l.strip() == ""]

  creators = []
  for line in md:
    assert line.startswith("- ")
    line = line[2:]
    split = line.split(",", 1)
    creator = {"name": split[0].strip()}

    if len(split) == 2:
      creator["affiliation"] = split[1].strip()

    creators.append(creator)
	
  project = gitlab.projects.get("acts/acts-core")
  milestones = project.milestones.list()
  milestone = find_milestone(version, milestones)
  mrs_grouped, issues_grouped = collect_milestone(milestone)

  assert milestone.state == "closed"

  tag = project.tags.get(format_version(version))
  print(tag)
  tag_date = dateutil.parser.parse(tag.commit["created_at"]).date().strftime("%Y-%m-%d")

  description = f'Milestone: <a href="{milestone.web_url}">%{milestone.title}</a> <br/> Merge requested accepted for this version: \n <ul>\n'

  for mr in sum(mrs_grouped.values(), []):
    description += f'<li><a href="{mr.web_url}">!{mr.iid} - {mr.title}</a></li>\n'

  description += "</ul>"

  data = {"metadata": {
    "title": f"Acts Project: {format_version(version)}",
    "upload_type": "software",
    "description": description,
    "creators": creators,
    "version": format_version(version),
    "publication_date": tag_date,
    "license": "MPL-2.0",
  }}
  zenodo.put(f"deposit/depositions/{draft_id}", data).json()


  with tempfile.TemporaryFile() as fh:
    r = requests.get(f"https://gitlab.cern.ch/acts/acts-core/-/archive/{format_version(version)}/acts-core-{format_version(version)}.zip", stream=True)
    r.raw.decode_content = True
    fh.write(r.raw.read())
    fh.seek(0)
    name = f"acts-core-{format_version(version)}.zip"
    zenodo.upload(draft_id, name, fh)

if "__main__" == __name__:
    main()


