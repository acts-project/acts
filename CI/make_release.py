#!/usr/bin/env python3

import re
import click
from sh import git
from pathlib import Path
import gitlab as gitlab_api
import sys
from datetime import datetime, date
from dateutil.relativedelta import relativedelta, FR
import jinja2
import humanize

from util import Spinner

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
        "--access-token", "-t", "gitlab", required=True, callback=gitlab_instance
    )(f)
    return f


def get_milestones(project, **kwargs):
    milestones = project.milestones.list(**kwargs)
    ms_dict = {}
    for ms in milestones:
        ms_dict[ms.title] = ms
    return ms_dict


def find_milestone(version, milestones):
    vstr = "{:d}.{:>0d}.{:>02d}".format(*version)
    ms_titles = (vstr, "v" + vstr)

    milestone = None
    for ms in milestones:
        #  print(ms.title, ms_title)
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
        if click.confirm(
            f"Do you want me to set due date of %{current_milestone.title} to {date.today()}? (is {current_milestone.due_date})"
        ):
            current_milestone.due_date = date.today().strftime(dfmt)
            current_milestone.save()

    dt = date.today()
    delta = relativedelta(weekday=FR(1))
    next_due = dt + delta
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
Bugfixes should be targeted at [`release/v0.12.X`](https://gitlab.cern.ch/acts/acts-core/tree/{{ release_branch  }}).

We will tag the next release `{{nv}}` on {{humanize.naturaldate(next_due)}} from [`%{{nm.title}}`](https://gitlab.cern.ch/acts/acts-core/-/milestones/{{nm.iid}}). 
This release can be cancelled if no merges happen before that date.
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


if "__main__" == __name__:
    main()
