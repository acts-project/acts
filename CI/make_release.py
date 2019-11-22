#!/usr/bin/env python3

import argparse
import re
import click
from sh import git
from pathlib import Path

from util import def_arguments, Spinner, gitlab

version_ex = re.compile(r"^v?(\d+)\.(\d{1,2})\.(\d{1,2})$")


def split_version(version):
    m = version_ex.match(version)
    assert m is not None, f"Version {version} is not in valid format"
    return [int(m.group(i)) for i in range(1, 4)]


def format_version(version):
    return "v{:d}.{:>2d}.{:>02d}".format(*version)


def main():
    p = argparse.ArgumentParser()
    p = def_arguments(p, gl=True)

    p.add_argument("--dry-run", "-s", action="store_true")
    p.add_argument("--verbose", "-v", action="store_true")

    p.add_argument("version")

    args = p.parse_args()

    gl = gitlab(args)

    project = gl.projects.get("acts/acts-core")

    version = split_version(args.version)
    ms_title = "{:d}.{:>0d}.{:>02d}".format(*version)

    milestones = project.milestones.list(state="active")
    milestone = None
    for ms in milestones:
        print(ms.title, ms_title)
        if ms.title == ms_title:
            milestone = ms
            break

    assert (
        milestone is not None
    ), f"Didn't find milestone for {args.version}. Is it closed already?"

    print(
        "Will make new release with version %s from milestone %s"
        % (format_version(version), milestone.title)
    )

    if not args.dry_run:
        if not click.confirm("Are you sure?"):
            return

    with Spinner(text="Updating local clone"):
        if not args.dry_run:
            git.fetch(all=True)
            git.checkout("master")
            git.pull()

    release_branch = "release/v{:d}.{:>02d}.X".format(*version)
    branches = (
        git("for-each-ref", "refs/heads", format="%(refname:short)").strip().split("\n")
    )

    if release_branch not in branches:
        print(
            "Release branch",
            release_branch,
            "doesn't exist. I will not attempt to create it. Please do that manually",
        )

    with Spinner(text=f"Checkout release branch {release_branch}"):
        if not args.dry_run:
            git.checkout(release_branch)
            git.pull()

    with Spinner(text=f"Merging master into {release_branch}"):
        if not args.dry_run:
            git.merge("master")

    version_file = Path() / "version_number"
    with Spinner(text=f"Bumping version to {format_version(version)}"):
        if not args.dry_run:
            with version_file.open("w") as fh:
                fh.write(".".join(map(str, version)))

    with Spinner(text="Committing bumped version"):
        if not args.dry_run:
            git.add(str(version_file))
            git.commit(message="Bump version to %s" % ".".join(map(str, version)))

    tag_name = format_version(version)
    print(
        f"Steps completed. I will now push {release_branch}, create tag {tag_name} and close milestone %{milestone.title}"
    )

    if not args.dry_run:
        if not click.confirm("Continue?"):
            return

    with Spinner(text=f"Pushing {release_branch}"):
        if not args.dry_run:
            git.push()

    with Spinner(text=f"Creating tag {tag_name} on {release_branch}"):
        if not args.dry_run:
            project.tags.create({"tag_name": tag_name, "ref": release_branch})

    with Spinner(text=f"Closing milestone {milestone.title}"):
        if not args.dry_run:
            milestone.state_event = "close"
            milestone.save()

    print("Done!")
    if args.dry_run:
        print("THIS WAS A DRY RUN!")


if "__main__" == __name__:
    main()
