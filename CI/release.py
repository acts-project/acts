#!/usr/bin/env python3

from pathlib import Path

import click
import github
from github import Github
import yaml

from util import Spinner


@click.command()
@click.option("--token", "-T", required=True, envvar="GITHUB_TOKEN")
@click.option("--repository", "-R", required=True, envvar="GITHUB_REPOSITORY")
@click.argument("tag_name")
@click.option("--retry", default=3)
@click.option("--draft/--publish", default=True)
def main(token, repository, tag_name, retry, draft):

    with (Path(__file__).parent.parent / ".labels.yml").open("r") as fh:
        labels = yaml.safe_load(fh)["labels"]

    gh = Github(token, retry=retry)
    repo = gh.get_repo(repository)

    with Spinner(f"Finding tag {tag_name}"):
        tag = None
        for t in repo.get_tags():
            if t.name == tag_name:
                tag = t
                break
        assert tag is not None, "Did not find tag"

    with Spinner(f"Loading milestone for tag {tag_name}"):
        tag_milestone = None
        for ms in repo.get_milestones(state="all"):
            if ms.title == tag_name:
                tag_milestone = ms
                break
        assert tag_milestone is not None, "Did not find milestone for tag"

    with Spinner(f"Getting PRs for milestone {tag_milestone.title}"):
        
        prs = list(
            gh.search_issues(
              "", milestone=tag_milestone.title, repo=repository, type="pr", **{"is": "merged"}
            )
        )

    assert not any(
        [pr.state == "open" for pr in prs]
    ), "PRs assigned to milestone that are still open!"

    click.echo("Have " + click.style(str(len(prs)), bold=True) + " PRs, all closed.")

    body = ""

    groups = {l: [] for l in sorted(labels)}

    for pr in prs:
        pr_labels = [l.name for l in pr.labels]

        for label in labels:
            if label in pr_labels:
                groups[label].append(pr)
                break

    for group, prs in groups.items():
        if len(prs) == 0:
            continue
        name = group
        if name.lower() == "bug":
          name = "Bug Fixes"
        body += f"#### {name}:\n\n"
        for pr in prs:
            body += f"- {pr.title} [#{pr.number}]({pr.html_url})\n"
        body += "\n"

    body = body.strip()

    width, _ = click.get_terminal_size()

    print()
    click.secho(
        "\n".join([l.ljust(width) for l in [""] + body.split("\n") + [""]]),
        fg="black",
        bg="white",
    )
    print()

    release = None
    with Spinner("Getting release"):
        try:
            release = repo.get_release(tag.name)

        except github.UnknownObjectException:
            pass

    if release is not None:
        # existing release, update

        click.echo(
            "Existing release {} is at {}".format(
                click.style(release.title, bold=True),
                click.style(release.html_url, bold=True),
            )
        )
        if click.confirm(f"Update release {release.title}?"):
            with Spinner(f"Updating release {release.title}"):
                release.update_release(name=release.title, message=body)
            click.echo(
                "Updated release is at {}".format(
                    click.style(release.html_url, bold=True)
                )
            )

    else:
        # new release
        if click.confirm(f"Create release for tag {tag.name} (draft: {draft})?"):
            with Spinner(f"Creating release {tag.name}"):
                release = repo.create_git_release(
                    tag=tag.name, name=tag.name, message=body, draft=draft
                )
            click.echo(
                "Created release is at {}".format(
                    click.style(release.html_url, bold=True)
                )
            )

    if tag_milestone.state == "open":
        if click.confirm(f"Do you want me to close milestone {tag_milestone.title}?"):
            with Spinner(f"Closing milestone {tag_milestone.title}"):
                tag_milestone.edit(title=tag_milestone.title, state="closed")


main()
