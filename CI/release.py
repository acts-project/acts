#!/usr/bin/env python3

from pathlib import Path

import click
import github
from github import Github
import yaml
import re
from sh import git
from rich import print

from util import Spinner

default_branch_name = "master"

def get_current_branch():
  return git("rev-parse", "--abbrev-ref", "HEAD").strip()

def split_version(version):
    version_ex = re.compile(r"^v?(\d+)\.(\d{1,2})\.(\d{1,2})$")
    m = version_ex.match(version)
    assert m is not None, f"Version {version} is not in valid format"
    return tuple((int(m.group(i)) for i in range(1, 4)))

def format_version(version):
    return "v{:d}.{:>2d}.{:>02d}".format(*version)

def check_branch_exists(branch):
    with Spinner(f"Checking for {branch} branch"):
      all_branches = [l.strip() for l in git.branch(all=True, _tty_out=False).strip().split("\n")]
      for b in all_branches:
        if b.endswith(branch):
          return True
    return False

@click.group()
@click.option("--token", "-T", required=True, envvar="GITHUB_TOKEN")
@click.option("--repository", "-R", required=True, envvar="GITHUB_REPOSITORY")
@click.option("--retry", default=3)
@click.pass_context
def main(ctx, token, repository, retry):
    gh = Github(token, retry=retry)
    repo = gh.get_repo(repository)

    ctx.obj = gh, repo


def confirm(*args, yes=False, **kwargs):
  if yes == True:
    return True
  return click.confirm(*args, **kwargs)

@main.command()
@click.argument("tag_name")
@click.option("--remote", default="origin")
@click.option("--yes", "-y", is_flag=True, default=False)
@click.pass_obj
def tag(obj, tag_name, remote, yes):
  current_branch = get_current_branch()
  remote_url = git.remote("get-url", remote).strip()

  gh, repo = obj


  tag = split_version(tag_name)
  tag_name = format_version(tag)
  major, minor, fix = tag

  with Spinner(f"Checking for milestone for tag {tag_name}"):
      tag_milestone = None
      for ms in repo.get_milestones(state="all"):
          if ms.title == tag_name:
              tag_milestone = ms
              break
      assert tag_milestone is not None, "Did not find milestone for tag"

  release_branch_name = f"release/v{major}.{minor:>02}.X"

  with Spinner("Refreshing branches"):
    git.fetch(all=True, prune=True)

  if fix == 0:
    # new minor release
    with Spinner(f"Checking out and updating {default_branch_name}"):
      git.checkout(default_branch_name)
      git.pull()


    assert not check_branch_exists(release_branch_name), "For new minor: release branch CANNOT exist yet"

    with Spinner(f"Creating {release_branch_name}"):
      git.checkout("-b", release_branch_name)
  else:
    assert check_branch_exists(release_branch_name), "For new fix: release brunch MUST exist"

    with Spinner(f"Checking out {release_branch_name}"):
      git.checkout(release_branch_name)

  # we are not on release branch

  version_file = Path("version_number")
  assert version_file.exists(), "Version number file not found"

  current_version_string = version_file.read_text()
  print(f"Current version: [bold]{current_version_string}[/bold]")

  if fix == 0:
    assert current_version_string == "9.9.9", "Unexpected current version string found"
  else:
    assert current_version_string != f"{major}.{minor}.{fix-1}", "Unexpected current version string found"

  version_string = f"{major}.{minor}.{fix}"
  with Spinner(f"Bumping version number in '{version_file}' to '{version_string}'"):
    with version_file.open("w") as fh:
      fh.write(version_string)

  with Spinner("Comitting"):
    git.add(version_file)
    git.commit(m=f"Bump version number to {version_string}")

  with Spinner(f"Creating tag {tag_name}"):
    git.tag(tag_name)

  print(f"I will now: push tag [bold green]{tag_name}[/bold green] and branch [bold green]{release_branch_name}[/bold green] to [bold]{remote_url}[/bold]")
  if not confirm("Continue?", yes=yes):
    raise SystemExit("Aborting")

  with Spinner(f"Pushing branch {release_branch_name}"):
    git.push("-u", remote, release_branch_name)

  with Spinner(f"Pushing tag {tag_name}"):
    git.push(remote, tag_name)


@main.command()
@click.argument("tag_name")
@click.option("--draft/--publish", default=True)
@click.option("--yes", "-y", is_flag=True, default=False)
@click.pass_obj
def notes(obj, tag_name, draft, yes):
    gh, repo = obj

    with (Path(__file__).parent.parent / ".labels.yml").open("r") as fh:
        labels = yaml.safe_load(fh)["labels"]

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
              "", milestone=tag_milestone.title, repo=repo.full_name, type="pr", **{"is": "merged"}
            )
        )

    assert not any(
        [pr.state == "open" for pr in prs]
    ), "PRs assigned to milestone that are still open!"

    click.echo("Have " + click.style(str(len(prs)), bold=True) + " PRs, all closed.")

    body = ""

    groups = {l: [] for l in sorted(labels)}
    groups["Uncategorized"] = []

    for pr in prs:
        pr_labels = [l.name for l in pr.labels]

        assigned = False
        for label in labels:
            if label in pr_labels:
                groups[label].append(pr)
                assigned = True
                break
        if not assigned:
          groups["Uncategorized"].append(pr)

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
        if confirm(f"Update release {release.title}?", yes=yes):
            with Spinner(f"Updating release {release.title}"):
                release.update_release(name=release.title, message=body)
            click.echo(
                "Updated release is at {}".format(
                    click.style(release.html_url, bold=True)
                )
            )

    else:
        # new release
        if confirm(f"Create release for tag {tag.name} (draft: {draft})?", yes=yes):
            with Spinner(f"Creating release {tag.name}"):
                release = repo.create_git_release(
                    tag=tag.name, name=tag.name, message=body, draft=draft
                )
            click.echo(
                "Created release is at {}".format(
                    click.style(release.html_url, bold=True)
                )
            )
        else:
          print("Not creating a release")

    if tag_milestone.state == "open":
        if confirm(f"Do you want me to close milestone {tag_milestone.title}?", yes=yes):
            with Spinner(f"Closing milestone {tag_milestone.title}"):
                tag_milestone.edit(title=tag_milestone.title, state="closed")
        else:
          print("Not closing milestone")


main()
