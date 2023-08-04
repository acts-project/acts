#!/usr/bin/env python3
import os
import asyncio
from typing import List, Optional, Tuple
import re
from pathlib import Path
import sys
import http
import json
import yaml
import datetime
import typer
import base64

import aiohttp
from gidgethub.aiohttp import GitHubAPI
from gidgethub import InvalidField
import gidgethub
from semantic_release.history import angular_parser, get_new_version
from semantic_release.errors import UnknownCommitMessageStyleError
from semantic_release.history.logs import LEVELS
from semantic_release.history.parser_helpers import ParsedCommit
import sh
from dotenv import load_dotenv
import functools

load_dotenv()

git = sh.git

RETRY_COUNT = 10
RETRY_INTERVAL = 0.5  # seconds


def get_repo():
    repo = os.environ.get("GITHUB_REPOSITORY", None)
    if repo is not None:
        return repo

    origin = git.remote("get-url", "origin")
    _, loc = origin.split(":", 1)
    repo, _ = loc.split(".", 1)
    return repo


def get_current_version():
    raw = git.describe().split("-")[0]
    m = re.match(r"v(\d+\.\d+\.\d+)", raw)
    return m.group(1)


class Commit:
    sha: str
    message: str
    author: str

    def __init__(self, sha: str, message: str, author: str):
        self.sha = sha
        self.message = self._normalize(message)
        self.author = author

    @staticmethod
    def _normalize(message):
        message = message.replace("\r", "\n")
        return message

    def __str__(self):
        message = self.message.split("\n")[0]
        return f"Commit(sha='{self.sha[:8]}', message='{message}')"


_default_parser = angular_parser


def evaluate_version_bump(
    commits: List[Commit], commit_parser=_default_parser
) -> Optional[str]:
    """
    Adapted from: https://github.com/relekang/python-semantic-release/blob/master/semantic_release/history/logs.py#L22
    """
    bump = None

    changes = []
    commit_count = 0

    for commit in commits:
        commit_count += 1
        try:
            message = commit_parser(commit.message)
            changes.append(message.bump)
        except UnknownCommitMessageStyleError:
            print("Unknown commit message style!")

    if changes:
        level = max(changes)
        if level in LEVELS:
            bump = LEVELS[level]
        else:
            print(f"Unknown bump level {level}")

    return bump


def generate_changelog(commits, commit_parser=_default_parser) -> dict:
    """
    Modified from: https://github.com/relekang/python-semantic-release/blob/48972fb761ed9b0fb376fa3ad7028d65ff407ee6/semantic_release/history/logs.py#L78
    """
    changes: dict = {"breaking": []}

    for commit in commits:
        try:
            message: ParsedCommit = commit_parser(commit.message)
            if message.type not in changes:
                changes[message.type] = list()

            capital_message = (
                message.descriptions[0][0].upper() + message.descriptions[0][1:]
            )
            changes[message.type].append((commit.sha, capital_message, commit.author))

            if message.breaking_descriptions:
                for paragraph in message.breaking_descriptions:
                    changes["breaking"].append((commit.sha, paragraph, commit.author))
            elif message.bump == 3:
                changes["breaking"].append(
                    (commit.sha, message.descriptions[0], commit.author)
                )

        except UnknownCommitMessageStyleError:
            print("Unknown commit message style!")

    return changes


def markdown_changelog(version: str, changelog: dict, header: bool = False) -> str:
    output = f"## v{version}\n" if header else ""

    for section, items in changelog.items():
        if len(items) == 0:
            continue
        output += "\n### {0}\n".format(section.capitalize())

        for sha, msg, author in items:
            output += "* {} ({}) (@{})\n".format(msg, sha, author)

    return output


def update_zenodo(zenodo_file: Path, repo: str, next_version):
    data = json.loads(zenodo_file.read_text())
    data["title"] = f"{repo}: v{next_version}"
    data["version"] = f"v{next_version}"
    zenodo_file.write_text(json.dumps(data, indent=2))


def update_citation(citation_file: Path, next_version):
    with citation_file.open() as fh:
        data = yaml.safe_load(fh)
    data["version"] = f"v{next_version}"
    data["date-released"] = datetime.date.today().strftime("%Y-%m-%d")
    with citation_file.open("w") as fh:
        yaml.dump(data, fh, indent=2)


def make_sync(fn):
    @functools.wraps(fn)
    def wrapped(*args, **kwargs):
        loop = asyncio.get_event_loop()
        loop.run_until_complete(fn(*args, **kwargs))

    return wrapped


app = typer.Typer()


async def get_parsed_commit_range(
    start: str, end: str, repo: str, gh: GitHubAPI, edit: bool = False
) -> Tuple[List[Commit], List[Commit]]:
    commits_iter = gh.getiter(f"/repos/{repo}/commits?sha={start}")

    commits = []
    unparsed_commits = []

    try:
        async for item in commits_iter:
            commit_hash = item["sha"]
            commit_message = item["commit"]["message"]
            if commit_hash == end:
                break

            invalid_message = False
            try:
                _default_parser(commit_message)
                # if this succeeds, do nothing
            except UnknownCommitMessageStyleError:
                print("Unknown commit message style!")
                if not commit_message.startswith("Merge"):
                    invalid_message = True
            if (
                (invalid_message or edit)
                and sys.stdout.isatty()
                and False
                and typer.confirm(f"Edit effective message '{commit_message}'?")
            ):
                commit_message = typer.edit(commit_message)
                _default_parser(commit_message)

            commit = Commit(commit_hash, commit_message, item["author"]["login"])
            commits.append(commit)

            if invalid_message:
                unparsed_commits.append(commit)

            print("-", commit)
            if len(commits) > 200:
                raise RuntimeError(f"{len(commits)} are a lot. Aborting!")
        return commits, unparsed_commits
    except gidgethub.BadRequest:
        print(
            "BadRequest for commit retrieval. That is most likely because you forgot to push the merge commit."
        )
        return


@app.command()
@make_sync
async def make_release(
    token: str = typer.Argument(..., envvar="GH_TOKEN"),
    draft: bool = True,
    dry_run: bool = False,
    edit: bool = False,
):
    async with aiohttp.ClientSession(loop=asyncio.get_event_loop()) as session:
        gh = GitHubAPI(session, __name__, oauth_token=token)

        version_file = Path("version_number")
        current_version = version_file.read_text()

        tag_hash = str(git("rev-list", "-n", "1", f"v{current_version}").strip())
        print("current_version:", current_version, "[" + tag_hash[:8] + "]")

        sha = git("rev-parse", "HEAD").strip()
        print("sha:", sha)

        repo = get_repo()
        print("repo:", repo)

        commits, _ = await get_parsed_commit_range(
            start=sha, end=tag_hash, repo=repo, gh=gh, edit=edit
        )

        bump = evaluate_version_bump(commits)
        print("bump:", bump)
        if bump is None:
            print("-> nothing to do")
            return
        next_version = get_new_version(current_version, bump)
        print("next version:", next_version)
        next_tag = f"v{next_version}"

        changes = generate_changelog(commits)
        md = markdown_changelog(next_version, changes, header=False)

        print(md)

        if not dry_run:
            version_file.write_text(next_version)
            git.add(version_file)

            zenodo_file = Path(".zenodo.json")
            update_zenodo(zenodo_file, repo, next_version)
            git.add(zenodo_file)

            citation_file = Path("CITATION.cff")
            update_citation(citation_file, next_version)
            git.add(citation_file)

            git.commit(m=f"Bump to version {next_tag}")

            target_hash = str(git("rev-parse", "HEAD")).strip()
            print("target_hash:", target_hash)

            git.push()

            commit_ok = False
            print("Waiting for commit", target_hash[:8], "to be received")
            for _ in range(RETRY_COUNT):
                try:
                    url = f"/repos/{repo}/commits/{target_hash}"
                    await gh.getitem(url)
                    commit_ok = True
                    break
                except InvalidField:
                    print("Commit", target_hash[:8], "not received yet")
                    pass  # this is what we want
                await asyncio.sleep(RETRY_INTERVAL)

            if not commit_ok:
                print("Commit", target_hash[:8], "was not created on remote")
                sys.exit(1)

            print("Commit", target_hash[:8], "received")

            await gh.post(
                f"/repos/{repo}/releases",
                data={
                    "body": md,
                    "tag_name": next_tag,
                    "name": next_tag,
                    "draft": draft,
                    "target_commitish": target_hash,
                },
            )


async def get_release_branch_version(
    repo: str, target_branch: str, gh: GitHubAPI
) -> str:
    content = await gh.getitem(
        f"repos/{repo}/contents/version_number?ref={target_branch}"
    )
    assert content["type"] == "file"
    return base64.b64decode(content["content"]).decode("utf-8")


async def get_tag_hash(tag: str, repo: str, gh: GitHubAPI) -> str:
    async for item in gh.getiter(f"repos/{repo}/tags"):
        if item["name"] == tag:
            return item["commit"]["sha"]
    raise ValueError(f"Tag {tag} not found")


async def get_merge_commit_sha(pr: int, repo: str, gh: GitHubAPI) -> str:
    for _ in range(RETRY_COUNT):
        pull = await gh.getitem(f"repos/{repo}/pulls/{pr}")
        if pull["mergeable"] is None:
            # no merge commit yet, wait a bit
            await asyncio.sleep(RETRY_INTERVAL)
            continue
        if not pull["mergeable"]:
            raise RuntimeError("Pull request is not mergeable, can't continue")
        return pull["merge_commit_sha"]
    raise RuntimeError("Timeout waiting for pull request merge status")


async def get_tag(tag: str, repo: str, gh: GitHubAPI):
    async for item in gh.getiter(f"repos/{repo}/tags"):
        if item["name"] == tag:
            return item
    return None


async def get_release(tag: str, repo: str, gh: GitHubAPI):
    existing_release = None
    try:
        existing_release = await gh.getitem(f"repos/{repo}/releases/tags/v{tag}")
    except gidgethub.BadRequest as e:
        if e.status_code == http.HTTPStatus.NOT_FOUND:
            pass  # this is what we want
        else:
            raise e
    return existing_release


@app.command()
@make_sync
async def pr_action(
    fail: bool = False,
    pr: int = None,
    token: Optional[str] = typer.Option(None, envvar="GH_TOKEN"),
    repo: Optional[str] = typer.Option(None, envvar="GH_REPO"),
):

    print("::group::Information")

    context = os.environ.get("GITHUB_CONTEXT")

    if context is not None:
        context = json.loads(context)
        repo = context["repository"]
        token = context["token"]
    else:
        if token is None or repo is None:
            raise ValueError("No context, need token and repo")
        if pr is None:
            raise ValueError("No context, need explicit PR to run on")

    async with aiohttp.ClientSession(loop=asyncio.get_event_loop()) as session:
        gh = GitHubAPI(session, __name__, oauth_token=token)

        if pr is not None:
            pr = await gh.getitem(f"repos/{repo}/pulls/{pr}")
        else:
            pr = context["event"]["pull_request"]

        target_branch = pr["base"]["ref"]
        print("Target branch:", target_branch)
        sha = pr["head"]["sha"]
        print("Source hash:", sha)

        merge_commit_sha = await get_merge_commit_sha(
            pr["number"],
            repo,
            gh,
        )
        print("Merge commit sha:", merge_commit_sha)

        # Get current version from target branch
        current_version = await get_release_branch_version(repo, target_branch, gh)
        tag_hash = await get_tag_hash(f"v{current_version}", repo, gh)
        print("current_version:", current_version, "[" + tag_hash[:8] + "]")

        commits, unparsed_commits = await get_parsed_commit_range(
            start=merge_commit_sha, end=tag_hash, repo=repo, gh=gh
        )

        bump = evaluate_version_bump(commits)
        print("bump:", bump)
        next_version = get_new_version(current_version, bump)
        print("next version:", next_version)
        next_tag = f"v{next_version}"

        print("::endgroup::")

        changes = generate_changelog(commits)
        md = markdown_changelog(next_version, changes, header=False)

        body = ""
        title = f"Release: {current_version} -> {next_version}"

        existing_release = await get_release(next_tag, repo, gh)
        existing_tag = await get_tag(next_tag, repo, gh)

        body += f"# `v{current_version}` -> `v{next_version}`\n"

        exit_code = 0

        if existing_release is not None or existing_tag is not None:

            if current_version == next_version:
                body += (
                    "## :no_entry_sign: Merging this will not result in a new version (no `fix`, "
                    "`feat` or breaking changes). I recommend **delaying** this PR until more changes accumulate.\n"
                )
                print("::warning::Merging this will not result in a new version")

            else:
                exit_code = 1
                title = f":no_entry_sign: {title}"
                if existing_release is not None:
                    body += f"## :warning: **WARNING**: A release for '{next_tag}' already exists"
                    body += f"[here]({existing_release['html_url']})** :warning:"
                    print(f"::error::A release for tag '{next_tag}' already exists")
                else:
                    body += (
                        f"## :warning: **WARNING**: A tag '{next_tag}' already exists"
                    )
                    print(f"::error::A tag '{next_tag}' already exists")

                body += "\n"
                body += ":no_entry_sign: I recommend to **NOT** merge this and double check the target branch!\n\n"

        else:
            body += f"## Merging this PR will create a new release `v{next_version}`\n"

        if len(unparsed_commits) > 0:
            body += "\n" * 3
            body += "## :warning: This PR contains commits which are not parseable:"
            for commit in unparsed_commits:
                msg, _ = commit.message.split("\n", 1)
                body += f"\n - {msg} {commit.sha})"
            body += "\n **Make sure these commits do not contain changes which affect the bump version!**"

        body += "\n\n"

        body += "### Changelog"

        body += md

        print("::group::PR message")
        print(body)
        print("::endgroup::")

        await gh.post(pr["url"], data={"body": body, "title": title})

        if fail:
            sys.exit(exit_code)


if __name__ == "__main__":
    app()
