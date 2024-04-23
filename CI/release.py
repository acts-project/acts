#!/usr/bin/env python3
import os
import asyncio
from typing import List, Optional, Tuple
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
from semantic_release.enums import LevelBump
from semantic_release.version import Version
from semantic_release.commit_parser.angular import (
    AngularCommitParser,
    AngularParserOptions,
)
from semantic_release.commit_parser.token import ParseError, ParseResult
import gidgethub
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

    # needed for semantic_release duck typing
    @property
    def hexsha(self):
        return self.sha


_default_parser = AngularCommitParser(AngularParserOptions())


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

        message: ParseResult = commit_parser.parse(commit)
        if isinstance(message, ParseError):
            print("Unknown commit message style!")
        else:
            changes.append(message.bump)

    if changes:
        level = max(changes)
        if level in LevelBump:
            bump = level
        else:
            print(f"Unknown bump level {level}")

    return bump


def generate_changelog(commits, commit_parser=_default_parser) -> dict:
    """
    Modified from: https://github.com/relekang/python-semantic-release/blob/48972fb761ed9b0fb376fa3ad7028d65ff407ee6/semantic_release/history/logs.py#L78
    """
    changes: dict = {"breaking": []}

    for commit in commits:
        message: ParseResult = commit_parser.parse(commit)

        if isinstance(message, ParseError):
            print("Unknown commit message style!")
            continue

        if message.type not in changes:
            changes[message.type] = list()

        capital_message = (
            message.descriptions[0][0].upper() + message.descriptions[0][1:]
        )
        changes[message.type].append((commit.sha, capital_message, commit.author))

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

            commit = Commit(commit_hash, commit_message, item["author"]["login"])

            invalid_message = False
            message: ParseResult = _default_parser.parse(commit)

            if isinstance(message, ParseError):
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
    force_next_version: Optional[str] = typer.Option(None, "--next-version"),
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

        current_version_obj = Version(*(map(int, current_version.split("."))))
        next_version_obj = current_version_obj.bump(bump)
        next_version = f"{next_version_obj.major}.{next_version_obj.minor}.{next_version_obj.patch}"
        if force_next_version is not None:
            next_version = force_next_version
        print("next version:", next_version)
        next_tag = f"v{next_version}"

        changes = generate_changelog(commits)
        md = markdown_changelog(next_version, changes, header=False)

        print(md)

        if not dry_run:
            execute_bump(next_version)

            git.commit(m=f"Bump to version {next_tag}", no_verify=True)

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


def execute_bump(next_version: str):
    version_file = Path("version_number")

    version_file.write_text(next_version)
    git.add(version_file)

    repo = get_repo()
    zenodo_file = Path(".zenodo.json")
    update_zenodo(zenodo_file, repo, next_version)
    git.add(zenodo_file)

    citation_file = Path("CITATION.cff")
    update_citation(citation_file, next_version)
    git.add(citation_file)


@app.command()
def bump(
    next_version: str = typer.Argument(..., help="Format: X.Y.Z"), commit: bool = False
):
    execute_bump(next_version)
    next_tag = f"v{next_version}"

    if commit:
        git.commit(m=f"Bump to version {next_tag}", no_verify=True)


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


if __name__ == "__main__":
    app()
