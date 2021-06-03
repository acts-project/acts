#!/usr/bin/env python3
import os
import asyncio
import subprocess
from typing import List, Optional
import re
from pathlib import Path
import sys
import http

import aiohttp
from gidgethub.aiohttp import GitHubAPI
from gidgethub import InvalidField
import gidgethub
from semantic_release.history import angular_parser, get_new_version
from semantic_release.errors import UnknownCommitMessageStyleError
from semantic_release.history.logs import LEVELS
from semantic_release.history.parser_helpers import ParsedCommit
import sh
import click
from dotenv import load_dotenv

load_dotenv()

git = sh.git


def run(cmd):
    return subprocess.check_output(cmd).decode("utf-8").strip()


def get_repo():
    # origin = run(["git", "remote", "get-url", "origin"])
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

    def __init__(self, sha: str, message: str):
        self.sha = sha
        self.message = self._normalize(message)

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
        except UnknownCommitMessageStyleError as err:
            pass

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
            changes[message.type].append((commit.sha, capital_message))

            if message.breaking_descriptions:
                for paragraph in message.breaking_descriptions:
                    changes["breaking"].append((commit.sha, paragraph))
            elif message.bump == 3:
                changes["breaking"].append((commit.sha, message.descriptions[0]))

        except UnknownCommitMessageStyleError as err:
            pass

    return changes


def markdown_changelog(version: str, changelog: dict, header: bool = False) -> str:
    output = f"## v{version}\n" if header else ""

    for section, items in changelog.items():
        if len(items) == 0:
            continue
        output += "\n### {0}\n".format(section.capitalize())

        for item in items:
            output += "* {0} ({1})\n".format(item[1], item[0])

    return output


async def main(draft, dry_run):
    token = os.environ["GH_TOKEN"]
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

        commits_iter = gh.getiter(f"/repos/{repo}/commits?sha={sha}")

        commits = []

        try:
          async for item in commits_iter:
              commit_hash = item["sha"]
              commit_message = item["commit"]["message"]
              if commit_hash == tag_hash:
                  break

              try:
                  _default_parser(commit_message)
                  # if this succeeds, do nothing
              except UnknownCommitMessageStyleError as err:
                print("Unkown commit message style:")
                print(commit_message)
                if sys.stdout.isatty() and click.confirm("Edit effective message?"):
                  commit_message = click.edit(commit_message)
                  _default_parser(commit_message)

              commit = Commit(commit_hash, commit_message)
              commits.append(commit)
              print("-", commit)
        except gidgethub.BadRequest:
          print("BadRequest for commit retrieval. That is most likely because you forgot to push the merge commit.")
          return

        if len(commits) > 100:
            print(len(commits), "are a lot. Aborting!")
            sys.exit(1)

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
          git.commit(m=f"Bump to version {next_tag}")

          # git.tag(next_tag)
          target_hash = str(git("rev-parse", "HEAD")).strip()
          print("target_hash:", target_hash)

          git.push()

          commit_ok = False
          print("Waiting for commit", target_hash[:8], "to be received")
          for _ in range(10):
              try:
                url = f"/repos/{repo}/commits/{target_hash}"
                await gh.getitem(url)
                commit_ok = True
                break
              except InvalidField as e:
                  print("Commit", target_hash[:8], "not received yet")
                  pass # this is what we want
              await asyncio.sleep(0.5)

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


@click.command()
@click.option("--draft/--no-draft", default=True)
@click.option("--dry-run/--no-dry-run", default=False)
def main_sync(*args, **kwargs):
    loop = asyncio.get_event_loop()
    loop.run_until_complete(main(*args, **kwargs))


if __name__ == "__main__":
    main_sync()
