#!/usr/bin/env python3
from dataclasses import dataclass
from typing import List, Dict, Any
from pathlib import Path
import shlex
import asyncio
import functools
import os
import click

import typer
import gidgethub
from gidgethub.aiohttp import GitHubAPI
import aiohttp


def wrap_async(fn):
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        return asyncio.run(fn(*args, **kwargs))

    return wrapper


class CommandError(Exception):
    pass


@dataclass
class Context:
    pr: Dict[str, Any]
    sender: str
    github_token: str


@click.group()
def app():
    pass


@app.group()
def run_experiment():
    pass


@run_experiment.command()
@click.option("--revert-sha", "-r", multiple=True)
@click.pass_obj
@wrap_async
async def atlas(ctx: Context, revert_sha: List[str]):
    gitlab_trigger_token = os.environ["GITLAB_TRIGGER_TOKEN"]
    gitlab_trigger_url = os.environ["GITLAB_TRIGGER_URL"]
    async with aiohttp.ClientSession() as session:
        gh = GitHubAPI(session, "acts-commands", oauth_token=ctx.github_token)

        pr = ctx.pr

        head_clone_url = pr["head"]["repo"]["clone_url"]
        head_branch = pr["head"]["ref"]
        head_sha = pr["head"]["sha"]

        variable_summary = f"""
| Variable | Value |
|------|------|
| `ACTS_GIT_REPO` | {head_clone_url} |
| `ACTS_REF` | `{head_branch}` |
| `SOURCE_SHA` | {head_sha} |
| `REVERT_SHAS` | {",".join(revert_sha)} |
        """

        body = f"""
@{ctx.sender}
🟡 I'm going to trigger an ATLAS experiment pipeline for you:

{variable_summary}
        """
        comment = await gh.post(pr["comments_url"], data={"body": body})

        variables = {
            "ACTS_GIT_REPO": head_clone_url,
            "ACTS_REF": head_branch,
            "SOURCE_SHA": head_sha,
            "PR_URL": pr["url"],
            "REVERT_SHAS": ",".join(revert_sha),
            "REPORT_COMMENT_URL": comment["url"],
        }
        data = {
            "token": gitlab_trigger_token,
            "ref": "main",
            **{f"variables[{k}]": v for k, v in variables.items()},
        }
        print(gitlab_trigger_url)
        print(data)
        async with session.post(
            url=gitlab_trigger_url,
            data=data,
        ) as resp:
            if resp.status != 201:
                body = f"""
@{ctx.sender}
🔴 I'm sorry, I couldn't run your command because of an error:
```
{await resp.text()}
```
{variable_summary}
                """
                await gh.post(comment["url"], data={"body": body})

                return

            data = await resp.json()
            pipeline_url = data["web_url"]

        body = f"""
@{ctx.sender}
🟡 I triggered an ATLAS experiment [pipeline]({pipeline_url}) for you

{variable_summary}
        """
        await gh.post(comment["url"], data={"body": body})


async def get_author_in_team(gh: GitHubAPI, author: str, allow_team: str) -> bool:
    allow_org, allow_team = allow_team.split("/", 1)

    try:
        membership = await gh.getitem(
            f"/orgs/{allow_org}/teams/{allow_team}/memberships/{author}"
        )
        return True
    except gidgethub.BadRequest as e:
        if e.status_code != 404:
            raise e

    return False


async def preflight(
    token: str, pr_url: str, sender: str, repository: str, allow_team: str
):
    async with aiohttp.ClientSession() as session:
        gh = GitHubAPI(session, "acts-commands", oauth_token=token)

        if not await get_author_in_team(gh, sender, allow_team):
            raise RuntimeError(f"{sender} is not in {allow_team}")

        return await gh.getitem(pr_url)


async def report_error(token: str, pr: Dict[str, Any], sender: str, error: Exception):
    async with aiohttp.ClientSession() as session:
        gh = GitHubAPI(session, "acts-commands", oauth_token=token)

        body = f"""
@{sender}
🔴 I'm sorry, I couldn't run your command because of an error:
```
{error}
```
"""
        await gh.post(pr["comments_url"], data={"body": body})


def main(
    pr: str = typer.Option(),
    body: str = typer.Option(),
    sender: str = typer.Option(),
    repository: str = typer.Option(),
    allow_team: str = typer.Option("acts-project/ci-perms", envvar="ALLOW_TEAM"),
):
    if Path(body).exists():
        body = Path(body).read_text().strip()

    if len(body.split("\n")) > 1:
        raise typer.BadParameter("Body must be a single line")

    if not body.startswith("/"):
        raise typer.BadParameter("Body must start with a slash")
    body = body[1:]

    args = shlex.split(body)

    token = os.environ["GITHUB_TOKEN"]
    pr = asyncio.run(preflight(token, pr, sender, repository, allow_team))

    try:
        app(
            args,
            obj=Context(pr=pr, github_token=token, sender=sender),
            standalone_mode=False,
        )
    except (CommandError, click.exceptions.ClickException) as e:
        asyncio.run(report_error(token, pr, sender, e))


typer.run(main)
