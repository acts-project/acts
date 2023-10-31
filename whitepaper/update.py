#!/usr/bin/env python3
from pathlib import Path
import shutil
import base64
import subprocess
from subprocess import run
from typing import List
import re
import textwrap
import tempfile
from typing import Optional
import functools
import asyncio

import typer
import toml
import pydantic
import rich.console
import rich.panel
import rich.pretty
from jinja2 import Template
import aiohttp
from gidgethub.aiohttp import GitHubAPI
from fsspec.implementations.zip import ZipFileSystem


app = typer.Typer()


class MetaData(pydantic.BaseModel):
    authors: List[str]
    title: str
    description: str = ""


class WhitePaper(pydantic.BaseModel):
    repository: str
    metadata: MetaData | None = None

    @property
    def repo(self) -> str:
        m = re.match(r"^https://github.com/(.*/.*)$", self.repository)
        assert m is not None
        return m.group(1)

    @property
    def slug(self) -> str:
        owner, name = self.repo.split("/")
        return f"{owner}_{name}"


class Config(pydantic.BaseModel):
    whitepapers: List[WhitePaper]


def parse_metadata(content: str) -> MetaData:
    m = re.search(r"\\author{(.*)}", content)
    assert m is not None
    author_text = m.group(1)
    authors = [a.strip() for a in author_text.split(r"\and")]

    m = re.search(r"\\title{(.*)}", content)
    assert m is not None
    title = m.group(1).strip()

    return MetaData(authors=authors, title=title)


def which(cmd: str) -> Optional[Path]:
    try:
        exe = (
            run(["command", "-v", cmd], check=True, capture_output=True)
            .stdout.decode()
            .strip()
        )
        return Path(exe)
    except subprocess.CalledProcessError:
        return None


def coro(fn):
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        return asyncio.run(fn(*args, **kwargs))

    return wrapper


def make_titlepage_image(convert: Path, pdf: Path, output: Path):
    run(
        [
            convert,
            "-density",
            "300",
            f"{pdf}[0]",
            "-background",
            "white",
            "-alpha",
            "remove",
            "-alpha",
            "off",
            "-depth",
            "2",
            str(output),
        ],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )


async def get_file(gh: GitHubAPI, repo: str, file: str) -> str:
    r = await gh.getitem(f"/repos/{repo}/contents/{file}")
    return base64.b64decode(r["content"]).decode()


@app.command()
@coro
async def pull(
    config_file: Path = Path(__file__).parent / "whitepapers.toml",
    github_token: str = typer.Option(..., envvar="GITHUB_TOKEN"),
):
    console = rich.console.Console()
    convert = which("convert")
    if convert is None:
        print("convert (imagemagick) not found, please install")
        raise typer.Exit(1)

    async with aiohttp.ClientSession() as session:
        gh = GitHubAPI(session, "requester", oauth_token=github_token)

        with console.status("[bold green]Loading config...") as status:
            config = Config(**toml.load(config_file))

            for whp in config.whitepapers:
                status.update(f"Loading {whp.repository}...")

                repo = whp.repo
                slug = whp.slug

                latest_artifact = None

                async for artifact in gh.getiter(
                    f"/repos/{repo}/actions/artifacts", iterable_key="artifacts"
                ):
                    if artifact["workflow_run"]["head_branch"] == "main":
                        latest_artifact = artifact
                        break
                if latest_artifact is None:
                    print("No artifacts found for", whp.repository)
                    raise typer.Exit(1)

                status.update(f"\[{repo}] Downloading artifact...")

                r = await session.get(
                    latest_artifact["archive_download_url"],
                    headers={"Authorization": f"Token {github_token}"},
                )

                with tempfile.TemporaryFile() as fh, tempfile.NamedTemporaryFile(
                    suffix=".pdf"
                ) as pdf_fh:
                    async for data, _ in r.content.iter_chunks():
                        fh.write(data)
                    fh.seek(0)

                    zipfs = ZipFileSystem(fh)
                    files = zipfs.ls("/", detail=False)
                    if len(files) != 1:
                        print("Unexpected number of files in artifact", files)
                        raise typer.Exit(1)
                    if not files[0].endswith(".pdf"):
                        print("Unexpected file in artifact", files)
                        raise typer.Exit(1)

                    shutil.copyfileobj(zipfs.open(files[0]), pdf_fh)

                    status.update(f"\[{repo}] Converting PDF to PNG...")
                    title_page = config_file.parent / f"{slug}.png"
                    make_titlepage_image(convert, Path(pdf_fh.name), title_page)

                status.update(f"\[{repo}] Getting metadata...")

                metadata, abstract = await asyncio.gather(
                    get_file(gh, repo, "metadata.tex"),
                    get_file(gh, repo, "abstract.tex"),
                )

                abstract = textwrap.dedent(abstract).strip()
                abstract = "\n".join(
                    [
                        line
                        for line in abstract.split("\n")
                        if not line.strip().startswith("%")
                    ]
                ).strip()
                metadata = parse_metadata(metadata)
                metadata.description = abstract

                console.print(
                    rich.panel.Panel(rich.pretty.Pretty(metadata), title=whp.repository)
                )

                whp.metadata = metadata

            status.update("Updating config...")

            with open(config_file, "w") as f:
                toml.dump(config.model_dump(), f)


@app.command()
def render(config_file: Path = Path(__file__).parent / "whitepapers.toml"):
    config = Config(**toml.load(config_file))
    template_file = Path(__file__).parent / "template.md"
    docs_path = Path(__file__).parent.parent / "docs" / "whitepapers"
    docs_path.mkdir(parents=True, exist_ok=True)
    target_file = docs_path / "whitepapers.md"
    image_path = Path(__file__).parent.parent / "docs" / "figures" / "whitepapers"
    image_path.mkdir(parents=True, exist_ok=True)

    tpl = Template(template_file.read_text())

    target_file.write_text(
        tpl.render(config=config, image_path="../figures/whitepapers")
    )

    for whp in config.whitepapers:
        shutil.copyfile(
            config_file.parent / f"{whp.slug}.png",
            f"{whp.slug}.png",
        )


if "__main__" == __name__:
    app()
