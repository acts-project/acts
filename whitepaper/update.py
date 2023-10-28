#!/usr/bin/env python3
import base64
from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess
from subprocess import run
from typing import List
import re
import textwrap
import tempfile

import typer
import toml
import pydantic
import requests
import rich.console
import rich.panel
import rich.pretty
from jinja2 import Template


app = typer.Typer()


class MetaData(pydantic.BaseModel):
    authors: List[str]
    title: str
    description: str


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


API_URL = "https://api.github.com"


def parse_metadata(content: str) -> MetaData:
    m = re.search(r"\\author{(.*)}", content)
    assert m is not None
    author_text = m.group(1)
    authors = [a.strip() for a in author_text.split(r"\and")]

    m = re.search(r"\\title{(.*)}", content)
    assert m is not None
    title = m.group(1).strip()

    m = re.search(
        r"^%% BEGIN SHORT DESCRIPTION$(.*)^%% END SHORT DESCRIPTION$",
        content,
        re.MULTILINE | re.DOTALL,
    )
    assert m is not None
    description = m.group(1)
    description = textwrap.dedent(description).strip()

    return MetaData(authors=authors, title=title, description=description)


def which(cmd: str) -> Path:
    exe = (
        run(["command", "-v", cmd], check=True, capture_output=True)
        .stdout.decode()
        .strip()
    )
    return Path(exe)


@app.command()
def pull(config_file: Path = Path(__file__).parent / "config.toml"):
    console = rich.console.Console()

    git = which("git")
    latexmk = which("latexmk")
    convert = which("convert")

    with console.status("[bold green]Loading config...") as status:
        config = Config(**toml.load(config_file))

        for whp in config.whitepapers:
            status.update(f"Loading {whp.repository}...")

            repo = whp.repo
            slug = whp.slug

            r = requests.get(f"{API_URL}/repos/{repo}")
            repo_data = r.json()

            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir = Path(tmpdir)
                checkout = tmpdir / slug
                status.update(f"\[{repo}] Cloning repo {repo_data['clone_url']}...")
                run(
                    [git, "clone", repo_data["clone_url"], checkout],
                    check=True,
                    cwd=tmpdir,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )

                status.update(f"\[{repo}] Compiling PDF...")
                run(
                    [latexmk],
                    check=True,
                    cwd=checkout,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )

                build_dir = checkout / "build"
                output_file = build_dir / "main.pdf"
                assert output_file.exists(), "PDF not found"

                title_page = build_dir / "title_page.png"

                status.update(f"\[{repo}] Converting PDF to PNG...")
                run(
                    [
                        convert,
                        "-density",
                        "300",
                        f"{output_file}[0]",
                        "-background",
                        "white",
                        "-alpha",
                        "remove",
                        "-alpha",
                        "off",
                        "-depth",
                        "2",
                        str(title_page),
                    ],
                    cwd=tmpdir,
                    check=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                )
                assert title_page.exists(), "PNG not found"
                shutil.copyfile(title_page, config_file.parent / f"{slug}.png")

                metadata_file = checkout / "metadata.tex"

                metadata = parse_metadata(metadata_file.read_text())

            console.print(
                rich.panel.Panel(rich.pretty.Pretty(metadata), title=whp.repository)
            )

            whp.metadata = metadata

        status.update("Updating config...")

        with open(config_file, "w") as f:
            toml.dump(config.model_dump(), f)


@app.command()
def render(config_file: Path = Path(__file__).parent / "config.toml"):
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
