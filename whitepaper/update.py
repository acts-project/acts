#!/usr/bin/env python3
import base64
from dataclasses import dataclass
from pathlib import Path
from typing import List
import re
import textwrap

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


@app.command()
def pull(config_file: Path = Path(__file__).parent / "config.toml"):
    console = rich.console.Console()

    with console.status("[bold green]Loading config...") as status:
        config = Config(**toml.load(config_file))

        for whp in config.whitepapers:
            status.update(f"Loading {whp.repository}...")
            m = re.match(r"^https://github.com/(.*/.*)$", whp.repository)

            repo = m.group(1)

            r = requests.get(f"{API_URL}/repos/{repo}")
            repo_data = r.json()

            r = requests.get(f"{API_URL}/repos/{repo}/contents/metadata.tex")
            metadata_data = r.json()
            content = base64.b64decode(metadata_data["content"]).decode("utf-8")

            metadata = parse_metadata(content)
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
    target_file = Path(__file__).parent.parent / "docs" / "whitepapers.md"

    tpl = Template(template_file.read_text())

    target_file.write_text(tpl.render(config=config))

    #  for whp in config.whitepapers:


if "__main__" == __name__:
    app()
