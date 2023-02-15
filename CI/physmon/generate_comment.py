#!/usr/bin/env python3


### CHANGE BELOW WITH CAUTION ###

import sys
import os
from pathlib import Path, PurePath

import jinja2

template_source = (Path(__file__).parent / "comment_template.md").read_text()

macro_source = """
{% macro detail_block(title, check) %}
{% if exists(check) %}
<details>
  <summary><b>{{ title }}</b></summary>
  {{ caller() }}
</details>
{% else %}
<b>:x: {{ title }}</b>
{% endif %}
{% endmacro %}
"""

template_source = macro_source + template_source

template = jinja2.Template(template_source)

_, artifact_dir, outfile = sys.argv
artifact_dir = Path(artifact_dir)
outfile = Path(outfile)

artifact_url = os.environ["ARTIFACT_URL"]
pr_sha = os.environ["PR_SHA"]

print("artifact_url:", artifact_url)
print("pr_sha:", pr_sha)

has_errors = False


def exists(arg):
    file = artifact_dir / arg
    result = file.exists()
    print("Check exists:", file, "=>", result)
    global has_errors
    if not result:
        has_errors = True
    return result


def all_exist(*args):
    result = all(exists(a) for a in args)
    global has_errors
    if not result:
        has_errors = True
    return result


def make_url(title, target):
    if exists(target):
        url = artifact_url + "/" + target
        return f"[{title}]({url})"
    else:
        global has_errors
        has_errors = True
        return f":x: {title}"


def make_image(target, width):
    file = target
    if "?" in target:
        file = target[: target.index("?")]
    if exists(file):
        url = artifact_url + "/" + file
        return f'<img src="{url}?to_png=1" width="{width}"/>'
    else:
        global has_errors
        has_errors = True
        return f":framed_picture: {target} :x:"


kwargs = dict(
    url=artifact_url,
    commit=pr_sha,
    exists=exists,
    all_exist=all_exist,
    make_url=make_url,
    make_image=make_image,
    has_errors=has_errors,
)

# render once to fill `has_errors`
template.render(**kwargs)

kwargs["has_errors"] = has_errors

outfile.write_text(template.render(**kwargs))
