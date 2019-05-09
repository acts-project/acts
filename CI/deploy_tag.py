#!/usr/bin/env python3
import argparse
import os
import fs
from fs.osfs import OSFS
from markdown import markdown as mdlib
import textwrap
import json

from util import smtp, def_arguments, get_lxplus_fs, Spinner, gitlab
from release_notes import collect_milestone, make_release_notes, parse_version

sender_email = "Acts Bot <atsjenkins@cern.ch>"
receiver_email = os.getenv("TAG_DEPLOY_EMAIL_RECEIVERS")

from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart


def markdown(md):
    return mdlib(md, extensions=["pymdownx.extra"])


def main():
    p = argparse.ArgumentParser()
    p = def_arguments(p, acc=True, gl=True)

    p.add_argument("--doc-source", required=True)
    p.add_argument("--dry-run", action="store_true")
    p.add_argument(
        "--ref", default=os.getenv("CI_COMMIT_REF_NAME", None), required=True
    )
    p.add_argument(
        "--doc-root",
        default=os.getenv("DOC_WEBSITE_ROOT", "/eos/user/a/atsjenkins/www/ACTS/"),
    )
    p.add_argument(
        "--doc-public-url",
        default=os.getenv("DOC_WEBSITE_URL", "https://acts.web.cern.ch/ACTS/"),
    )

    args = p.parse_args()

    src_fs = OSFS(os.path.abspath(args.doc_source))

    www_fs = get_lxplus_fs(args).opendir(os.path.join(args.doc_root))

    if not www_fs.exists(args.ref):
        www_fs.makedirs(os.path.join(args.ref, "doc"))
    refdir = www_fs.opendir(os.path.join(args.ref, "doc"))

    # refdir = OSFS("/tmp/doctest")

    with Spinner(f"Publishing doc for {args.ref}"):
        if not args.dry_run:
            fs.copy.copy_dir(src_fs, ".", refdir, ".")

    doc_url = os.path.join(args.doc_public_url, args.ref, "doc")
    print("Doc is available at", doc_url)

    # write tag info json file
    if not args.dry_run:
        with www_fs.open("latest_release.json", "w") as f:
            json.dump({"subject": "release", "status": args.ref, "color": "yellow"}, f)

    gl = gitlab(args)
    project = gl.projects.get("acts/acts-core")

    version = parse_version(args.ref)
    with Spinner(text="Loading milestone"):
        milestones = project.milestones.list(all=True)
        milestone = None
        for ms in milestones:
            if ms.title == version:
                milestone = ms
                break

    relnotes = make_release_notes(milestone, *collect_milestone(milestone))

    message = MIMEMultipart("alternative")
    message["Subject"] = f"New Acts release: {args.ref}"
    message["From"] = sender_email
    message["To"] = receiver_email

    text = """
    Dear Acts enthusiasts,

    a new tag '{ref}' of the Acts project has been created.

    You can get the source code from git using:

    git clone https://gitlab.cern.ch/acts/acts-core.git
    cd acts-core/
    git checkout {ref}

    or download a tarball with the source from

    https://gitlab.cern.ch/acts/acts-core/-/archive/{ref}/acts-core-{ref}.tar.gz

    The documentation is deployed at
    https://acts.web.cern.ch/ACTS/{ref}/doc/index.html

    Cheers,
    your friendly Acts robot
    """
    text = textwrap.dedent(text).format(ref=args.ref, relnotes=relnotes)

    md = """
    Dear Acts enthusiasts,

    a new tag of the Acts project has been created.

    ---

    # {ref}
    [![](https://badgen.net/badge/release/{ref}/yellow)](https://gitlab.cern.ch/acts/acts-core/tags/{ref})
    {relnotes}

    ---

    You can get the source code from git using:

    ```bash
    git clone https://gitlab.cern.ch/acts/acts-core.git
    cd acts-core/
    git checkout {ref}
    ```

    or download a tarball with the source from

    https://gitlab.cern.ch/acts/acts-core/-/archive/{ref}/acts-core-{ref}.tar.gz

    The documentation is deployed at

    https://acts.web.cern.ch/ACTS/{ref}/doc/index.html

    Cheers,<br/>
    your friendly Acts robot
    """

    md = textwrap.dedent(md).format(
        ref=args.ref, relnotes=relnotes, commit=os.getenv("CI_COMMIT_SHA")
    )

    html = """\
    <html>
      <body>
        {text}
      </body>
    </html>
    """.format(
        text=markdown(textwrap.dedent(md))
    )

    # print(html)

    part1 = MIMEText(text, "plain")
    part2 = MIMEText(html, "html")

    message.attach(part1)
    message.attach(part2)

    with Spinner("Sending email"):
        if not args.dry_run:
            with smtp(args) as server:
                server.sendmail(sender_email, receiver_email, message.as_string())


if "__main__" == __name__:
    main()
