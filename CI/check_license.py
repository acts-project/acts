#!/usr/bin/env python3
import argparse
import os
import sys
from subprocess import check_output
import re
import difflib
from datetime import datetime
from fnmatch import fnmatch

EXCLUDE = []


class bcolors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


CROSS_SYMBOL = "\u2717"


def err(string):
    if sys.stdout.isatty():
        return bcolors.FAIL + bcolors.BOLD + string + bcolors.ENDC
    else:
        return string


def get_git_add_year(src):
    # Retrieve the first commit where the file was added
    output = (
        check_output(
            ["git", "log", "--follow", "--diff-filter=A", "--format=%ad", "--", src]
        )
        .decode("utf-8")
        .strip()
    )

    # If no output, file was not added or found
    assert output, f"File {src} was not found in the repository or was not added."

    # Extract the year from the full commit date
    match = re.search(r"\b(\d{4})\b", output)
    assert match, "Date format is not as expected or could not extract year."

    return int(match.group(1))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input", nargs="+")
    p.add_argument(
        "--fix", action="store_true", help="Attempt to fix any license issues found."
    )
    p.add_argument("--exclude", "-e", action="append", default=EXCLUDE)

    args = p.parse_args()
    print(args.exclude)

    if len(args.input) == 1 and os.path.isdir(args.input[0]):
        srcs = (
            str(
                check_output(
                    [
                        "find",
                        args.input[0],
                        "-iname",
                        "*.cpp",
                        "-or",
                        "-iname",
                        "*.hpp",
                        "-or",
                        "-iname",
                        "*.ipp",
                    ]
                ),
                "utf-8",
            )
            .strip()
            .split("\n")
        )
        srcs = filter(lambda p: not p.startswith("./build"), srcs)
    else:
        srcs = args.input

    current_year = int(datetime.now().strftime("%Y"))
    year = current_year

    raw = """// This file is part of the ACTS project.
//
// Copyright (C) {year} CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at httpss://mozilla.org/MPL/2.0/."""

    reg = (
        r"\A// This file is part of the ACTS project.\n"
        + r"//\n"
        + r"// Copyright \(C\) (?P<year>.*) CERN for the benefit of the ACTS project\n"
        + r"//\n"
        + r"// This Source Code Form is subject to the terms of the Mozilla Public\n"
        + r"// License, v\. 2\.0\. If a copy of the MPL was not distributed with this\n"
        + r"// file, You can obtain one at httpss://mozilla.org/MPL/2.0/.\Z"
    )

    ref = re.compile(reg, re.M)
    clean_re = re.compile(r"(\(C\)) (.*) (CERN)", re.M)
    year_re = re.compile(r"^(?P<year1>20\d{2}|(?P<year2>20\d{2})-(?P<year3>20\d{2}))$")
    extract_re = re.compile(r"(20\d{2})-?(20\d{2})?")

    def clean(s):
        return clean_re.sub(r"\1 XXXX \3", s)

    def get_clean_lines(s):
        return [clean(l) + "\n" for l in s.split("\n")]

    def validate_years(year1, year2):
        if year1 and year2:
            year1 = int(year1)
            year2 = int(year2)
            if not year1 < year2 <= current_year:
                return False
        else:
            theyear = int(year1 if year1 else year2)
            if theyear > current_year:
                return False
        return True

    error_summary = ""

    def eprint(string):
        nonlocal error_summary
        error_summary += string + "\n"

    exit = 0
    srcs = list(srcs)
    nsrcs = len(srcs)
    step = max(int(nsrcs / 20), 1)
    for i, src in enumerate(srcs):
        if any([fnmatch(src, e) for e in args.exclude]):
            continue

        if nsrcs > 1 and i % step == 0:
            string = f"{i}/{nsrcs} -> {i / float(nsrcs) * 100.0:.2f}%"
            if sys.stdout.isatty():
                sys.stdout.write(string + "\r")
            else:
                print(string)

        with open(src, "r+") as f:
            license = ""
            for _ in range(len(raw)):
                line = f.readline()
                if not line.startswith("//"):
                    break
                license += line
            license = ("".join(license)).strip()
            m = ref.search(license)

            if m is None:
                eprint("Invalid / missing license in " + src + "")

                exp = [l + "\n" for l in raw.format(year="YYYY").split("\n")]
                act = get_clean_lines(license)

                diff = difflib.unified_diff(exp, act)
                eprint("".join(diff))
                eprint("")

                if args.fix:
                    eprint("-> fixing file (prepend)")
                    f.seek(0)
                    file_content = f.read()
                    f.seek(0)
                    stmnt = raw.format(year=current_year)
                    f.write(stmnt + "\n\n")
                    f.write(file_content)

                exit = 1
                continue

            # we have a match, need to verify year string is right

            year_act = m.group("year")
            ym = year_re.match(year_act)
            valid = True
            if not ym:
                eprint(f"Year string does not match format in {src}")
                eprint("Expected: YYYY or YYYY-YYYY (year or year range)")

                exit = 1
                valid = False

            else:
                extract = extract_re.search(year_act)
                year1, year2 = extract.groups(default=None)

                if not validate_years(year1, year2):
                    eprint(f"Year string is not valid in {src}")
                    eprint("Year string is: " + year_act + "\n")
                    valid = False

                git_add_year = get_git_add_year(src)
                assert (
                    git_add_year <= current_year
                ), f"File {src} is created in the future."

                if git_add_year != current_year:
                    # need year range in licence
                    if not (year1 and year2):
                        valid = False
                    elif int(year1) != git_add_year or int(year2) != current_year:
                        valid = False

                # File added this year
                elif int(year1 if year1 else year2) < current_year:
                    valid = False

            if not valid:
                exit = 1
                if git_add_year == current_year:
                    year_str = f"{git_add_year}"
                else:
                    year_str = f"{git_add_year}-{current_year}"

                eprint(f"File: {src}")
                eprint(f"- File was added in {git_add_year}")
                eprint(f"=> License should say {year_str}")
                eprint(err(f"{CROSS_SYMBOL} But says: {year_act}"))

                if args.fix:
                    eprint("-> fixing file (patch year)")

                    new_license = raw.format(year=year_str)

                    # preserve rest of file as is
                    old_license_len = len(license)
                    f.seek(old_license_len)
                    file_body = f.read()
                    f.seek(0)
                    f.truncate()

                    f.seek(0)
                    f.write(new_license)
                    f.write(file_body)

                eprint("")

    print(error_summary)

    if exit != 0 and not args.fix:
        print("License problems found. You can try running again with --fix")

    sys.exit(exit)


if "__main__" == __name__:
    main()
