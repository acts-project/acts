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


class CommitInfo:
    date = None
    year = None
    author = None
    subject = None
    body = None


def check_git_dates(src):
    output = (
        check_output(["git", "log", "--format={{{%an|%ad|%s|%b}}}", "--", src])
        .decode("utf-8")
        .strip()
    )

    # find single outputs
    commits = re.findall(r"{{{((?:.|\n)*?)}}}", output)
    commits = [c for c in commits if "[ignore-license]" not in c]
    commits = [c.split("|") for c in commits]

    # print(output)
    mod = commits[0]
    add = commits[-1]

    madd = re.match(r".*\d{2}:\d{2}:\d{2} (\d{4})", add[1])
    assert madd != None, "Regex did not match git log output"
    mmod = re.match(r".*\d{2}:\d{2}:\d{2} (\d{4})", mod[1])
    assert mmod != None, "Regex did not match git log output"

    addcommit = CommitInfo()
    addcommit.date = add[1]
    addcommit.year = int(madd.group(1))
    addcommit.author = add[0]
    addcommit.subject = add[2]
    addcommit.body = add[3]

    modcommit = CommitInfo()
    modcommit.date = mod[1]
    modcommit.year = int(mmod.group(1))
    modcommit.author = mod[0]
    modcommit.subject = mod[2]
    modcommit.body = mod[3]

    return addcommit, modcommit


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument(
        "--fix", action="store_true", help="Attempt to fix any license issues found."
    )
    p.add_argument(
        "--check-years",
        action="store_true",
        help="Check the license year info using git info for each file.",
    )
    p.add_argument(
        "--fail-year-mismatch",
        action="store_true",
        help="Fail if year in license statement is not valid.",
    )
    p.add_argument("--exclude", "-e", action="append", default=EXCLUDE)

    args = p.parse_args()
    print(args.exclude)

    if os.path.isdir(args.input):
        srcs = (
            str(
                check_output(
                    [
                        "find",
                        args.input,
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
        srcs = [args.input]

    year = int(datetime.now().strftime("%Y"))

    raw = """// This file is part of the Acts project.
//
// Copyright (C) {year} CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/."""

    reg = (
        r"\A// This file is part of the Acts project.\n"
        + r"//\n"
        + r"// Copyright \(C\) (?P<year>.*) CERN for the benefit of the Acts project\n"
        + r"//\n"
        + r"// This Source Code Form is subject to the terms of the Mozilla Public\n"
        + r"// License, v\. 2\.0\. If a copy of the MPL was not distributed with this\n"
        + r"// file, You can obtain one at http://mozilla.org/MPL/2.0/.\Z"
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
            if year1 >= year2:
                return False
            if year1 > year or year2 > year:
                return False
        else:
            theyear = int(year1 if year1 else year2)
            if theyear > year:
                return False
        return True

    error_summary = ""
    info_summary = ""

    def eprint(*args):
        nonlocal error_summary
        error_summary += " ".join(map(str, args)) + "\n"

    def year_print(*pargs):
        nonlocal error_summary
        nonlocal info_summary
        string = " ".join(map(str, pargs)) + "\n"
        if args.fail_year_mismatch:
            error_summary += string
        else:
            info_summary += string

    exit = 0
    srcs = list(srcs)
    nsrcs = len(srcs)
    step = int(nsrcs / 20)
    for i, src in enumerate(srcs):
        if any([fnmatch(src, e) for e in args.exclude]):
            continue

        if nsrcs > 1 and i % step == 0:
            string = "{}/{} -> {:.2f}%".format(i, nsrcs, i / float(nsrcs) * 100.0)
            if sys.stdout.isatty():
                sys.stdout.write(string + "\r")
            else:
                print(string)

        with open(src, "r+") as f:
            license = ""
            for x in range(len(raw)):
                line = f.readline()
                if not line.startswith("//"):
                    break
                license += line
            license = ("".join(license)).strip()
            m = ref.search(license)

            if m == None:
                eprint("Invalid / missing license in " + src + "")

                exp = [l + "\n" for l in raw.format(year="XXXX").split("\n")]
                act = get_clean_lines(license)

                diff = difflib.unified_diff(exp, act)
                eprint("".join(diff))
                eprint()

                if args.fix:
                    eprint("-> fixing file (prepend)")
                    f.seek(0)
                    file_content = f.read()
                    f.seek(0)
                    stmnt = raw.format(year=year)
                    f.write(stmnt + "\n\n")
                    f.write(file_content)

                exit = 1
            else:
                # we have a match, need to verify year string is right

                if args.check_years:
                    git_add_commit, git_mod_commit = check_git_dates(src)
                    git_add_year = git_add_commit.year
                    git_mod_year = git_mod_commit.year
                year_act = m.group("year")
                ym = year_re.match(year_act)
                valid = True
                if not ym:
                    eprint("Year string does not match format in {}".format(src))
                    eprint("Expected: YYYY or YYYY-YYYY (year or year range)")
                    eprint("Actual:   {}\n".format(year_act))

                    if args.fix:
                        extract = extract_re.search(year_act)
                        year1 = extract.group(1)
                        year2 = extract.group(2)

                    exit = 1
                    valid = False

                else:
                    extract = extract_re.search(year_act)
                    year1 = extract.group(1)
                    year2 = extract.group(2)

                    if not validate_years(year1, year2):
                        eprint("Year string is not valid in {}".format(src))
                        eprint("Year string is: " + year_act + "\n")
                        exit = 1
                        valid = False

                    if args.check_years:
                        if git_add_year != git_mod_year:
                            # need year range in licence
                            if not (year1 and year2):
                                year_print("File: {}".format(src))
                                # year_print("o File was modified in a different year ({}) than it was added ({})."
                                # .format(git_mod_year, git_add_year))
                                year_print(
                                    "- File was added in {}".format(git_add_year)
                                )
                                year_print(
                                    "- File was modified on {} by {}:\n{}".format(
                                        git_mod_commit.date,
                                        git_mod_commit.author,
                                        git_mod_commit.subject + git_mod_commit.body,
                                    )
                                )
                                year_print(
                                    "=> License should say {}-{}".format(
                                        git_add_year, git_mod_year
                                    )
                                )

                                act_year = year1 if year1 else year2
                                year_print(
                                    err(
                                        "{} But says: {}".format(CROSS_SYMBOL, act_year)
                                    )
                                )

                                if args.fail_year_mismatch:
                                    exit = 1
                                    year_print("\n")
                                else:
                                    year_print("This is not treated as an error\n")
                                valid = False
                            else:
                                if (
                                    int(year1) != git_add_year
                                    or int(year2) != git_mod_year
                                ):
                                    year_print("File: {}".format(src))
                                    year_print(
                                        "Year range {}-{} does not match range from git {}-{}".format(
                                            year1, year2, git_add_year, git_mod_year
                                        )
                                    )
                                    year_print(
                                        "- File was added in {}".format(git_add_year)
                                    )
                                    year_print(
                                        "- File was modified on {} by {}:\n{}".format(
                                            git_mod_commit.date,
                                            git_mod_commit.author,
                                            git_mod_commit.subject
                                            + git_mod_commit.body,
                                        )
                                    )
                                    year_print(
                                        "=> License should say {}-{}".format(
                                            git_add_year, git_mod_year
                                        )
                                    )
                                    year_print(
                                        err(
                                            "{} But says: {}-{}".format(
                                                CROSS_SYMBOL, year1, year2
                                            )
                                        )
                                    )
                                    if args.fail_year_mismatch:
                                        exit = 1
                                        year_print("\n")
                                    else:
                                        year_print("This is not treated as an error\n")
                                    valid = False

                        else:
                            if int(year1) < git_mod_year:
                                year_print("File: {}".format(src))
                                year_print(
                                    "- Year {} does not match git modification year {}".format(
                                        year1, git_mod_year
                                    )
                                )
                                year_print(
                                    "- File was modified on {} by {}:\n{}".format(
                                        git_mod_commit.date,
                                        git_mod_commit.author,
                                        git_mod_commit.subject + git_mod_commit.body,
                                    )
                                )
                                year_print(
                                    "=> License should say {}".format(git_mod_year)
                                )
                                year_print(
                                    err("{} But says: {}".format(CROSS_SYMBOL, year1))
                                )
                                if args.fail_year_mismatch:
                                    exit = 1
                                    year_print("\n")
                                else:
                                    year_print("This is not treated as an error\n")
                                valid = False

                if args.fix and not valid:
                    eprint("-> fixing file (patch year)")
                    year_str = "2016-{}".format(year)
                    if args.check_years:
                        if git_add_year == git_mod_year:
                            year_str = "{}".format(git_add_year)
                        else:
                            year_str = "{}-{}".format(git_add_year, git_mod_year)
                    new_license = raw.format(year=year_str)

                    # preserve rest of file as is
                    if args.check_years:
                        old_license_len = len(license)
                        f.seek(old_license_len)
                        file_body = f.read()
                        f.seek(0)
                        f.truncate()

                    f.seek(0)
                    f.write(new_license)

                    if args.check_years:
                        f.write(file_body)

    print("\n--- INFO ---\n")
    print(info_summary)
    print("\n--- ERROR ---\n")
    print(error_summary)

    if exit != 0 and not args.fix:
        print("License problems found. You can try running again with --fix")

    sys.exit(exit)


if "__main__" == __name__:
    main()
