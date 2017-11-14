#!/usr/bin/env python3
import argparse
import os
import sys
from subprocess import check_output
import re
import difflib
from datetime import datetime


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("--fix", action="store_true")

    args = p.parse_args()

    if os.path.isdir(args.input):
        srcs = str(check_output(["find", args.input, "-iname", "*.cpp", "-or", "-iname", "*.hpp", "-or", "-iname", "*.ipp"]), "utf-8").strip().split("\n")
        srcs = filter(lambda p: not p.startswith("./build"), srcs)
    else:
        srcs = [args.input]

    year = int(datetime.now().strftime("%Y"))

    raw = """// This file is part of the ACTS project.
//
// Copyright (C) {year} ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/."""

    reg = (
        r"\A// This file is part of the ACTS project.\n"
        +r"//\n"
        +r"// Copyright \(C\) (?P<year>.*) ACTS project team\n"
        +r"//\n"
        +r"// This Source Code Form is subject to the terms of the Mozilla Public\n"
        +r"// License, v\. 2\.0\. If a copy of the MPL was not distributed with this\n"
        +r"// file, You can obtain one at http://mozilla.org/MPL/2.0/.\Z"
        )

    ref = re.compile(reg, re.M)
    clean_re = re.compile(r"(\(C\)) (.*) (ACTS)", re.M)
    year_re = re.compile(r"^(?P<year1>20\d{2}|(?P<year2>20\d{2})-(?P<year3>20\d{2}))$")
    extract_re = re.compile(r"(20\d{2})-?(20\d{2})?")
    
    def clean(s):
        return clean_re.sub(r"\1 XXXX \3", s)
    def get_clean_lines(s):
        return [clean(l)+"\n" for l in s.split("\n")]
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

    exit = 0
    for src in srcs:
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
                sys.stderr.write("Invalid / missing license in "+src+"\n")
                
                exp = [l+"\n" for l in raw.format(year="XXXX").split("\n")]
                act = get_clean_lines(license)

                diff = difflib.unified_diff(exp, act)
                sys.stderr.writelines(diff)
                sys.stderr.write("\n")

                if args.fix:
                    print("-> fixing file (prepend)")
                    f.seek(0)
                    file_content = f.read()
                    f.seek(0)
                    f.write(rawstr+"\n\n")
                    f.write(file_content)

                exit = 1
            else:
                # we have a match, need to verify year string is right
                year_act = m.group("year")
                ym = year_re.match(year_act)
                valid = True
                if not ym:
                    sys.stderr.write("Year string does not match format in {}\n".format(src))
                    sys.stderr.write("Expected: YYYY or YYYY-YYYY (or year or year range)\n")
                    sys.stderr.write("Actual:   {}\n\n".format(year_act))
                    
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
                        sys.stderr.write("Year string is not valid in {}\n".format(src))
                        sys.stderr.write("Year string is: "+year_act+"\n\n")
                        exit = 1
                        valid = False
                    
                if args.fix and not valid:
                    print("-> fixing file (patch year)")
                    year_str = "2016-{}".format(year)
                    new_license = raw.format(year=year_str)
                    # only license year is invalid, license exist, we can overwrite
                    f.seek(0)
                    f.write(new_license)

    if exit != 0 and not args.fix:
        sys.stderr.write("License problems found. You can try running again with --fix\n")
    sys.exit(exit)



if "__main__" == __name__:
    main()
