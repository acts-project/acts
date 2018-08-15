#!/usr/bin/env python3

"""
This produces structured analysis results:
It takes the JSON input produced by make_report.py
and applies the limits specifiec in --limitfile.
It then prints a summary and exits with 0 if all warning
categories are below their limits, or 1 if any of the
limits is exceeded.
"""

import argparse
import yaml
import os
import json
import sys
from fnmatch import fnmatch
import codecs
from tabulate import tabulate
from operator import itemgetter

if sys.stdout.encoding != 'UTF-8':
    sys.stdout = codecs.getwriter('utf-8')(sys.stdout.buffer, 'strict')
if sys.stderr.encoding != 'UTF-8':
    sys.stderr = codecs.getwriter('utf-8')(sys.stderr.buffer, 'strict')

from codereport import CodeReport, ReportItem



def analysis(limitfile, results, verbose=False, md=False):
    output = ""
    
    class bcolors:
        HEADER = '\033[95m'
        OKBLUE = '\033[94m'
        OKGREEN = '\033[92m'
        WARNING = '\033[93m'
        FAIL = '\033[91m'
        ENDC = '\033[0m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'

    if not md:
        CROSS_SYMBOL = "\u2716"
        CHECK_SYMBOL = "\u2714"
    else:
        CROSS_SYMBOL = "\u274C"
        CHECK_SYMBOL = "\u2705"

    def colored(s, attr):
        if md: return s
        if sys.stdout.isatty():
            return attr+s+bcolors.ENDC
        else:
            return s

    def green(s):
        return colored(s, bcolors.OKGREEN)
    def orange(s):
        return colored(s, bcolors.WARNING)
    def red(s):
        return colored(s, bcolors.FAIL)

    with open(limitfile, "r") as f:
        config = yaml.load(f)


    limits  = {}
    counts = {}
    for key in config["limits"]:
        limit = config["limits"][key]
        limits[key] = int(limit)
        counts[key] = 0

    with open(results, "r") as f:
        items = json.load(f)

    items = [ReportItem(**item) for item in items]

    codes = {}

    for item in items:
        if item.code in config["ignore"]:
            continue
        if not item.code in codes:
            codes[item.code] = 0
        codes[item.code] += 1
        for pattern in limits.keys():

            if fnmatch(item.code, pattern):
                counts[pattern] += 1

    if verbose:
        output += tabulate(reversed(sorted(list(codes.items()), key=itemgetter(1))), headers=("code", "count"), tablefmt="psql")
        output += "\n\n"

    exit = 0
    lines = []
    for key, count in counts.items():
        limit = limits[key]

        key_s = key
        if md:
            key_s = "`%s`"%key_s

        if count < limit:
            line = (CHECK_SYMBOL, key_s, count, limit)
            lines.append(map(green, map(str, line)))
        elif count == limit:
            line = (CHECK_SYMBOL, key_s, count, limit)
            lines.append(map(orange, map(str, line)))
        else:
            line = (CROSS_SYMBOL, key_s, count, limit)
            lines.append(map(red, map(str, line)))
            exit = 1

    result_status = ("Failed") if exit == 1 else ("Accepted")

    if md:
        output += "# Static analysis results: %s\n\n" % result_status
    else:
        output += "Results: %s\n\n" % result_status

    output += tabulate(lines, headers=("ok", "pattern", "count", "limit"), tablefmt="pipe")

    if exit == 1:
        if not md:
            output += red("\n\n=> Failed rules")

    return exit, output


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--limitfile", required=True,
                   help="The input limit file (yml)")
    p.add_argument("--itemfile", required=True,
                   help="The input item file containing the warnings (json)")
    p.add_argument("--verbose", "-v", action="store_true",
                   help="More output")
    p.add_argument("--markdown", "-md", action="store_true",
                   help="Produce MD output instead of terminal ooutput")

    args = p.parse_args()

    assert os.path.exists(args.limitfile)
    assert os.path.exists(args.itemfile)

    exit, string = analysis(args.limitfile, args.itemfile, args.verbose, args.markdown)
    print(string)
    sys.exit(exit)

if "__main__" == __name__:
    main()

