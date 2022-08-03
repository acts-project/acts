#!/usr/bin/env python3

"""
This script merges two clang-tidy reports files
"""

import argparse
import json

from item import Item, ItemCollection


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "inputfiles", help="The input file containing the warnings", nargs="+"
    )
    p.add_argument(
        "--output",
        required=True,
        help="The resulting JSON file",
    )
    args = p.parse_args()

    items = set()

    for infile in args.inputfiles:
        with open(infile, "r", encoding="utf-8") as f:
            items = items | set(ItemCollection(__root__=json.load(f)).__root__)

    items = list(items)

    print("Write to", args.output)
    with open(args.output, "w+") as jf:
        jf.write(ItemCollection(__root__=items).json(indent=2))


if "__main__" == __name__:
    main()
