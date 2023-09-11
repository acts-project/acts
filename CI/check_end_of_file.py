#!/usr/bin/env python3

import sys
import argparse


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input", nargs="+")
    p.add_argument("--fix", action="store_true")
    p.add_argument("--reject-multiple-newlines", action="store_true")

    args = p.parse_args()

    failed = []

    for filename in args.input:
        with open(filename) as f:
            lines = f.readlines()

        if not lines[-1].endswith("\n"):
            print(f"Missing newline at end of file: {filename}")
            if args.fix:
                with open(filename, "a") as f:
                    f.write("\n")
            else:
                failed.append(filename)
        elif args.reject_multiple_newlines and lines[-1] == "\n":
            print(f"Multiple newlines at end of file: {filename}")
            if args.fix:
                while lines[-1] == "\n":
                    lines.pop(-1)
                with open(filename, "w") as f:
                    f.write("".join(lines))
            else:
                failed.append(filename)

    if failed:
        print(f"failed for files: {' '.join(failed)}")
        return 1

    print("success")
    return 0


if "__main__" == __name__:
    sys.exit(main())
