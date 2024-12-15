#!/usr/bin/env python3

import os
import argparse
import json
import urllib.request
import urllib.error
import re
import subprocess
from typing import Tuple, Dict, Optional


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tag", type=str, required=True, help="Tag to use")
    parser.add_argument("--arch", type=str, required=True, help="Architecture to use")
    parser.add_argument(
        "--compiler-binary",
        type=str,
        default=os.environ.get("CXX"),
        help="Compiler to use (defaults to CXX environment variable if set)",
    )
    parser.add_argument(
        "--compiler",
        type=str,
        default=None,
        help="Compiler to use (defaults to compiler binary if set)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output file to write lockfile to",
    )
    args = parser.parse_args()

    print("Fetching lockfiles for tag:", args.tag)
    print("Architecture:", args.arch)

    base_url = f"https://api.github.com/repos/acts-project/ci-dependencies/releases/tags/{args.tag}"

    data = json.loads(fetch_github(base_url))

    lockfiles = parse_assets(data)

    print("Available lockfiles:")
    for arch, compilers in lockfiles.items():
        print(f"> {arch}:")
        for c, (n, _) in compilers.items():
            print(f"  - {c}: {n}")

    if args.arch not in lockfiles:
        print(f"No lockfile found for architecture {args.arch}")
        exit(1)

    if args.compiler_binary is not None:
        compiler = determine_compiler_version(args.compiler_binary)
        print("Compiler:", args.compiler_binary, f"{compiler}")
    elif args.compiler is not None:
        if not re.match(r"^([\w-]+)@(\d+\.\d+\.\d+)$", args.compiler):
            print(f"Invalid compiler format: {args.compiler}")
            exit(1)
        compiler = args.compiler
        print("Compiler:", f"{compiler}")
    else:
        compiler = None

    lockfile = select_lockfile(lockfiles, args.arch, compiler)

    print("Selected lockfile:", lockfile)

    if args.output:
        with open(args.output, "wb") as f:
            f.write(fetch_github(lockfile))


def parse_assets(data: Dict) -> Dict[str, Dict[str, Tuple[str, str]]]:
    lockfiles: Dict[str, Dict[str, Tuple[str, str]]] = {}

    for asset in data["assets"]:
        url = asset["browser_download_url"]

        name = asset["name"]
        if not name.endswith(".lock") or not name.startswith("spack_"):
            continue

        m = re.match(r"spack_(.*(?:aarch64|x86_64))(?:_(.*))?\.lock", name)
        if m is None:
            continue

        arch, compiler = m.groups()
        compiler = compiler if compiler else "default"
        lockfiles.setdefault(arch, {})[compiler] = (name, url)

    return lockfiles


def select_lockfile(
    lockfiles: Dict[str, Dict[str, Tuple[str, str]]], arch: str, compiler: Optional[str]
):
    # Default to the default lockfile
    _, lockfile = lockfiles[arch]["default"]

    if compiler is None:
        return lockfile

    # Extract compiler family and version
    compiler_family = compiler.split("@")[0]

    # Find all matching compiler families
    matching_compilers = {
        comp: ver
        for comp, ver in lockfiles[arch].items()
        if comp != "default" and comp.split("@")[0] == compiler_family
    }

    if matching_compilers:
        if compiler in matching_compilers:
            # Exact match found
            _, lockfile = matching_compilers[compiler]
        else:
            # Find highest version of same compiler family
            highest_version = max(
                matching_compilers.keys(),
                key=lambda x: [int(v) for v in x.split("@")[1].split(".")],
            )
            _, lockfile = matching_compilers[highest_version]

    return lockfile


def determine_compiler_version(binary: str):
    try:
        result = subprocess.run([binary, "--version"], capture_output=True, text=True)

        line = result.stdout.split("\n", 1)[0]
        print(line)
        if "clang" in line:
            compiler = "clang"
            if "Apple" in line:
                compiler = "apple-clang"
        elif "gcc" in line or "GCC" in line or "g++" in line:
            compiler = "gcc"
        else:
            print(f"Unknown compiler: {binary}")
            exit(1)

        m = re.search(r"(\d+\.\d+\.\d+)", line)
        if m is None:
            print(f"Failed to determine version for compiler: {binary}")
            exit(1)
        (version,) = m.groups()
        return f"{compiler}@{version}"

    except (subprocess.SubprocessError, FileNotFoundError):
        print(f"Failed to determine version for compiler: {binary}")
        exit(1)


def fetch_github(base_url):
    headers = {}
    if "GITHUB_TOKEN" in os.environ:
        headers["Authorization"] = f"token {os.environ['GITHUB_TOKEN']}"

    try:
        req = urllib.request.Request(base_url, headers=headers)
        with urllib.request.urlopen(req) as response:
            return response.read()
    except urllib.error.URLError as e:
        print(f"Failed to fetch from {base_url}: {e}")
        exit(1)
    except json.JSONDecodeError as e:
        print(f"Failed to parse JSON response: {e}")
        exit(1)


if __name__ == "__main__":
    main()
