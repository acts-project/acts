#!/usr/bin/env python3
import hashlib
from pathlib import Path
import sys
from typing import Optional
import argparse

import uproot
import numpy as np
import awkward as ak


def hash_root_file(
    path: Path, tree_name: Optional[str] = None, ordering_invariant: bool = True
) -> str:
    rf = uproot.open(path)

    gh = hashlib.sha256()

    for tree_name in sorted(rf.keys()):
        gh.update(tree_name.encode("utf8"))

        tree = rf[tree_name]
        keys = list(sorted(tree.keys()))

        branches = tree.arrays(library="ak")

        if not ordering_invariant:

            h = hashlib.sha256()
            for name in keys:
                h.update(name.encode("utf8"))
                arr = branches[name]
                arr = ak.flatten(arr, axis=None)
                arr = np.array(arr)
                h.update(arr.tobytes())
            gh.update(h.digest())

        else:
            items = np.array([])

            for entry in branches:
                h = hashlib.md5()
                for col in keys:
                    try:
                        h.update(np.array(entry[col]).tobytes())
                    except ValueError:
                        h.update(("%s" % entry[col]).encode("utf8"))
                items = np.append(items, h.digest())

            items.sort()

            h = hashlib.sha256()
            h.update("".join(keys).encode("utf8"))
            h.update(items.tobytes())

            gh.update(h.digest())
    return gh.hexdigest()


def assert_root_hash(file: Path, refhash: str):
    __tracebackhide__ = True
    act_hash = hash_root_file(file)
    assert refhash == act_hash, f"{refhash} != {act_hash}"


if "__main__" == __name__:
    p = argparse.ArgumentParser(
        description="Calculate a hash of the numeric content of a root file"
    )

    p.add_argument(
        "input_file", type=Path, help="The input ROOT file to calculate a hash for"
    )
    p.add_argument(
        "--tree-name",
        help="Explicit tree name to use. This needs to be specified if the ROOT file contains more than one tree.",
    )
    p.add_argument(
        "--no-ordering-invariant",
        "-n",
        action="store_true",
        help="Calculate a hash that is not invariant under reordering of entries? (faster than invariant)",
    )

    args = p.parse_args()

    print(
        hash_root_file(
            path=args.input_file,
            tree_name=args.tree_name,
            ordering_invariant=not args.no_ordering_invariant,
        )
    )
