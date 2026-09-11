#!/bin/python3


# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2023 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0


import argparse
import re
import collections
import pathlib


class InstructionCounter:
    """
    Class for counting the use of certain instructions in translation units.

    The instructions and translation units are counted, not linked. So this
    class cannot reconstruct that a particular instruction was emitted in a
    particular translation unit. But I don't think that's necessary at this
    time."""

    def __init__(self):
        """
        Initialize the counter.

        This creates some empty integer dicts with default value zero."""
        self.instructions = collections.defaultdict(int)
        self.translations = collections.defaultdict(int)

    def add(self, instr, trans):
        """Register the occurance of an instruction in a translation unit."""
        self.instructions[instr] += 1
        self.translations[trans] += 1


def oxford_join(lst):
    """
    Format a list of strings in a human-readable way using an Oxford comma.

    This function takes ["a", "b", "c"] to the string "a, b, and c"."""
    if not lst:
        return ""
    elif len(lst) == 1:
        return str(lst[0])
    elif len(lst) == 2:
        return f"{str(lst[0])} and {str(lst[1])}"
    return f"{', '.join(lst[:-1])}, and {lst[-1]}"


def run(files, source, build):
    """
    Perform a search for FP64 instructions in a list of files.

    This function takes a list of file paths as well as the root path of the
    source code and the build path. These are necessary because the paths
    reported in the PTX emitted by NVCC are relative to the build directory. In
    order to get the paths relative to the GitHub root directory (which GitHub
    Action Commands require), we need to do some path magic."""
    # Create a dictionary of counters. The keys in this dictionary are line
    # information tuples (source file name and line) and the values are counter
    # objects which count how many times that line generates each instruction,
    # and how many times it generates instructions in a given translation unit.
    counter = collections.defaultdict(InstructionCounter)

    # Resolve the source and build paths if they are relevant. Since these are
    # constant, we can move this operation out of the loop.
    source_path = source.resolve()
    build_path = build.resolve()

    # Iterate over the list of files that we are given by the user. We do this
    # multi-file analysis so we can analyse the mapping of shared source code
    # to multiple translation units.
    for n in files:
        # Read the PTX file and split it into multiple lines. This is NOT a
        # proper parsing of PTX and could break, but works for now.
        with open(n, "r") as f:
            lines = f.read().split("\n")

        # At the beginning of the file, the line data is unknown.
        linedata = None

        # Iterate over the source lines in the PTX.
        for l in lines:
            if m := re.match(r"^//(?P<file>[/\w\-. ]*):(?P<line>\d+)", l):
                # If the line of the form "//[filename]:[line] [code]", we
                # parse the file name and line number, then update the line
                # data. Any subsequent instructions will be mapped onto this
                # source line.
                linedata = (m.group("file"), int(m.group("line")))
            elif m := re.match(
                r"^\s*(?P<instruction>(?:[a-z][a-z0-9]*)(?:\.[a-z][a-z0-9]+)*)", l
            ):
                # If the line is of the form "    [instruction] [operands]", we
                # parse the instruction. The operands are irrelevant.
                if "f64" in m.group("instruction"):
                    # PTX has the pleasant property that all instructions
                    # explicitly specify their operand types (Intel x86 syntax
                    # could learn from this), so if "f64" is contained in the
                    # instruction it will be a double-precision operator. We
                    # now proceed to compute the real path of the line that
                    # produced this instruction.
                    real_path = (build_path / linedata[0]).resolve()
                    if linedata is not None:
                        # If the line data is not none, we have a line to link
                        # this instruction to. We compute the relative path of
                        # the source file to the root of the source directory,
                        # and add the result to the counting dictionary.
                        try:
                            counter[
                                (real_path.relative_to(source_path), linedata[1])
                            ].add(m.group("instruction"), n)
                        except ValueError:
                            pass
                    else:
                        # If we do not have line data, we register an FP64
                        # instruction of unknown origin.
                        counter[None].add(m.group("instruction"), n)

    # After we complete our analysis, we print some output to stdout which will
    # be parsed by GitHub Actions. For the syntax, please refer to
    # https://docs.github.com/en/actions/using-workflows/workflow-commands-for-github-actions
    for dt in counter:
        instrs = oxford_join(
            [f"{counter[dt].instructions[i]} × `{i}`" for i in counter[dt].instructions]
        )
        units = oxford_join(
            [f"`{pathlib.Path(f).name}`" for f in counter[dt].translations]
        )
        details = (
            f"Instruction(s) generated are {instrs} in translation unit(s) {units}."
        )

        # Handle the cases where the source line information is unknown and
        # known, respectively.
        if dt is None:
            print(
                f"::warning title=FP64 instructions emitted in unknown locations::{details}"
            )
        else:
            print(
                f"::warning file={dt[0]},line={dt[1]},title=FP64 instructions emitted::{details}"
            )


if __name__ == "__main__":
    # Construct an argument parser, asking the user for a set of files, as well
    # as their source and build directories, in a fashion similar to what CMake
    # does.
    parser = argparse.ArgumentParser(
        description="Find unwanted 64-bit float operations in annotated PTX."
    )

    parser.add_argument("files", type=str, help="PTX file to use", nargs="+")
    parser.add_argument(
        "--source", "-S", type=pathlib.Path, help="source directory", required=True
    )
    parser.add_argument(
        "--build", "-B", type=pathlib.Path, help="build directory", required=True
    )

    args = parser.parse_args()

    # Finally, run the analysis!
    run(args.files, args.source, args.build)
