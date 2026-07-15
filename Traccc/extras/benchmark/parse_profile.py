# SPDX-PackageName = "traccc, a part of the ACTS project"
# SPDX-FileCopyrightText: CERN
# SPDX-License-Identifier: MPL-2.0

import logging
import argparse
import pathlib

from traccc_bench_tools.parse_profile import parse_profile_ncu, parse_profile_csv
from traccc_bench_tools import types


log = logging.getLogger("traccc_parse_profile")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "input",
        type=pathlib.Path,
        help="the input .ncu-rep or .csv file",
    )

    parser.add_argument(
        "output",
        type=pathlib.Path,
        help="the output .csv file",
    )

    parser.add_argument(
        "--num-sm",
        help="number of SMs in the modelled GPU",
        type=int,
        required=True,
        dest="num_sm",
    )

    parser.add_argument(
        "--num-threads-per-sm",
        help="number of thread slots per SM in the modelled GPU",
        type=int,
        required=True,
        dest="num_threads_per_sm",
    )

    parser.add_argument(
        "--event-delimiter-kernel",
        help="name of kernel that delimits events",
        type=str,
        default="ccl_kernel",
        dest="event_delimiter_kernel",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        help="enable verbose output",
        action="store_true",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if (args.verbose or False) else logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    spec = types.GpuSpec(getattr(args, "num_sm"), getattr(args, "num_threads_per_sm"))

    df = None

    if args.input.suffix == ".ncu-rep":
        log.info("Input file is an NSight Compute profile file")
        df = parse_profile_ncu(
            args.input,
            spec,
            event_marker_kernel=getattr(args, "event_delimiter_kernel", None),
        )
    elif args.input.suffix == ".csv":
        log.info("Input file is a CSV file")
        df = parse_profile_csv(
            args.input,
            spec,
            event_marker_kernel=getattr(args, "event_delimiter_kernel", None),
        )
    else:
        log.fatal("Input file has unknown extension %s", args.input.suffix)
        exit(1)

    df.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
