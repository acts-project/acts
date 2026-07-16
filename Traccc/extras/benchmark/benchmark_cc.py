# SPDX-PackageName = "traccc, a part of the ACTS project"
# SPDX-FileCopyrightText: CERN
# SPDX-License-Identifier: MPL-2.0

import argparse
import sys
import csv
import git
import pathlib
import shutil
import logging
import tempfile
import subprocess
import os
import time

from traccc_bench_tools import parse_profile, types, build, profile, benchmark


ALL_COMPUTE_CAPABILITIES = [
    "50",
    "52",
    "60",
    "61",
    "70",
    "75",
    "80",
    "86",
    "89",
    "90",
    "100",
    "120",
]


log = logging.getLogger("traccc_benchmark")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "repo",
        type=pathlib.Path,
        help="the traccc build repository",
    )

    parser.add_argument(
        "db",
        type=pathlib.Path,
        help="the CSV database",
    )

    parser.add_argument(
        "data",
        type=pathlib.Path,
        help="the traccc dataset to use",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        help="enable verbose output",
        action="store_true",
    )

    parser.add_argument(
        "-j",
        "--parallel",
        help="number of threads to use for compilation",
        default=1,
        type=int,
    )

    parser.add_argument(
        "-e",
        "--events",
        help="number of events to process per commit",
        type=int,
        default=10,
    )

    parser.add_argument(
        "--cc",
        help="Compute Capability of the modelled GPU",
        type=str,
        required=True,
        dest="cc",
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
        "--ncu-wrapper",
        help="wrapper to use around the ncu command",
        type=str,
        dest="ncu_wrapper",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if (args.verbose or False) else logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    log.info(
        "Using GPU with %d SMs and %d thread slots per SM (%d thread slots total)",
        getattr(args, "num_sm"),
        getattr(args, "num_threads_per_sm"),
        getattr(args, "num_sm") * getattr(args, "num_threads_per_sm"),
    )

    repo = git.Repo(args.repo)

    log.info("Using git repository at %s", repo.git_dir)

    if repo.is_dirty():
        e = "Repository is dirty; please clean it before use!"
        log.fatal(e)
        raise RuntimeError(e)

    if "TRACCC_TEST_DATA_DIR" not in os.environ:
        e = 'Environment variable "TRACCC_TEST_DATA_DIR" is not set; aborting!'
        log.fatal(e)
        raise RuntimeError(e)

    for exec in ["ncu", "g++", "nvcc", "cmake"]:
        if shutil.which(exec) is None:
            e = 'Executable "%s" is not available; aborting' % exec
            log.fatal(e)
            raise RuntimeError(e)

    old_commit_hash = repo.head.object.hexsha
    log.info("Current commit hash is %s", old_commit_hash)

    if args.db.is_file():
        log.info('Database file "%s" already exists; creating a backup', args.db)
        shutil.copy(str(args.db), str(args.db) + ".bak")
        results = []
        with open(args.db, "r") as f:
            reader = csv.DictReader(f)
            for i in reader:
                results.append(i)
        log.info("Database contained %d pre-existing results", len(results))
    else:
        log.info('Database file "%s" does not exist; starting from scratch', args.db)
        results = []

    known_ccs = set(x["cc"] for x in results)

    log.info(
        "Currently have pre-existing results for %d Compute Capabilities",
        len(known_ccs),
    )

    for progress_id, x in enumerate(ALL_COMPUTE_CAPABILITIES):
        log.info(
            "Considering Compute Capability %s.%s (%u out of %u)",
            x[:-1],
            x[-1],
            progress_id,
            len(ALL_COMPUTE_CAPABILITIES),
        )

        # Skip commits which we have already benchmarked
        if x in known_ccs:
            log.info(
                "Compute capability %s.%s is already know; skipping", x[:-1], x[-1]
            )
            continue

        if int(x) > int(args.cc):
            log.info(
                "Compute Capability %s.%s is newer than target %s.%s; skipping",
                x[:-1],
                x[-1],
                args.cc[:-1],
                args.cc[-1],
            )
            continue

        try:
            log.info("Running benchmark for Compute Capability %s.%s", x[:-1], x[-1])

            result_df = benchmark.run_benchmark(
                source_dir=args.repo,
                data_dir=args.data,
                commit=repo.head.object,
                gpu_spec=types.GpuSpec(
                    getattr(args, "num_sm"), getattr(args, "num_threads_per_sm")
                ),
                parallel=args.parallel,
                events=args.events,
                ncu_wrapper=getattr(args, "ncu_wrapper", None),
                cc=x,
            )

            for y in result_df.iloc:
                results.append(
                    {
                        "cc": str(x),
                        "kernel": y["Kernel Name"],
                        "throughput": y["ThroughputMean"],
                        "throughput_dev": y["ThroughputStd"],
                        "rec_throughput": y["RecThroughputMean"],
                        "rec_throughput_dev": y["RecThroughputStd"],
                    }
                )

        except Exception as e:
            log.exception(e)
        except KeyboardInterrupt as e:
            log.info("Received keyboard interrupt; skipping to post-processing")
            break

    log.info("Gathered a total of %d results (incl. pre-existing)", len(results))

    log.info("Writing data to %s", args.db)
    with open(args.db, "w") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "cc",
                "kernel",
                "throughput",
                "throughput_dev",
                "rec_throughput",
                "rec_throughput_dev",
            ],
        )
        writer.writeheader()

        for i in results:
            writer.writerow(i)

    log.info("Processing complete; goodbye!")


if __name__ == "__main__":
    main()
