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


log = logging.getLogger("traccc_benchmark")


KNOWN_BROKEN_COMMITS = [
    "9444f505d62d1213da7b3a502da3a233d524d264",
]


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
        "-f",
        "--from",
        type=str,
        help="the first commit in range (exclusive)",
        required=True,
    )

    parser.add_argument(
        "-t",
        "--to",
        type=str,
        help="the last commit in range (inclusive)",
        default="HEAD",
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

    commit_range = repo.iter_commits(
        "%s...%s" % (getattr(args, "from"), getattr(args, "to"))
    )

    commit_list = list(commit_range)
    commit_str_list = list(str(x) for x in commit_list)

    log.info(
        "Examining a total of %d commits between %s and %s",
        len(commit_list),
        getattr(args, "from"),
        getattr(args, "to"),
    )

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

    known_commits = set(x["commit"] for x in results)

    log.info("Currently have pre-existing results for %d commits", len(known_commits))

    commits_to_skip = set()

    for progress_id, x in enumerate(commit_list):
        log.info(
            "Considering commit %s (%u out of %u)", x, progress_id, len(commit_list)
        )

        parents = x.parents

        # If this commit has a parent with a single parent, and that
        # grandparent is also a parent of this current commit, we have a
        # single-commit merge which we don't need to measure separately.
        if len(parents) == 2:
            for p1 in parents:
                if len(p1.parents) == 1 and p1.parents[0] in parents:
                    log.info("Commit %s is a triangle commit; adding to skip list", p1)
                    commits_to_skip.add(str(p1))

        # Skip commits which we have already benchmarked
        if str(x) in known_commits:
            log.info("Commit %s is already know; skipping", x)
            continue

        if str(x) in commits_to_skip:
            log.info("Commit %s is marked for skipping", x)
            continue

        if str(x) in KNOWN_BROKEN_COMMITS:
            log.info("Commit %s is known to be broken; skipping", x)
            continue

        try:
            log.info("Running benchmark for commit %s", x)

            log.info("Running checkout step")

            start_time = time.time()
            repo.git.checkout(x)
            end_time = time.time()

            log.info("Completed checkout step in %.1f seconds", end_time - start_time)

            result_df = benchmark.run_benchmark(
                source_dir=args.repo,
                data_dir=args.data,
                commit=x,
                gpu_spec=types.GpuSpec(
                    getattr(args, "num_sm"), getattr(args, "num_threads_per_sm")
                ),
                parallel=args.parallel,
                events=args.events,
                ncu_wrapper=getattr(args, "ncu_wrapper", None),
                cc=getattr(args, "cc", None),
            )

            for y in result_df.iloc:
                results.append(
                    {
                        "commit": str(x),
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
    output_results = sorted(
        [x for x in results if x["commit"] in commit_str_list],
        key=lambda x: -commit_str_list.index(x["commit"]),
    )
    log.info("Keeping a total of %d results after pruning", len(output_results))

    log.info("Writing data to %s", args.db)
    with open(args.db, "w") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "commit",
                "kernel",
                "throughput",
                "throughput_dev",
                "rec_throughput",
                "rec_throughput_dev",
            ],
        )
        writer.writeheader()

        for i in output_results:
            writer.writerow(i)

    log.info("Checking out repository to previous commit %s", old_commit_hash)
    repo.git.checkout(old_commit_hash)
    log.info("Processing complete; goodbye!")


if __name__ == "__main__":
    main()
