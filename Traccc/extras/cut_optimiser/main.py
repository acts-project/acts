import argparse
import json
import csv
import pathlib
import logging
import functools
import operator
import os
import shutil
import itertools
import random
import tempfile
import time
import subprocess
import re


import traccc_bench_tools.parse_profile
import traccc_bench_tools.types


log = logging.getLogger("traccc_cut_optimiser")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "exec",
        type=pathlib.Path,
        help="the traccc ",
    )

    parser.add_argument(
        "db",
        type=pathlib.Path,
        help="the CSV database file",
    )

    parser.add_argument(
        "parameters",
        type=pathlib.Path,
        help="the JSON parameter space file",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        help="enable verbose output",
        action="store_true",
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

    parser.add_argument(
        "--timeout",
        help="timeout in seconds for benchmarks",
        default=120,
        type=int,
    )

    parser.add_argument(
        "--random",
        help="enable random search rather than grid search",
        action="store_true",
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

    gpu_spec = traccc_bench_tools.types.GpuSpec(
        getattr(args, "num_sm"), getattr(args, "num_threads_per_sm")
    )

    with open(args.parameters, "r") as f:
        params = json.load(f)

    parameter_space = params["parameters"]
    parameter_names = sorted(parameter_space)

    log.info(
        "Running optimisation for %d parameters: %s",
        len(params["parameters"]),
        ", ".join(parameter_names),
    )

    csv_field_names = parameter_names + [
        "success",
        "rec_throughput",
        "efficiency",
        "fake_rate",
        "duplicate_rate",
    ]

    if args.db.is_file():
        log.info('Database file "%s" already exists; creating a backup', args.db)
        shutil.copy(str(args.db), str(args.db) + ".bak")
        results = {}
        with open(args.db, "r") as f:
            reader = csv.DictReader(f)
            for i in reader:
                key = (i[k] for k in parameter_names)
                if set(i.keys()) != set(csv_field_names):
                    raise ValueError(
                        "Input database has unexpected columns or is missing columns"
                    )
                if i["success"] != 0:
                    results[key] = {
                        "rec_throughput": i["rec_throughput"],
                        "efficiency": i["efficiency"],
                        "fake_rate": i["fake_rate"],
                        "duplicate_rate": i["duplicate_rate"],
                    }
                else:
                    results[key] = None
        log.info("Database contained %d pre-existing results", len(results))
    else:
        log.info('Database file "%s" does not exist; starting from scratch', args.db)
        results = {}

    log.info("Starting optimisation with %d known results", len(results))

    log.info(
        "Total parameter space has size %d",
        functools.reduce(operator.mul, [len(x) for x in parameter_space.values()], 1),
    )

    if "TRACCC_TEST_DATA_DIR" not in os.environ:
        e = 'Environment variable "TRACCC_TEST_DATA_DIR" is not set; aborting!'
        log.fatal(e)
        raise RuntimeError(e)

    for exec in ["ncu"]:
        if shutil.which(exec) is None:
            e = 'Executable "%s" is not available; aborting' % exec
            log.fatal(e)
            raise RuntimeError(e)

    space_product = itertools.product(
        *[parameter_space[name] for name in parameter_names]
    )
    random_retry_count = 0
    it_count = 0

    while True:
        if args.random:
            param_dict = {n: random.choice(vs) for n, vs in parameter_space.items()}
            param_list = tuple(param_dict.values())

            if param_list in results:
                random_retry_count += 1
                if random_retry_count < 10:
                    log.info(
                        "Configuration is already known; continuing with random abort counter %d",
                        random_retry_count,
                    )
                    continue
                else:
                    log.info(
                        "Configuration is already known and max retry counter is reached; aborting"
                    )
                    break
            else:
                random_retry_count = 0
        else:
            try:
                param_vector = next(space_product)
                assert len(param_vector) == len(parameter_names)
                param_dict = {parameter_names[i]: p for i, p in enumerate(param_vector)}
                param_list = tuple(param_dict.values())
            except StopIteration:
                log.info("Design space has been exhausted; exiting")
                break
            if param_list in results:
                log.info(
                    "Configuration %s is already known; continuing", str(param_dict)
                )
                continue

        try:
            log.info("Running benchmark for parameters %s", str(param_dict))

            with tempfile.TemporaryDirectory() as tmpdirname:
                tmppath = pathlib.Path(tmpdirname)

                log.info('Created temporary directory "%s"', str(tmppath))

                log.info("Running benchmark step")

                start_time = time.time()

                profile_args = [
                    "ncu",
                    "--import-source",
                    "no",
                    "--section",
                    "LaunchStats",
                    "--section",
                    "Occupancy",
                    "--metrics",
                    "gpu__time_duration.sum",
                    "-f",
                    "-o",
                    "profile",
                    str(args.exec.resolve()),
                    "--input-directory=%s" % params["input"]["event_dir"],
                    "--digitization-file=%s" % params["input"]["digitization_file"],
                    "--conditions-file=%s" % params["input"]["conditions_file"],
                    "--detector-file=%s" % params["input"]["detector_file"],
                    "--grid-file=%s" % params["input"]["grid_file"],
                    "--material-file=%s" % params["input"]["material_file"],
                    "--input-events=1",
                    "--check-performance",
                    "--use-acts-geom-source=on",
                ]

                for k, v in params["config"].items():
                    profile_args.append("--%s=%s" % (k, str(v)))

                for k, v in param_dict.items():
                    profile_args.append("--%s=%s" % (k, str(v)))

                if args.ncu_wrapper is not None:
                    profile_args = args.ncu_wrapper.split() + profile_args

                try:
                    result = subprocess.run(
                        profile_args,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        cwd=tmppath,
                        check=True,
                        timeout=args.timeout,
                    )
                except subprocess.CalledProcessError as e:
                    log.warning("Process failed to execute; continuing")
                    results[param_list] = None
                    continue
                except subprocess.TimeoutExpired as e:
                    log.warning("Process timed out; marking as failure and continuing")
                    results[param_list] = None
                    continue

                stdout = result.stdout.decode("utf-8")

                if match := re.search(
                    r"Total track efficiency was (\d+(?:\.\d*)?)%", stdout
                ):
                    result_efficiency = float(match.group(1)) / 100.0
                else:
                    raise ValueError("Efficiency could not be parsed from stdout!")

                if match := re.search(
                    r"Total track duplicate rate was (\d+(?:\.\d*)?)", stdout
                ):
                    result_duplicate_rate = float(match.group(1))
                else:
                    raise ValueError("Dupliacate rate could not be parsed from stdout!")

                if match := re.search(
                    r"Total track fake rate was (\d+(?:\.\d*)?)", stdout
                ):
                    result_fake_rate = float(match.group(1))
                else:
                    raise ValueError("Fake rate could not be parsed from stdout!")

                log.info(
                    "Physics performance was %.1f%% efficiency, %.1f fake rate, %.1f duplicate rate",
                    result_efficiency * 100.0,
                    result_duplicate_rate,
                    result_fake_rate,
                )

                end_time = time.time()

                log.info(
                    "Completed benchmark step in %.1f seconds", end_time - start_time
                )

                log.info("Running profile processing step")

                start_time = time.time()

                result_df = traccc_bench_tools.parse_profile.parse_profile_ncu(
                    tmppath / "profile.ncu-rep",
                    gpu_spec,
                    event_marker_kernel="count_grid_capacities",
                )

                total_rec_throughput = result_df["RecThroughputMean"].sum()

                log.info(
                    "Compute performance was a combined reciprocal throughput of %.1fms",
                    1000.0 * total_rec_throughput,
                )

                end_time = time.time()

                log.info(
                    "Completed profile processing step in %.1f seconds",
                    end_time - start_time,
                )

                results[param_list] = {
                    "rec_throughput": total_rec_throughput,
                    "efficiency": result_efficiency,
                    "fake_rate": result_fake_rate,
                    "duplicate_rate": result_duplicate_rate,
                }

        except Exception as e:
            log.exception(e)
            results[param_list] = None
        except KeyboardInterrupt as e:
            log.info("Received keyboard interrupt; skipping to post-processing")
            break

    log.info("Gathered a total of %d results (incl. pre-existing)", len(results))

    log.info("Writing data to %s", args.db)
    with open(args.db, "w") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=csv_field_names,
        )
        writer.writeheader()

        for k, v in results.items():

            param_dict = {parameter_names[i]: p for i, p in enumerate(k)}
            if v is None:
                writer.writerow(
                    dict(
                        **param_dict,
                        success=0,
                        rec_throughput=0.0,
                        efficiency=0.0,
                        fake_rate=0.0,
                        duplicate_rate=0.0,
                    )
                )
            else:
                writer.writerow(
                    dict(
                        **param_dict,
                        **v,
                        success=1,
                    )
                )

    log.info("Processing complete; goodbye!")


if __name__ == "__main__":
    main()
