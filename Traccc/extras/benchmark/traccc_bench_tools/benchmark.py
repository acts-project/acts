# SPDX-PackageName = "traccc, a part of the ACTS project"
# SPDX-FileCopyrightText: CERN
# SPDX-License-Identifier: MPL-2.0

import tempfile
import pathlib
import typing
import logging
import time

from . import types, build, profile, parse_profile

log = logging.getLogger("traccc_benchmark")


def run_benchmark(
    source_dir: pathlib.Path,
    data_dir: pathlib.Path,
    commit,
    gpu_spec: types.GpuSpec,
    parallel: int = 1,
    events: int = 1,
    ncu_wrapper: str = None,
    cc: typing.Union[str, tuple[str, str]] = None,
):
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmppath = pathlib.Path(tmpdirname)

        log.info('Created temporary directory "%s"', str(tmppath))

        build_dir = tmppath / "build"

        log.info('Building traccc into build directory "%s"', build_dir)

        log.info("Running configuration step")

        start_time = time.time()

        build.configure(source_dir, build_dir, commit, cc=cc)

        end_time = time.time()

        log.info(
            "Completed configuration step in %.1f seconds",
            end_time - start_time,
        )

        log.info("Running build step with %d thread(s)", parallel)

        start_time = time.time()

        build.build(build_dir, parallel=parallel)

        end_time = time.time()

        log.info("Completed build step in %.1f seconds", end_time - start_time)

        log.info("Running benchmark step")

        start_time = time.time()

        profile.run_profile(
            build_dir,
            data_dir=data_dir,
            commit=commit,
            events=events,
            ncu_wrapper=ncu_wrapper,
        )

        end_time = time.time()

        log.info("Completed benchmark step in %.1f seconds", end_time - start_time)

        log.info("Running data processing step")

        start_time = time.time()

        result_df = parse_profile.parse_profile_ncu(
            build_dir / "profile.ncu-rep",
            gpu_spec,
        )

        end_time = time.time()

        log.info(
            "Completed data processing step in %.1f seconds",
            end_time - start_time,
        )

        return result_df
