# SPDX-PackageName = "traccc, a part of the ACTS project"
# SPDX-FileCopyrightText: CERN
# SPDX-License-Identifier: MPL-2.0

import subprocess
import pathlib
import logging
import os

from . import git

log = logging.getLogger("traccc_benchmark")


DETERMINISTIC_ORDER_COMMIT = "7e7f17ccd2e2b0db8971655773b351a365ee1cfc"
DETERMINISTIC_DEFAULT_COMMIT = "80e34dc2195fc3897f1d7d97a71c0701edaa6aa1"
BOOLEAN_FLAG_COMMIT = "380fc78ba63a79ed5c8f19d01d57636aa31cf4fd"
INHOMOGENEOUS_BFIELD_COMMIT = "3654a64d5fe06509e6bf8be332f5aae7b8ff2da9"


def run_profile(
    build_dir: pathlib.Path, data_dir: str, commit, events=1, ncu_wrapper=None
):
    profile_args = [
        "ncu",
        "--import-source",
        "no",
        "--section LaunchStats",
        "--section Occupancy",
        "--metrics gpu__time_duration.sum",
        "-f",
        "-o",
        build_dir / "profile",
        build_dir / "bin" / "traccc_throughput_st_cuda",
        "--input-directory=%s" % data_dir,
        "--digitization-file=geometries/odd/odd-digi-geometric-config.json",
        "--conditions-file=geometries/odd/odd-digi-geometric-config.json",
        "--detector-file=geometries/odd/odd-detray_geometry_detray.json",
        "--grid-file=geometries/odd/odd-detray_surface_grids_detray.json",
        "--input-events=%d" % events,
        "--cold-run-events=0",
        "--processed-events=%d" % events,
    ]

    if ncu_wrapper is not None:
        profile_args = ncu_wrapper.split() + profile_args

    if git.is_parent_of(commit, DETERMINISTIC_DEFAULT_COMMIT):
        log.info(
            "Commit is a child of (or is) %s; deterministic processing is default",
            DETERMINISTIC_DEFAULT_COMMIT[:8],
        )
    elif git.is_parent_of(commit, DETERMINISTIC_ORDER_COMMIT):
        log.info(
            "Commit is a child of (or is) %s; enabling deterministic processing",
            DETERMINISTIC_ORDER_COMMIT[:8],
        )
        profile_args.append("--deterministic")
    else:
        log.info(
            "Commit is not a child of %s; event order is random",
            DETERMINISTIC_ORDER_COMMIT[:8],
        )

    if git.is_parent_of(commit, BOOLEAN_FLAG_COMMIT):
        log.info(
            "Commit is a child of (or is) %s; using explicit boolean flags",
            BOOLEAN_FLAG_COMMIT[:8],
        )
        profile_args.append("--use-acts-geom-source=1")
        profile_args.append("--use-detray-detector=1")
    else:
        log.info(
            "Commit is not a child of %s; using implicit boolean flags",
            BOOLEAN_FLAG_COMMIT[:8],
        )
        profile_args.append("--use-acts-geom-source")
        profile_args.append("--use-detray-detector")

    if git.is_parent_of(commit, INHOMOGENEOUS_BFIELD_COMMIT):
        log.info(
            "Commit is a child of (or is) %s; enabling inhomogeneous magnetic field",
            INHOMOGENEOUS_BFIELD_COMMIT[:8],
        )
        profile_args.append("--read-bfield-from-file")
        profile_args.append("--bfield-file=geometries/odd/odd-bfield.cvf")
    else:
        log.info(
            "Commit is not a child of %s; disabling inhomogeneous magnetic field",
            INHOMOGENEOUS_BFIELD_COMMIT[:8],
        )

    subprocess.run(
        profile_args,
        stdout=subprocess.DEVNULL,
        check=True,
    )
