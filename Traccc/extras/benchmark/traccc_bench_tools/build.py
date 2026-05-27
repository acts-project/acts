# SPDX-PackageName = "traccc, a part of the ACTS project"
# SPDX-FileCopyrightText: CERN
# SPDX-License-Identifier: MPL-2.0

import pathlib
import logging
import subprocess
import typing
import os

from . import git


log = logging.getLogger("traccc_benchmark")


SPACK_LIBS_COMMIT = "069cc80b845c16bf36430fdc90130f0306b47f3e"
COVFIE_V0_15_1_HASH_MISMATCH_BEGIN = "2757ac69fc61bf21fb7959d6fb30d003d6128e44"
COVFIE_V0_15_1_HASH_MISMATCH_END = "00f027941a40157bed2c3f94beecd82df2b34544"


def configure(
    source_dir: pathlib.Path, build_dir: pathlib.Path, commit, cc: str = None
):
    config_args = [
        "cmake",
        "-S",
        source_dir,
        "-B",
        build_dir,
        "-DTRACCC_BUILD_CUDA=ON",
        "-DCMAKE_BUILD_TYPE=Release",
        "-DTRACCC_USE_ROOT=OFF",
    ]

    if cc is not None:
        config_args.append("-DCMAKE_CUDA_ARCHITECTURES=%s" % cc)

    if git.is_parent_of(commit, SPACK_LIBS_COMMIT):
        log.info(
            "Commit is a child of (or is) %s; enabling Spack libraries",
            SPACK_LIBS_COMMIT[:8],
        )
        config_args.append("-DTRACCC_USE_SPACK_LIBS=ON")
    else:
        log.info(
            "Commit is not a child of %s; disabling Spack libraries",
            SPACK_LIBS_COMMIT[:8],
        )
        config_args.append("-DTRACCC_USE_SYSTEM_ACTS=ON")
        config_args.append("-DTRACCC_USE_SYSTEM_TBB=ON")

    if git.is_parent_of(
        commit, COVFIE_V0_15_1_HASH_MISMATCH_BEGIN
    ) and not git.is_parent_of(commit, COVFIE_V0_15_1_HASH_MISMATCH_END):
        log.info(
            "Commit is a child of (or is) %s but not %s; removing covfie v0.15.1 hash",
            COVFIE_V0_15_1_HASH_MISMATCH_BEGIN[:8],
            COVFIE_V0_15_1_HASH_MISMATCH_END[:8],
        )
        config_args.append(
            "-DTRACCC_COVFIE_SOURCE='URL;https://github.com/acts-project/covfie/archive/refs/tags/v0.15.1.tar.gz'"
        )

    subprocess.run(
        config_args,
        check=True,
        stdout=subprocess.DEVNULL,
    )


def build(build_dir: pathlib.Path, parallel: int = 1):
    subprocess.run(
        [
            "cmake",
            "--build",
            build_dir,
            "--",
            "-j",
            str(parallel),
            "traccc_throughput_st_cuda",
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
