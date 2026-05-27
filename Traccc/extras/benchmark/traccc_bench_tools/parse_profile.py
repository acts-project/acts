# SPDX-PackageName = "traccc, a part of the ACTS project"
# SPDX-FileCopyrightText: CERN
# SPDX-License-Identifier: MPL-2.0

import pathlib
import pandas
import functools
import numpy
import tempfile
import subprocess
import operator
import logging

from .types import GpuSpec
from .kernels import simplify_name, map_name
from .utils import harmonic_sum


log = logging.getLogger("traccc_bench_tools.parse_profile")


def parse_triple(triple):
    assert triple[0] == "(" and triple[-1] == ")"
    x, y, z = triple[1:-1].split(", ")
    return int(x), int(y), int(z)


def parse_profile_ncu(file: pathlib.Path, gpu_spec: GpuSpec, **kwargs):
    with tempfile.TemporaryDirectory() as tmpdirname:
        tmppath = pathlib.Path(tmpdirname)

        profile_file = tmppath / "profile.csv"

        with open(profile_file, "w") as f:
            subprocess.run(
                [
                    "ncu",
                    "-i",
                    file,
                    "--csv",
                    "--print-units",
                    "base",
                ],
                stdout=f,
            )

        return parse_profile_csv(profile_file, gpu_spec, **kwargs)


def parse_profile_csv(file: pathlib.Path, gpu_spec: GpuSpec, event_marker_kernel=None):
    if event_marker_kernel is None:
        event_marker_kernel = "ccl_kernel"

    df = pandas.read_csv(file)

    ndf = df[df["Metric Name"] == "gpu__time_duration.sum"][
        ["ID", "Kernel Name", "Block Size", "Grid Size", "Metric Value", "Metric Unit"]
    ]

    assert (ndf["Metric Unit"] == "ns").all()

    ndf["ThreadsPerBlock"] = df["Block Size"].apply(
        lambda x: functools.reduce(operator.mul, (parse_triple(x)))
    )
    ndf["BlocksPerGrid"] = df["Grid Size"].apply(
        lambda x: functools.reduce(operator.mul, (parse_triple(x)))
    )
    ndf["TotalThreads"] = ndf["ThreadsPerBlock"] * ndf["BlocksPerGrid"]
    ndf = ndf.drop(
        ["Block Size", "Grid Size", "ThreadsPerBlock", "BlocksPerGrid", "Metric Unit"],
        axis=1,
    )
    ndf["Metric Value"] = ndf["Metric Value"].apply(lambda x: int(x.replace(",", "")))
    ndf["Kernel Name"] = ndf["Kernel Name"].apply(simplify_name)
    ndf["Kernel Name"] = ndf["Kernel Name"].apply(map_name)

    curr_evt_id = None
    evt_ids = []
    for x in ndf.iloc:
        if x["Kernel Name"] == event_marker_kernel:
            if curr_evt_id is None:
                curr_evt_id = 0
            else:
                curr_evt_id += 1
        evt_ids.append(curr_evt_id)

    num_unidentified_event_ids = sum(1 for x in evt_ids if x is None)

    if num_unidentified_event_ids > 0:
        log.warning(
            "Event ID could not be determined for %d kernels",
            num_unidentified_event_ids,
        )

    ndf["EventID"] = evt_ids

    thr_occ = df[df["Metric Name"] == "Theoretical Occupancy"]

    ndf = ndf.merge(
        thr_occ[["ID", "Metric Value"]],
        on="ID",
        how="left",
        validate="one_to_one",
        suffixes=("", "R"),
    )

    ndf["Occupancy"] = ndf["Metric ValueR"].apply(lambda x: float(x) / 100.0)

    ndf["k"] = ndf["TotalThreads"] / (
        gpu_spec.n_sm * gpu_spec.n_threads_per_sm * ndf["Occupancy"]
    )

    ndf["Latency"] = ndf["Metric Value"] / 1e9

    ndf["Throughput"] = (numpy.ceil(ndf["k"]) / ndf["k"]) / (ndf["Metric Value"] / 1e9)
    ndf["RecThroughput"] = 1.0 / ndf["Throughput"]

    ndf = ndf.drop(["Metric ValueR"], axis=1)

    ndf = ndf.groupby(["Kernel Name", "EventID"], as_index=False).agg(
        {"Throughput": harmonic_sum, "RecThroughput": "sum", "Latency": "sum"},
    )

    ndf = ndf.groupby("Kernel Name", as_index=False).agg(
        ThroughputMean=("Throughput", "mean"),
        ThroughputStd=("Throughput", "std"),
        RecThroughputMean=("RecThroughput", "mean"),
        RecThroughputStd=("RecThroughput", "std"),
        LatencyMean=("Latency", "mean"),
        LatencyStd=("Latency", "std"),
    )

    return ndf.fillna(0)
