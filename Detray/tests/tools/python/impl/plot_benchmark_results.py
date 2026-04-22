# Detray library, part of the ACTS project (R&D line)
#
# (c) 2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# detray includes
import plotting

# python includes
from collections import namedtuple
import json
import itertools
import math
import numpy as np
import os
import pandas as pd
import sys

# Plot types for benchmarks
benchmark_plots = namedtuple(
    "benchmark_plots",
    "latency throughput, weak_scaling, strong_scaling",
    defaults=[None, None, None, None],
)

# How to label plots
label_data = namedtuple("label_data", "title label x_axis y_axis")

# Define labels for google benchmark data collections
label_dict = {
    "latency": label_data(
        "Propagation Latency", "", "No. tracks", r"t $[\mathrm{ms}]$"
    ),
    "throughput": label_data(
        "Propagation Throughout", "", "No. tracks", r"Prop. rate $[\mathrm{MHz}]$"
    ),
    "weak_scaling": label_data(
        "Propagation Weak Scaling", "", "No. threads", "Efficiency"
    ),
    "strong_scaling": label_data(
        "Propagation Strong Scaling", "", "No. threads", "Speedup"
    ),
}

# Common options
ldg_loc = "upper left"


""" Read google benchmark data from json file """


def read_benchmark_data(logging, input_path, benchmark_file):

    file_path = os.path.join(input_path, benchmark_file)
    with open(file_path, "r") as file:
        logging.debug(f"Reading file '{file_path}'")

        results = json.load(file)

        context = results["context"]
        data = pd.DataFrame(results["benchmarks"])

        return context, data

    logging.error(f"Could not find file: {benchmark_file}")

    return None, None


""" Adds a column 'x' to the data frame that contains the number of tracks """


def add_track_multiplicity_column(df):

    assert "_TRACKS" in str(
        df["run_name"][0]
    ), "Benchmark case name not correctly formatted: (BM_PROPAGATION_<detector name>_<#tracks>_TRACKS)"

    # The number of tracks is the second last part of the benchmark name
    find_track_multiplicity = lambda n: (int(n.split("_")[-3]))

    # Add new column based on benchmark case name
    df["x"] = df["run_name"].apply(find_track_multiplicity)


""" Read the benchmark data and prepare it for plotting """


def prepare_benchmark_data(logging, input_dir, file):

    # Convert benchmark timings to 'ms'
    unit_conversion = {"ns": 10**-6, "um": 10**-3, "ms": 1, "s": 10**3}

    # Read the data part into a pandas frame
    context, data = read_benchmark_data(logging, input_dir, file)

    if context is None or data is None:
        logging.warning(f"Failed to read data in file: {file}")
        sys.exit(1)

    # Add the number of tracks per benchmark case as new column 'x'
    # A column called 'x' is expected by the 'plot_benchmark' method
    add_track_multiplicity_column(data)

    # Convert timings to 'ms'
    bench_time_unit = data["time_unit"][0]
    to_milliseconds = lambda x: (x * unit_conversion[bench_time_unit])

    data["real_time"] = data["real_time"].apply(to_milliseconds)
    data["cpu_time"] = data["cpu_time"].apply(to_milliseconds)

    # Convert from Hz to MHz
    data["TracksPropagated"] = data["TracksPropagated"] / 1000000

    return context, data


""" Filter the data frame for a specific type of data """


def filter_benchmark_data(df, data_type):

    assert len(df["x"]) != 0, "Data frame has to provide column 'x'"
    assert len(df[data_type]) != 0, f"Data frame has to provide column '{data_type}'"

    # Filter the relevant data from the frame
    median = lambda data_frame: (data_frame["aggregate_name"] == "mean")
    stddev = lambda data_frame: (data_frame["aggregate_name"] == "stddev")

    data, n_tracks = plotting.filter_data(
        data=df, filter=median, variables=[data_type, "x"]
    )

    err = plotting.filter_data(data=df, filter=stddev, variables=[data_type])

    return n_tracks, data, err


"""
Plot the benchmark latency and throughout for different hardware backends and
algebra plugins
"""


def plot_benchmark_case(
    plot_factory,
    x,
    y,
    label,
    y_error=[],
    plot_type="latency",
    title="",
    marker=".",
    plot=None,
    xaxis_format=None,
    yaxis_format=None,
    log_scale=10,
):
    if plot is None:
        # Create new plot
        lgd_ops = plotting.legend_options(
            loc=ldg_loc, horiz_anchor=1.0, vert_anchor=1.02
        )

        labels = label_dict[plot_type]
        x_axis_opts = plotting.axis_options(
            label=labels.x_axis,
            log_scale=log_scale,
            tick_positions=x,
            label_format=xaxis_format,
        )
        y_axis_opts = plotting.axis_options(
            label=labels.y_axis, log_scale=log_scale, label_format=yaxis_format
        )

        # Plot the propagation latency against the number of tracks
        plot_data = plot_factory.graph(
            x=x,
            y=y,
            y_errors=y_error,
            x_axis=x_axis_opts,
            y_axis=y_axis_opts,
            title=title,
            label=label,
            lgd_ops=lgd_ops,
            marker=marker,
            figsize=(18, 8),
        )
    else:
        # Add new data to exiting plot
        plot_data = plot_factory.add_graph(
            plot=plot,
            x=x,
            y=y,
            y_errors=y_error,
            label=label,
            marker=marker,
            color=None,
        )

    return plot_data


""" Plot the data of all benchmark files given in 'data_files' """


def plot_benchmark_data(
    logging,
    input_dir,
    det_name,
    file_list,
    label_list,
    title,
    plot_series_name,
    plot_factory,
    out_format,
):

    # Cycle through marker styles per plot
    marker_styles = ["o", "x", "*", "v", "s", "^", "<", ">"]
    marker_style_cycle = itertools.cycle(marker_styles)

    # Save the different plots per hardware backend
    plots = benchmark_plots()

    # Go through all benchmark data files in the list and make a comparison plot
    for i, file in enumerate(file_list):
        # Get the data for the next benchmark case
        _, df = prepare_benchmark_data(logging, input_dir, file)
        marker = next(marker_style_cycle)

        n_tracks, latency, latency_sigma = filter_benchmark_data(df, "real_time")
        _, throughput, throughput_sigma = filter_benchmark_data(df, "TracksPropagated")

        # Initialize plots
        if i == 0:
            # Plot the data against the number of tracks
            latency_plot = plot_benchmark_case(
                plot_factory=plot_factory,
                plot_type="latency",
                label=label_list[i],
                x=n_tracks,
                y=latency,
                y_error=latency_sigma,
                marker=marker,
                title=title,
            )

            throughput_plot = plot_benchmark_case(
                plot_factory=plot_factory,
                plot_type="throughput",
                label=label_list[i],
                x=n_tracks,
                y=throughput,
                y_error=throughput_sigma,
                marker=marker,
                title=title,
            )

            plots = benchmark_plots(latency=latency_plot, throughput=throughput_plot)

        # Add new data to plots
        else:
            plot_benchmark_case(
                plot_factory=plot_factory,
                plot_type="latency",
                label=label_list[i],
                x=n_tracks,
                y=latency,
                y_error=latency_sigma,
                marker=marker,
                plot=plots.latency,
            )

            plot_benchmark_case(
                plot_factory=plot_factory,
                plot_type="throughput",
                label=label_list[i],
                x=n_tracks,
                y=throughput,
                y_error=throughput_sigma,
                marker=marker,
                plot=plots.throughput,
            )

    # Write to disk
    plot_factory.write_plot(
        plots.latency, f"{det_name}_{plot_series_name}_latency", out_format
    )

    plot_factory.write_plot(
        plots.throughput, f"{det_name}_{plot_series_name}_throughput", out_format
    )


""" Plot weak and strong scaling data """


def plot_scaling_data(
    logging,
    input_dir,
    det_name,
    file_list,
    label_list,
    title,
    plot_factory,
    out_format,
    n_threads,
    n_cores,
):
    # Cycle through marker styles per plot
    marker_styles = ["o", "x", "*", "v", "s", "^", "<", ">"]
    marker_style_cycle = itertools.cycle(marker_styles)

    # Save the different plots per algebra plugin
    plots = benchmark_plots()

    # Go through all benchmark data files in the list and make a comparison plot
    for i, file in enumerate(file_list):
        # Get the data for the next benchmark case
        _, df = prepare_benchmark_data(logging, input_dir, file)
        marker = next(marker_style_cycle)

        # Filter the relevant data
        _, latency, latency_sigma = filter_benchmark_data(df, "real_time")

        # Split the data set
        weak_sc_latency = latency[: len(n_threads)]
        strong_sc_latency = latency[len(n_threads) :]

        # Calculate speedups and efficiencies
        weak_sc_efficiency = weak_sc_latency[0] / weak_sc_latency
        strong_sc_speedup = strong_sc_latency[0] / strong_sc_latency

        weak_sc_stddev = latency_sigma[: len(n_threads)]
        strong_sc_stddev = latency_sigma[len(n_threads) :]

        # Gaussian error propagation x/y_i
        err_prob = lambda x, y, err_x, err_y: np.sqrt(
            np.power(err_x / y, 2) + np.power((x * err_y) / np.power(y, 2), 2)
        )

        weak_sc_stddev = err_prob(
            weak_sc_latency[0], weak_sc_latency, weak_sc_stddev[0], weak_sc_stddev
        )
        strong_sc_stddev = err_prob(
            strong_sc_latency[0],
            strong_sc_latency,
            strong_sc_stddev[0],
            strong_sc_stddev,
        )

        # Initialize plots
        if i == 0:
            # Plot the data against the number of tracks
            weak_sc_plot = plot_benchmark_case(
                plot_factory=plot_factory,
                plot_type="weak_scaling",
                label=label_list[i],
                x=n_threads,
                y=weak_sc_efficiency,
                y_error=weak_sc_stddev,
                marker=marker,
                title=title,
                log_scale=2,
                xaxis_format="{x:3.0f}",
                yaxis_format="{x:3.2f}",
            )

            strong_sc_plot = plot_benchmark_case(
                plot_factory=plot_factory,
                plot_type="strong_scaling",
                label=label_list[i],
                x=n_threads,
                y=strong_sc_speedup,
                y_error=strong_sc_stddev,
                marker=marker,
                title=title,
                log_scale=2,
                xaxis_format="{x:3.0f}",
                yaxis_format="{x:3.0f}",
            )

            plots = benchmark_plots(
                weak_scaling=weak_sc_plot, strong_scaling=strong_sc_plot
            )

        # Add new data to plots
        else:
            plot_benchmark_case(
                plot_factory=plot_factory,
                plot_type="weak_scaling",
                label=label_list[i],
                x=n_threads,
                y=weak_sc_efficiency,
                y_error=weak_sc_stddev,
                marker=marker,
                plot=plots.weak_scaling,
            )

            plot_benchmark_case(
                plot_factory=plot_factory,
                plot_type="strong_scaling",
                label=label_list[i],
                x=n_threads,
                y=strong_sc_speedup,
                y_error=strong_sc_stddev,
                marker=marker,
                plot=plots.strong_scaling,
            )

    # Ideal weak scaling
    plot_factory.add_graph(
        plot=plots.weak_scaling,
        x=n_threads,
        y=[1] * len(n_threads),
        marker="",
        color="r",
        label="ideal scaling",
    )

    plot_factory.vertical_line(
        plot_data=plots.weak_scaling,
        x=n_cores,
        y=2 * min(weak_sc_efficiency),
        color="black",
        label="no. cores",
    )

    # Ideal strong scaling
    plot_factory.add_graph(
        plot=plots.strong_scaling,
        x=n_threads,
        y=n_threads,
        marker="",
        color="r",
        label="ideal scaling",
    )

    plot_factory.vertical_line(
        plot_data=plots.strong_scaling, x=n_cores, color="black", label="no. cores"
    )

    # Write to disk
    plot_factory.write_plot(plots.weak_scaling, f"{det_name}_weak_scaling", out_format)

    plot_factory.write_plot(
        plots.strong_scaling, f"{det_name}_strong_scaling", out_format
    )
