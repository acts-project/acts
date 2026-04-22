# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# detray includes
import plotting

# python includes
import math
import numpy as np
import os
import pandas as pd

# Common plot labels
label_eta = r"$\eta$"
label_phi = r"$\phi\,\mathrm{[rad]}$"
label_thickness_x0 = r"thickness / $X_0$"
label_path_x0 = r"path length / $X_0$"
label_thickness_l0 = r"thickness / $\Lambda_0$"
label_path_l0 = r"path length / $\Lambda_0$"

# Common options
ldg_loc = "upper center"

""" Read the material scan data from file and prepare data frame """


def read_material_data(inputdir, logging, det_name, read_cuda):

    # Input data directory
    data_dir = os.fsencode(inputdir)

    material_scan_file = cpu_material_trace_file = cuda_material_trace_file = ""

    # Find the data files by naming convention
    for file in os.listdir(data_dir):
        filename = os.fsdecode(file)

        if filename.find(det_name + "_material_scan") != -1:
            material_scan_file = inputdir + "/" + filename
        elif filename.find(det_name + "_navigation_material_trace_cpu") != -1:
            cpu_material_trace_file = inputdir + "/" + filename
        elif (
            read_cuda
            and filename.find(det_name + "_navigation_material_trace_cuda") != -1
        ):
            cuda_material_trace_file = inputdir + "/" + filename

    df_scan = pd.read_csv(material_scan_file, float_precision="round_trip")
    df_cpu_trace = pd.read_csv(cpu_material_trace_file, float_precision="round_trip")
    df_cuda_trace = pd.DataFrame({})
    if read_cuda:
        df_cuda_trace = pd.read_csv(
            cuda_material_trace_file, float_precision="round_trip"
        )

    return df_scan, df_cpu_trace, df_cuda_trace


""" Calculate edges of bins to plot the material data """


def get_n_bins(df):
    # Find the number of ray directions
    row_count = df.groupby(df["eta"]).count()
    y_bins = row_count["phi"].max()
    x_bins = int(len(df["eta"]) / y_bins)
    assert (
        len(df["eta"]) == x_bins * y_bins
    ), "Could not infer the number of rays correctly"

    # Get the axis spacing
    x_range = np.max(df["eta"]) - np.min(df["eta"])
    x_binning = np.linspace(
        np.min(df["eta"]) - 0.5 * x_range / x_bins,
        np.max(df["eta"]) + 0.5 * x_range / x_bins,
        x_bins + 1,
    )

    y_range = np.max(df["phi"]) - np.min(df["phi"])
    y_binning = np.linspace(
        np.min(df["phi"]) - 0.5 * y_range / y_bins,
        np.max(df["phi"]) + 0.5 * y_range / y_bins,
        y_bins + 1,
    )

    return x_binning, y_binning


""" Calculate the binwise errors: Standard Error on the Mean """


def get_errors(df, n, name):
    # Number of entries per bin
    errors = []
    # Iterate over data per bin
    for i in range(0, n):
        # Project the next n rows of the data frame (the bin content is the
        # mean of the material along phi)
        bin_data = df.iloc[i * n : (i + 1) * n]
        # Calculate the error on the mean: stddev/sqrt(n)
        errors.append(np.std(bin_data[name], axis=0) / math.sqrt(n))

    return errors


""" Plot the material thickenss vs phi and eta in units of X_0 """


def X0_vs_eta_phi(df, label, detector, plot_factory, out_format="pdf"):

    # Histogram bin edges
    x_binning, y_binning = get_n_bins(df)

    # Plot the thickness of every material slab in units of X_0
    hist_data = plot_factory.hist2D(
        x=df["eta"],
        y=df["phi"],
        z=df["mat_tX0"],
        label=detector,
        x_axis=plotting.axis_options(label=label_eta),
        y_axis=plotting.axis_options(label=label_phi),
        z_axis=plotting.axis_options(label=label_thickness_x0),
        x_bins=x_binning,
        y_bins=y_binning,
        figsize=(9, 7),
        show_stats=False,
    )

    plot_factory.write_plot(hist_data, detector + "_" + label + "_t_X0_map", out_format)

    # Plot path length through material of the respective ray in units of X_0
    hist_data = plot_factory.hist2D(
        x=df["eta"],
        y=df["phi"],
        z=df["mat_sX0"],
        label=detector,
        x_axis=plotting.axis_options(label=label_eta),
        y_axis=plotting.axis_options(label=label_phi),
        z_axis=plotting.axis_options(label=label_path_x0),
        x_bins=x_binning,
        y_bins=y_binning,
        figsize=(9, 7),
        show_stats=False,
    )

    plot_factory.write_plot(hist_data, detector + "_" + label + "_s_X0_map", out_format)


""" Plot the material thickenss vs phi and eta in units of L_0 """


def L0_vs_eta_phi(df, label, detector, plot_factory, out_format="pdf"):

    # Histogram bin edges
    x_binning, y_binning = get_n_bins(df)

    # Plot the thickness of every material slab in units of L_0
    hist_data = plot_factory.hist2D(
        x=df["eta"],
        y=df["phi"],
        z=df["mat_tL0"],
        label=detector,
        x_axis=plotting.axis_options(label=label_eta),
        y_axis=plotting.axis_options(label=label_phi),
        z_axis=plotting.axis_options(label=label_thickness_l0),
        x_bins=x_binning,
        y_bins=y_binning,
        figsize=(9, 7),
        show_stats=False,
    )

    plot_factory.write_plot(hist_data, detector + "_" + label + "_t_L0_map", out_format)

    # Plot path length through material of the respective ray in units of L_0
    hist_data = plot_factory.hist2D(
        x=df["eta"],
        y=df["phi"],
        z=df["mat_sL0"],
        label=detector,
        x_axis=plotting.axis_options(label=label_eta),
        y_axis=plotting.axis_options(label=label_phi),
        z_axis=plotting.axis_options(label=label_path_l0),
        x_bins=x_binning,
        y_bins=y_binning,
        figsize=(9, 7),
        show_stats=False,
    )

    plot_factory.write_plot(hist_data, detector + "_" + label + "_s_L0_map", out_format)


""" Plot the material thickness in units of X_0 vs eta """


def X0_vs_eta(df, label, detector, plot_factory, out_format="pdf"):
    # Where to place the legend box
    lgd_ops = plotting.legend_options(
        loc=ldg_loc, horiz_anchor=0.725, vert_anchor=1.145
    )

    # Histogram bin edges
    x_binning, y_binning = get_n_bins(df)

    # Same number of entries in every bin as per uniform ray scan
    n_phi = len(y_binning) - 1

    hist_data = plot_factory.hist1D(
        x=df["eta"],
        w=df["mat_tX0"] / n_phi,
        errors=get_errors(df, n_phi, "mat_tX0"),
        normalize=False,
        label=rf"{detector}",
        x_axis=plotting.axis_options(label=label_eta),
        y_axis=plotting.axis_options(label=label_thickness_x0),
        bins=x_binning,
        show_stats=False,
        figsize=(9, 7),
        lgd_ops=lgd_ops,
    )

    plot_factory.write_plot(hist_data, detector + "_" + label + "_t_X0", out_format)

    hist_data = plot_factory.hist1D(
        x=df["eta"],
        w=df["mat_sX0"] / n_phi,
        errors=get_errors(df, n_phi, "mat_sX0"),
        normalize=False,
        label=rf"{detector}",
        x_axis=plotting.axis_options(label=label_eta),
        y_axis=plotting.axis_options(label=label_path_x0),
        bins=x_binning,
        show_stats=False,
        figsize=(9, 7),
        lgd_ops=lgd_ops,
    )

    plot_factory.write_plot(hist_data, detector + "_" + label + "_s_X0", out_format)


""" Plot the material thickness in units of L_0 vs eta """


def L0_vs_eta(df, label, detector, plot_factory, out_format="pdf"):
    # Where to place the legend box
    lgd_ops = plotting.legend_options(
        loc=ldg_loc, horiz_anchor=0.725, vert_anchor=1.145
    )

    # Histogram bin edges
    x_binning, y_binning = get_n_bins(df)

    # Same number of entries in every bin as per uniform ray scan
    n_phi = len(y_binning) - 1

    hist_data = plot_factory.hist1D(
        x=df["eta"],
        w=df["mat_tL0"] / n_phi,
        errors=get_errors(df, n_phi, "mat_tL0"),
        normalize=False,
        label=rf"{detector}",
        x_axis=plotting.axis_options(label=label_eta),
        y_axis=plotting.axis_options(label=label_thickness_l0),
        bins=x_binning,
        show_stats=False,
        figsize=(9, 7),
        lgd_ops=lgd_ops,
    )

    plot_factory.write_plot(hist_data, detector + "_" + label + "_t_L0", out_format)

    hist_data = plot_factory.hist1D(
        x=df["eta"],
        w=df["mat_sL0"] / n_phi,
        errors=get_errors(df, n_phi, "mat_sL0"),
        normalize=False,
        label=rf"{detector}",
        x_axis=plotting.axis_options(label=label_eta),
        y_axis=plotting.axis_options(label=label_path_l0),
        bins=x_binning,
        show_stats=False,
        figsize=(9, 7),
        lgd_ops=lgd_ops,
    )

    plot_factory.write_plot(hist_data, detector + "_" + label + "_s_L0", out_format)


""" Compare two material distributions """


def compare_mat(df_truth, df_rec, label, detector, plot_factory, out_format="pdf"):
    # Where to place the legend box
    lgd_ops = plotting.legend_options(loc=ldg_loc, horiz_anchor=0.5, vert_anchor=1.29)

    # Histogram bin edges
    x_binning, y_binning = get_n_bins(df_truth)

    # Same number of entries in every bin as per uniform ray scan
    n_phi = len(y_binning) - 1

    truth_data = plot_factory.hist1D(
        x=df_truth["eta"],
        w=df_truth["mat_sX0"] / n_phi,
        errors=get_errors(df_truth, n_phi, "mat_sX0"),
        normalize=False,
        label=rf"{detector}: scan",
        x_axis=plotting.axis_options(label=label_eta),
        y_axis=plotting.axis_options(label=label_path_x0),
        bins=x_binning,
        show_stats=False,
        figsize=(10, 10),
        layout="tight",
        lgd_ops=lgd_ops,
    )

    # Add recorded data for comparison
    rec_data = plot_factory.add_hist(
        old_hist=truth_data,
        x=df_rec["eta"].to_numpy(dtype=np.double),
        w=df_rec["mat_sX0"] / n_phi,
        errors=get_errors(df_rec, n_phi, "mat_sX0"),
        normalize=False,
        label=rf"{detector}: navigator",
    )

    # Add a ratio plot to hist_data
    ratio_data = plot_factory.add_ratio(
        nom=truth_data,
        denom=rec_data,
        label="scan/navigation",
        set_log=False,
        show_error=True,
    )

    plot_factory.write_plot(
        ratio_data, detector + "_" + label + "_comparison_s_X0", out_format
    )
