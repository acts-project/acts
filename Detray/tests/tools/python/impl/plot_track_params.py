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
import numpy as np
import pandas as pd
import math
import sys

# Common options
lgd_loc = "upper right"

""" Read track position data """


def read_track_data(file, logging):
    if file:
        # Preserve floating point precision
        df = pd.read_csv(file, float_precision="round_trip")
        logging.debug(df)

        return df
    else:
        logging.error("Could not find navigation data file: " + file)
        sys.exit(1)


""" Plot the distributions of track parameter data """


def plot_track_params(opts, detector, track_type, plot_factory, out_format, df):

    detector_name = detector.replace(" ", "_")

    fig_size = (8.5, 8.5)

    # Configure the plot legend
    lgd_ops = plotting.legend_options(
        loc=lgd_loc,
        ncol=4,
        colspacing=0.8,
        handletextpad=0.005,
        horiz_anchor=1.02,
        vert_anchor=1.12,
    )

    # Configure the y axes for all hists
    y_axis_opts = plotting.axis_options(label="")

    # Plot the track data
    def __create_plot(x, bins, suffix, x_axis, y_axis=y_axis_opts):
        hist_data = plot_factory.hist1D(
            x=x,
            bins=bins,
            x_axis=x_axis,
            y_axis=y_axis,
            lgd_ops=lgd_ops,
            show_stats=False,
            figsize=fig_size,
        )

        plot_factory.write_plot(
            hist_data, f"{detector_name}_{track_type}_{suffix}", out_format
        )

    # Plot the charge
    x_axis_opts = plotting.axis_options(label=r"$q\,\mathrm{[e]}$")
    __create_plot(x=df["q"], bins=4, x_axis=x_axis_opts, suffix="charge_dist")

    # Plot the total momentum
    p = np.sqrt(np.square(df["px"]) + np.square(df["py"]) + np.square(df["pz"]))
    x_axis_opts = plotting.axis_options(label=r"$p_{tot}\,\mathrm{[GeV]}$")

    __create_plot(
        x=p,
        bins=1 if np.std(p) < 1e-10 else 100,
        x_axis=x_axis_opts,
        y_axis=y_axis_opts._replace(log_scale=10),
        suffix="p_dist",
    )

    # Plot the transverse momentum
    pT = np.sqrt(np.square(df["px"]) + np.square(df["py"]))
    x_axis_opts = plotting.axis_options(label=r"$p_{T}\,\mathrm{[GeV]}$")

    __create_plot(
        x=pT,
        bins=1 if np.std(pT) < 1e-10 else 100,
        x_axis=x_axis_opts,
        suffix="pT_dist",
    )

    # Plot the x-origin
    x_axis_opts = plotting.axis_options(label=r"$x\,\mathrm{[mm]}$")
    __create_plot(x=df["x"], bins=100, x_axis=x_axis_opts, suffix="x_origin")

    # Plot the y-origin
    x_axis_opts = plotting.axis_options(label=r"$y\,\mathrm{[mm]}$")
    __create_plot(x=df["y"], bins=100, x_axis=x_axis_opts, suffix="y_origin")

    # Plot the z-origin
    x_axis_opts = plotting.axis_options(label=r"$z\,\mathrm{[mm]}$")
    __create_plot(x=df["z"], bins=100, x_axis=x_axis_opts, suffix="z_origin")

    # Plot the phi angle of the track direction
    phi = np.arctan2(df["py"], df["px"])
    x_axis_opts = plotting.axis_options(label=r"$\varphi\,\mathrm{[rad]}$")

    __create_plot(x=phi, bins=100, x_axis=x_axis_opts, suffix="dir_phi")

    # Plot the theta value of the track direction
    theta = np.arctan2(pT, df["pz"])
    x_axis_opts = plotting.axis_options(label=r"$\theta\,\mathrm{[rad]}$")

    __create_plot(x=theta, bins=100, x_axis=x_axis_opts, suffix="dir_theta")

    # Plot the eta value of the track direction
    eta = np.arctanh(df["pz"] / p)
    x_axis_opts = plotting.axis_options(label=r"$\eta$")

    __create_plot(x=eta, bins=100, x_axis=x_axis_opts, suffix="dir_eta")


""" Plot the track positions of two data sources - rz view """


def compare_track_pos_xy(
    opts,
    detector,
    scan_type,
    plot_factory,
    out_format,
    df1,
    label1,
    color1,
    df2,
    label2,
    color2,
):

    n_rays = np.max(df1["track_id"]) + 1
    tracks = "rays" if scan_type == "ray" else "helices"

    # Reduce data to the requested z-range (50mm tolerance)
    min_z = opts.z_range[0]
    max_z = opts.z_range[1]
    assert min_z < max_z, "xy plotting range: min z must be smaller that max z"
    pos_range = lambda data: ((data["z"] > min_z) & (data["z"] < max_z))

    first_x, first_y = plotting.filter_data(
        data=df1, filter=pos_range, variables=["x", "y"]
    )

    second_x, second_y = plotting.filter_data(
        data=df2, filter=pos_range, variables=["x", "y"]
    )

    # Plot the xy coordinates of the filtered track positions
    lgd_ops = plotting.legend_options(
        loc="upper center",
        ncol=4,
        colspacing=0.4,
        handletextpad=0.005,
        horiz_anchor=0.5,
        vert_anchor=1.095,
    )

    hist_data = plot_factory.scatter(
        figsize=(10, 10),
        x=first_x,
        y=first_y,
        x_axis=plotting.axis_options(label=r"$x\,\mathrm{[mm]}$"),
        y_axis=plotting.axis_options(label=r"$y\,\mathrm{[mm]}$"),
        label=label1,
        color=color1,
        alpha=1.0,
        show_stats=lambda x, y: f"{n_rays} {tracks}",
        lgd_ops=lgd_ops,
    )

    # Compare against second data set
    plot_factory.highlight_region(hist_data, second_x, second_y, color2, label2)

    detector_name = detector.replace(" ", "_")
    l1 = label1.replace(" ", "_").replace("(", "").replace(")", "")
    l2 = label2.replace(" ", "_").replace("(", "").replace(")", "")

    # Need a very high dpi to reach a good coverage of the individual points
    plot_factory.write_plot(
        hist_data,
        f"{detector_name}_{scan_type}_track_pos_{l1}_{l2}_xy",
        out_format,
        dpi=600,
    )


""" Plot the track positions of two data sources - rz view """


def compare_track_pos_rz(
    opts,
    detector,
    scan_type,
    plot_factory,
    out_format,
    df1,
    label1,
    color1,
    df2,
    label2,
    color2,
):

    n_rays = np.max(df1["track_id"]) + 1
    tracks = "rays" if scan_type == "ray" else "helices"

    first_x, first_y, first_z = plotting.filter_data(
        data=df1, variables=["x", "y", "z"]
    )

    second_x, second_y, second_z = plotting.filter_data(
        data=df2, variables=["x", "y", "z"]
    )

    # Plot the xy coordinates of the filtered intersections points
    lgd_ops = plotting.legend_options(
        loc="upper center",
        ncol=4,
        colspacing=0.8,
        handletextpad=0.005,
        horiz_anchor=0.5,
        vert_anchor=1.168,
    )

    hist_data = plot_factory.scatter(
        figsize=(12, 6),
        x=first_z,
        y=np.hypot(first_x, first_y),
        x_axis=plotting.axis_options(label=r"$z\,\mathrm{[mm]}$"),
        y_axis=plotting.axis_options(label=r"$r\,\mathrm{[mm]}$"),
        label=label1,
        color=color1,
        alpha=1.0,
        show_stats=lambda x, y: f"{n_rays} {tracks}",
        lgd_ops=lgd_ops,
    )

    # Compare against second data set
    plot_factory.highlight_region(
        hist_data, second_z, np.hypot(second_x, second_y), color2, label2
    )

    detector_name = detector.replace(" ", "_")
    l1 = label1.replace(" ", "_").replace("(", "").replace(")", "")
    l2 = label2.replace(" ", "_").replace("(", "").replace(")", "")

    # Need a very high dpi to reach a good coverage of the individual points
    plot_factory.write_plot(
        hist_data,
        f"{detector_name}_{scan_type}_track_pos_{l1}_{l2}_rz",
        out_format,
        dpi=600,
    )


""" Plot the absolute track positions distance """


def plot_track_pos_dist(
    opts, detector, scan_type, plot_factory, out_format, df1, label1, df2, label2
):

    dist = np.sqrt(
        np.square(df1["x"] - df2["x"])
        + np.square(df1["y"] - df2["y"])
        + np.square(df1["z"] - df2["z"])
    )

    dist_outlier = math.sqrt(3 * opts.outlier**2)

    # Remove outliers
    filter_dist = np.absolute(dist) < dist_outlier
    filtered_dist = dist[filter_dist]

    if not np.all(filter_dist == True):
        print("\nRemoved outliers (dist):")
        for i, d in enumerate(dist):
            if math.fabs(d) > dist_outlier:
                track_id = (df1["track_id"].to_numpy())[i]
                print(f"track {track_id}: {d}")

    # Where to place the legend box
    lgd_ops = plotting.legend_options(
        loc=lgd_loc,
        ncol=4,
        colspacing=0.8,
        handletextpad=0.005,
        horiz_anchor=1.02,
        vert_anchor=1.24,
    )

    # Plot the xy coordinates of the filtered intersections points
    hist_data = plot_factory.hist1D(
        x=filtered_dist,
        bins=100,
        x_axis=plotting.axis_options(label=r"$d\,\mathrm{[mm]}$"),
        y_axis=plotting.axis_options(label="", log_scale=10),
        figsize=(8.5, 8.5),
        lgd_ops=lgd_ops,
    )

    detector_name = detector.replace(" ", "_")
    l1 = label1.replace(" ", "_").replace("(", "").replace(")", "")
    l2 = label2.replace(" ", "_").replace("(", "").replace(")", "")
    plot_factory.write_plot(
        hist_data, f"{detector_name}_{scan_type}_dist_{l1}_{l2}", out_format
    )


""" Plot the track position residual for the given variable """


def plot_track_pos_res(
    opts, detector, scan_type, plot_factory, out_format, df1, label1, df2, label2, var
):

    tracks = "rays" if scan_type == "ray" else "helices"

    assert len(df1[var]) == len(df2[var])
    res = df1[var] - df2[var]

    # Remove outliers
    filter_res = np.absolute(res) < opts.outlier
    filtered_res = res[filter_res]

    u_out = o_out = int(0)
    if not np.all(filter_res == True):
        print(f"\nRemoved outliers ({var}):")
        for i, r in enumerate(res):
            if math.fabs(r) > opts.outlier:
                track_id = (df1["track_id"].to_numpy())[i]
                print(f"track {track_id}: {df1[var][i]} - {df2[var][i]} = {r}")

                if r < 0.0:
                    u_out = u_out + 1
                else:
                    o_out = o_out + 1

    lgd_ops = plotting.legend_options(
        loc=lgd_loc,
        ncol=4,
        colspacing=0.01,
        handletextpad=0.0005,
        horiz_anchor=1.02,
        vert_anchor=1.28,
    )

    # Plot the residuals as a histogram and fit a gaussian to it
    hist_data = plot_factory.hist1D(
        x=filtered_res,
        figsize=(9, 9),
        bins=100,
        x_axis=plotting.axis_options(
            label=r"$\mathrm{res}" + rf"\,{var}" + r"\,\mathrm{[mm]}$"
        ),
        lgd_ops=lgd_ops,
        u_outlier=u_out,
        o_outlier=o_out,
    )

    mu, sig = plot_factory.fit_gaussian(hist_data)
    if mu is None or sig is None:
        print(rf"WARNING: fit failed (res ({tracks}): {label1} - {label2} )")

    detector_name = detector.replace(" ", "_")
    l1 = label1.replace(" ", "_").replace("(", "").replace(")", "")
    l2 = label2.replace(" ", "_").replace("(", "").replace(")", "")

    plot_factory.write_plot(
        hist_data, f"{detector_name}_{scan_type}_track_res_{var}_{l1}_{l2}", out_format
    )
