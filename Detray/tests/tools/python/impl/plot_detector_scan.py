# Detray library, part of the ACTS project (R&D line)
#
# (c) 2023-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# detray includes
import plotting

# python includes
import numpy as np
import pandas as pd
import math
import os
import sys

# Common options
lgd_loc = "upper right"

""" Read the detector scan data from files and prepare data frames """


def read_detector_scan_data(intersection_file, track_param_file, logging):
    if intersection_file:
        inters_df = pd.read_csv(intersection_file, float_precision="round_trip")
        trk_param_df = pd.read_csv(track_param_file, float_precision="round_trip")
        scan_df = pd.concat([inters_df, trk_param_df], axis=1)

        logging.debug(scan_df)
    else:
        logging.warning("Could not find ray scan data: " + intersection_file)
        scan_df = pd.DataFrame({})

    return scan_df


""" Read intersection data """


def read_intersection_data(file, logging):
    if file:
        # Preserve floating point precision
        df = pd.read_csv(file, float_precision="round_trip")
        logging.debug(df)

        return df
    else:
        logging.error("Could not find intersection data file: " + file)
        sys.exit(1)


""" Plot the intersection points of the detector with the rays - xy view """


def plot_intersection_points_xy(
    opts, df, detector, scan_type, plot_factory, out_format="png"
):

    n_rays = np.max(df["track_id"]) + 1
    tracks = "rays" if scan_type == "ray" else "helices"

    # Reduce data to the requested z-range (50mm tolerance)
    min_z = opts.z_range[0]
    max_z = opts.z_range[1]
    assert min_z < max_z, "xy plotting range: min z must be smaller that max z"
    sensitive_range = lambda data: (
        (data["z"] > min_z) & (data["z"] < max_z) & (data["type"] == 1)
    )
    portal_range = lambda data: (
        (data["z"] > min_z) & (data["z"] < max_z) & (data["type"] == 0)
    )
    passive_range = lambda data: (
        (data["z"] > min_z) & (data["z"] < max_z) & (data["type"] == 2)
    )

    senstive_x, senstive_y = plotting.filter_data(
        data=df, filter=sensitive_range, variables=["x", "y"]
    )

    # Plot the xy coordinates of the filtered intersections points
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
        x=senstive_x,
        y=senstive_y,
        x_axis=plotting.axis_options(label=r"$x\,\mathrm{[mm]}$"),
        y_axis=plotting.axis_options(label=r"$y\,\mathrm{[mm]}$"),
        label="sensitives",
        color="C5",
        show_stats=lambda x, y: f"{n_rays} {tracks}",
        lgd_ops=lgd_ops,
    )

    # Portal surfaces
    if not opts.hide_portals:
        portal_x, portal_y = plotting.filter_data(
            data=df, filter=portal_range, variables=["x", "y"]
        )

        plot_factory.highlight_region(hist_data, portal_x, portal_y, "C0", "portals")

    # Passive surfaces
    if not opts.hide_passives:
        passive_x, passive_y = plotting.filter_data(
            data=df, filter=passive_range, variables=["x", "y"]
        )

        plot_factory.highlight_region(hist_data, passive_x, passive_y, "C2", "passives")

    # Set aspect ratio
    hist_data.axes.set_aspect("equal")

    detector_name = detector.replace(" ", "_")
    plot_factory.write_plot(
        hist_data, f"{detector_name}_{scan_type}_scan_xy", out_format
    )


""" Plot the intersection points of the detector with the rays - rz view """


def plot_intersection_points_rz(
    opts, df, detector, scan_type, plot_factory, out_format="png"
):

    n_rays = np.max(df["track_id"]) + 1
    tracks = "rays" if scan_type == "ray" else "helices"

    # Reduce data to the requested z-range
    sensitive_range = lambda data: (data["type"] == 1)
    portal_range = lambda data: (data["type"] == 0)
    passive_range = lambda data: (data["type"] == 2)

    sensitive_x, sensitive_y, sensitive_z = plotting.filter_data(
        data=df, filter=sensitive_range, variables=["x", "y", "z"]
    )

    # Plot the xy coordinates of the filtered intersections points
    lgd_ops = plotting.legend_options(
        loc="upper center",
        ncol=4,
        colspacing=0.8,
        handletextpad=0.005,
        horiz_anchor=0.5,
        vert_anchor=1.165,
    )

    hist_data = plot_factory.scatter(
        figsize=(12, 6),
        x=sensitive_z,
        y=np.hypot(sensitive_x, sensitive_y),
        x_axis=plotting.axis_options(label=r"$z\,\mathrm{[mm]}$"),
        y_axis=plotting.axis_options(label=r"$r\,\mathrm{[mm]}$"),
        label="sensitives",
        color="C5",
        show_stats=lambda x, y: f"{n_rays} {tracks}",
        lgd_ops=lgd_ops,
    )

    # Portal surfaces
    if not opts.hide_portals:
        portal_x, portal_y, portal_z = plotting.filter_data(
            data=df, filter=portal_range, variables=["x", "y", "z"]
        )

        plot_factory.highlight_region(
            hist_data, portal_z, np.hypot(portal_x, portal_y), "C0", "portals"
        )

    # Passive surfaces
    if not opts.hide_passives:
        passive_x, passive_y, passive_z = plotting.filter_data(
            data=df, filter=passive_range, variables=["x", "y", "z"]
        )

        plot_factory.highlight_region(
            hist_data, passive_z, np.hypot(passive_x, passive_y), "C2", "passives"
        )

    detector_name = detector.replace(" ", "_")
    plot_factory.write_plot(
        hist_data, f"{detector_name}_{scan_type}_scan_rz", out_format
    )


""" Plot the intersection local position residual for the given variable """


def plot_intersection_pos_res(
    opts, detector, plot_factory, scan_type, df1, label1, df2, label2, var, out_format
):

    tracks = "rays" if scan_type == "ray" else "helices"

    # Filter the relevant data from the frame (sensitive = 1, hole = 15)
    is_sensitive = lambda data_frame: (
        (data_frame["type"] == 1) | (data_frame["type"] == 15)
    )

    var_truth, track_ids = plotting.filter_data(
        data=df1, filter=is_sensitive, variables=[var, "track_id"]
    )
    var_nav = plotting.filter_data(data=df2, filter=is_sensitive, variables=[var])

    assert len(var_truth) == len(var_nav)
    res = var_truth - var_nav

    # Remove outliers (happens when comparing a hole with a valid intersection)
    filter_res = np.absolute(res) < opts.outlier
    filtered_res = res[filter_res]

    u_out = o_out = int(0)
    if not np.all(filter_res == True):
        print(f"\nRemoved outliers ({var}):")
        for i, r in enumerate(res):
            if math.fabs(r) > opts.outlier:
                print(f"track {track_ids[i]}: {var_truth[i]} - {var_nav[i]} = {r}")

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
    var_for_label = var.replace("_", "\,")
    hist_data = plot_factory.hist1D(
        x=filtered_res,
        figsize=(9, 9),
        bins=100,
        x_axis=plotting.axis_options(
            label=r"$\mathrm{res}" + rf"~{var_for_label}" + r"\,\mathrm{[mm]}$"
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
        hist_data, f"{detector_name}_{scan_type}_intr_res_{var}_{l1}_{l2}", out_format
    )
