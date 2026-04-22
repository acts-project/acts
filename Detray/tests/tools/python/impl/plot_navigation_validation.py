# Detray library, part of the ACTS project (R&D line)
#
# (c) 2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

from .plot_detector_scan import (
    read_detector_scan_data,
    read_intersection_data,
    plot_intersection_points_xy,
    plot_intersection_points_rz,
    plot_intersection_pos_res,
)
from .plot_track_params import (
    read_track_data,
    compare_track_pos_xy,
    compare_track_pos_rz,
    plot_track_pos_dist,
    plot_track_pos_res,
)

# python includes
import pandas as pd
import os

""" Construct a string in the format of the output file names and compare it"""


def __comp_filename(filename, det_name, file_stem, p_min="", p_max=""):
    check_file_stem = f"{det_name}_{file_stem}" in filename
    check_p_min = f"_{p_min}" in filename
    check_p_max = f"_{p_max}" in filename

    return check_file_stem and check_p_min and check_p_max


""" Read the detector scan data from files and prepare data frames """


def read_scan_data(logging, inputdir, det_name, p_min, p_max):

    # Input data directory
    data_dir = os.fsencode(inputdir)

    ray_scan_intersections_file = ray_scan_track_param_file = ""
    helix_scan_intersections_file = helix_scan_track_param_file = ""

    # Find the data files by naming convention
    for file in os.listdir(data_dir):
        filename = os.fsdecode(file)

        if __comp_filename(filename, det_name, "ray_scan_intersections"):
            ray_scan_intersections_file = inputdir + "/" + filename
        elif __comp_filename(filename, det_name, "ray_scan_track_parameters"):
            ray_scan_track_param_file = inputdir + "/" + filename
        elif __comp_filename(
            filename, det_name, "helix_scan_intersections", p_min, p_max
        ):
            helix_scan_intersections_file = inputdir + "/" + filename
        elif __comp_filename(
            filename, det_name, "helix_scan_track_parameters", p_min, p_max
        ):
            helix_scan_track_param_file = inputdir + "/" + filename

    # Read ray scan data
    ray_scan_df = read_detector_scan_data(
        ray_scan_intersections_file, ray_scan_track_param_file, logging
    )

    # Read helix scan data
    helix_scan_df = read_detector_scan_data(
        helix_scan_intersections_file, helix_scan_track_param_file, logging
    )

    return ray_scan_df, helix_scan_df


""" Read the recorded track positions from files and prepare data frames """


def read_navigation_intersection_data(
    logging, inputdir, det_name, p_min, p_max, read_cuda
):

    # Input data directory
    data_dir = os.fsencode(inputdir)

    ray_truth_file = ray_data_file = ray_data_cuda_file = ""
    helix_truth_file = helix_data_file = helix_data_cuda_file = ""

    # Find the data files by naming convention
    for file in os.listdir(data_dir):
        filename = os.fsdecode(file)

        if read_cuda and __comp_filename(
            filename, det_name, "ray_navigation_intersections_cuda"
        ):
            ray_data_cuda_file = inputdir + "/" + filename
        elif __comp_filename(filename, det_name, "ray_navigation_intersections"):
            ray_data_file = inputdir + "/" + filename
        elif __comp_filename(filename, det_name, "ray_truth_intersections"):
            ray_truth_file = inputdir + "/" + filename
        elif read_cuda and __comp_filename(
            filename, det_name, "helix_navigation_intersections_cuda", p_min, p_max
        ):
            helix_data_cuda_file = inputdir + "/" + filename
        elif __comp_filename(
            filename, det_name, "helix_navigation_intersections", p_min, p_max
        ):
            helix_data_file = inputdir + "/" + filename
        elif __comp_filename(
            filename, det_name, "helix_truth_intersections", p_min, p_max
        ):
            helix_truth_file = inputdir + "/" + filename

    ray_df = read_intersection_data(ray_data_file, logging)
    ray_truth_df = read_intersection_data(ray_truth_file, logging)
    helix_df = read_intersection_data(helix_data_file, logging)
    helix_truth_df = read_intersection_data(helix_truth_file, logging)

    ray_cuda_df = helix_cuda_df = pd.DataFrame({})
    if read_cuda:
        ray_cuda_df = read_intersection_data(ray_data_cuda_file, logging)
        helix_cuda_df = read_intersection_data(helix_data_cuda_file, logging)

    return ray_df, ray_truth_df, ray_cuda_df, helix_df, helix_truth_df, helix_cuda_df


""" Read the recorded track positions from files and prepare data frames """


def read_navigation_track_data(logging, inputdir, det_name, p_min, p_max, read_cuda):

    # Input data directory
    data_dir = os.fsencode(inputdir)

    ray_truth_file = ray_data_file = ray_data_cuda_file = ""
    helix_truth_file = helix_data_file = helix_data_cuda_file = ""

    # Find the data files by naming convention
    for file in os.listdir(data_dir):
        filename = os.fsdecode(file)

        if read_cuda and __comp_filename(
            filename, det_name, "ray_navigation_track_params_cuda"
        ):
            ray_data_cuda_file = inputdir + "/" + filename
        elif __comp_filename(filename, det_name, "ray_navigation_track_params"):
            ray_data_file = inputdir + "/" + filename
        elif __comp_filename(filename, det_name, "ray_truth_track_params"):
            ray_truth_file = inputdir + "/" + filename
        elif read_cuda and __comp_filename(
            filename, det_name, "helix_navigation_track_params_cuda", p_min, p_max
        ):
            helix_data_cuda_file = inputdir + "/" + filename
        elif __comp_filename(
            filename, det_name, "helix_navigation_track_params", p_min, p_max
        ):
            helix_data_file = inputdir + "/" + filename
        elif __comp_filename(
            filename, det_name, "helix_truth_track_params", p_min, p_max
        ):
            helix_truth_file = inputdir + "/" + filename

    ray_df = read_track_data(ray_data_file, logging)
    ray_truth_df = read_track_data(ray_truth_file, logging)
    helix_df = read_track_data(helix_data_file, logging)
    helix_truth_df = read_track_data(helix_truth_file, logging)

    ray_cuda_df = helix_cuda_df = pd.DataFrame({})
    if read_cuda:
        ray_cuda_df = read_track_data(ray_data_cuda_file, logging)
        helix_cuda_df = read_track_data(helix_data_cuda_file, logging)

    return ray_df, ray_truth_df, ray_cuda_df, helix_df, helix_truth_df, helix_cuda_df


""" Plot the intersection data gathered during the detector scan """


def plot_detector_scan_data(
    args, det_name, plot_factory, data_type, df_scan, out_format="png"
):

    # Plot truth scan
    plot_intersection_points_xy(
        args, df_scan, det_name, data_type, plot_factory, out_format
    )
    plot_intersection_points_rz(
        args, df_scan, det_name, data_type, plot_factory, out_format
    )


""" Plot the intersection data gathered during the navigation validation """


def plot_navigation_intersection_data(
    args, det_name, plot_factory, data_type, df_scan, df_nav, label, out_format="png"
):
    # Plot the residuals in local 0 and local 1
    plot_intersection_pos_res(
        args,
        det_name,
        plot_factory,
        data_type,
        df_scan,
        "truth",
        df_nav,
        label,
        "loc_0",
        out_format,
    )
    plot_intersection_pos_res(
        args,
        det_name,
        plot_factory,
        data_type,
        df_scan,
        "truth",
        df_nav,
        label,
        "loc_1",
        out_format,
    )


""" Plot the track data gathered during the navigation validation """


def plot_navigation_track_data(
    args,
    det_name,
    plot_factory,
    data_type,
    df_truth,
    truth_name,
    df_ref,
    ref_name,
    out_format="png",
):

    # xy
    compare_track_pos_xy(
        args,
        det_name,
        data_type,
        plot_factory,
        out_format,
        df_truth,
        truth_name,
        "r",
        df_ref,
        ref_name,
        "darkgrey",
    )
    # rz
    compare_track_pos_rz(
        args,
        det_name,
        data_type,
        plot_factory,
        out_format,
        df_truth,
        truth_name,
        "r",
        df_ref,
        ref_name,
        "darkgrey",
    )

    # Absolute distance
    plot_track_pos_dist(
        args,
        det_name,
        data_type,
        plot_factory,
        out_format,
        df_truth,
        truth_name,
        df_ref,
        ref_name,
    )

    # Residuals
    plot_track_pos_res(
        args,
        det_name,
        data_type,
        plot_factory,
        out_format,
        df_truth,
        truth_name,
        df_ref,
        ref_name,
        "x",
    )
    plot_track_pos_res(
        args,
        det_name,
        data_type,
        plot_factory,
        out_format,
        df_truth,
        truth_name,
        df_ref,
        ref_name,
        "y",
    )
    plot_track_pos_res(
        args,
        det_name,
        data_type,
        plot_factory,
        out_format,
        df_truth,
        truth_name,
        df_ref,
        ref_name,
        "z",
    )
