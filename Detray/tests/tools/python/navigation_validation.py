# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# detray imports
from impl import (
    read_scan_data,
    read_navigation_intersection_data,
    read_navigation_track_data,
)
from impl import (
    plot_detector_scan_data,
    plot_navigation_intersection_data,
    plot_navigation_track_data,
)
from impl import plot_track_params
from options import (
    common_options,
    detector_io_options,
    random_track_generator_options,
    propagation_options,
    plotting_options,
)
from options import (
    parse_common_options,
    parse_detector_io_options,
    parse_plotting_options,
)
from plotting import pyplot_factory as plt_factory
from utils import read_detector_name, get_p_range
from utils import add_track_generator_args, add_propagation_args, add_detector_io_args

# python imports
import argparse
import os
import subprocess
import sys
import json


def __main__():

    # ---------------------------------------------------------------arg parsing

    descr = "Detray Navigation Validation"

    # Define options
    parent_parsers = [
        common_options(descr),
        detector_io_options(),
        random_track_generator_options(),
        propagation_options(),
        plotting_options(),
    ]

    parser = argparse.ArgumentParser(description=descr, parents=parent_parsers)

    parser.add_argument(
        "--bindir",
        "-bin",
        help=("Directory containing the validation executables"),
        default="./bin",
        type=str,
    )
    parser.add_argument(
        "--datadir",
        "-data",
        help=("Directory containing the data files"),
        default="./validation_data",
        type=str,
    )
    parser.add_argument(
        "--cuda",
        help=("Run the CUDA navigation validation."),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--sycl",
        help=("Run the SYCL navigation validation."),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--overlaps_tol",
        "-ot",
        help=("Tolerance for considering surfaces to be overlapping [mm]"),
        default=0.0001,
        type=float,
    )
    parser.add_argument(
        "--z_range",
        "-zrng",
        nargs=2,
        help=("z range for the xy-view [mm]."),
        default=[-50, 50],
        type=float,
    )
    parser.add_argument(
        "--hide_portals",
        help=("Hide portal surfaces in plots."),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--hide_passives",
        help=("Hide passive surfaces in plots."),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--outlier",
        "-out",
        help=("Threshold for outliers in residual plots [mm]."),
        default=1,
        type=float,
    )

    # Parse options
    args = parser.parse_args()

    logging = parse_common_options(args, descr)
    parse_detector_io_options(args, logging)
    _, out_dir, out_format = parse_plotting_options(args, logging)

    # IO path for data files
    datadir = args.datadir.strip("/")

    # Check bin path
    bindir = args.bindir.strip("/")
    cpu_validation = bindir + "/detray_navigation_validation"
    cuda_validation = bindir + "/detray_navigation_validation_cuda"

    if not os.path.isdir(bindir) or not os.path.isfile(cpu_validation):
        logging.error(f"Navigation validation binaries were not found! ({args.bindir})")
        sys.exit(1)

    # -----------------------------------------------------------------------run

    # Pass on the options for the validation tools
    args_list = [
        "--data_dir",
        datadir,
        "--overlaps_tol",
        str(args.overlaps_tol),
    ]

    # Add parsed options to argument list
    add_detector_io_args(args_list, args)
    add_track_generator_args(args_list, args)
    add_propagation_args(args_list, args)

    logging.debug(args_list)

    # Run the host validation and produce the truth data
    logging.debug("Running CPU validation")
    subprocess.run([cpu_validation, "--write_scan_data"] + args_list)

    # Run the device validation (if it has been built)
    if args.cuda and os.path.isfile(cuda_validation):
        logging.debug("Running CUDA validation")
        subprocess.run([cuda_validation] + args_list)

    elif args.cuda:
        logging.error("Could not find CUDA navigation validation executable")

    if args.sycl:
        logging.error("SYCL validation is not implemented")

    # ----------------------------------------------------------------------plot

    logging.info("Generating data plots...\n")

    det_name = read_detector_name(args.geometry_file, logging)
    logging.debug("Detector: " + det_name)

    # Check the data path (should have been created when running the validation)
    if not os.path.isdir(datadir):
        logging.error(f"Data directory was not found! ({args.datadir})")
        sys.exit(1)

    plot_factory = plt_factory(out_dir, logging)

    # Read the truth data
    p_min, p_max = get_p_range(args, logging)
    ray_scan_df, helix_scan_df = read_scan_data(
        logging, datadir, det_name, p_min, p_max
    )

    # Plot detector scan data
    plot_detector_scan_data(
        args, det_name, plot_factory, "ray", ray_scan_df, out_format
    )
    plot_detector_scan_data(
        args, det_name, plot_factory, "helix", helix_scan_df, out_format
    )

    # Read the recorded intersection data
    (
        ray_nav_intr_df,
        ray_nav_intr_truth_df,
        ray_nav_intr_cuda_df,
        helix_nav_intr_df,
        helix_nav_intr_truth_df,
        helix_nav_intr_cuda_df,
    ) = read_navigation_intersection_data(
        logging, datadir, det_name, p_min, p_max, args.cuda
    )

    # Plot intersection data
    label_cpu = "navigation (CPU)"
    label_cuda = "navigation (CUDA)"

    plot_navigation_intersection_data(
        args,
        det_name,
        plot_factory,
        "ray",
        ray_nav_intr_truth_df,
        ray_nav_intr_df,
        label_cpu,
        out_format,
    )

    plot_navigation_intersection_data(
        args,
        det_name,
        plot_factory,
        "helix",
        helix_nav_intr_truth_df,
        helix_nav_intr_df,
        label_cpu,
        out_format,
    )

    if args.cuda:
        plot_navigation_intersection_data(
            args,
            det_name,
            plot_factory,
            "ray",
            ray_nav_intr_truth_df,
            ray_nav_intr_cuda_df,
            label_cuda,
            out_format,
        )

        plot_navigation_intersection_data(
            args,
            det_name,
            plot_factory,
            "helix",
            helix_nav_intr_truth_df,
            helix_nav_intr_cuda_df,
            label_cuda,
            out_format,
        )

    # Plot distributions of track parameter values
    # Only take initial track parameters from generator
    ray_intial_trk_df = ray_scan_df.drop_duplicates(subset=["track_id"])
    helix_intial_trk_df = helix_scan_df.drop_duplicates(subset=["track_id"])
    plot_track_params(
        args, det_name, "helix", plot_factory, out_format, helix_intial_trk_df
    )
    plot_track_params(
        args, det_name, "ray", plot_factory, out_format, ray_intial_trk_df
    )

    # Read the recorded track data
    (
        ray_nav_df,
        ray_truth_df,
        ray_nav_cuda_df,
        helix_nav_df,
        helix_truth_df,
        helix_nav_cuda_df,
    ) = read_navigation_track_data(logging, datadir, det_name, p_min, p_max, args.cuda)

    # Plot track data
    plot_navigation_track_data(
        args,
        det_name,
        plot_factory,
        "ray",
        ray_truth_df,
        "truth",
        ray_nav_df,
        label_cpu,
        out_format,
    )

    plot_navigation_track_data(
        args,
        det_name,
        plot_factory,
        "helix",
        helix_truth_df,
        "truth",
        helix_nav_df,
        label_cpu,
        out_format,
    )

    if args.cuda:
        # Truth vs. Device
        plot_navigation_track_data(
            args,
            det_name,
            plot_factory,
            "ray",
            ray_truth_df,
            "truth",
            ray_nav_cuda_df,
            label_cuda,
            out_format,
        )

        plot_navigation_track_data(
            args,
            det_name,
            plot_factory,
            "helix",
            helix_truth_df,
            "truth",
            helix_nav_cuda_df,
            label_cuda,
            out_format,
        )

        # Host vs. Device
        plot_navigation_track_data(
            args,
            det_name,
            plot_factory,
            "ray",
            ray_nav_df,
            label_cpu,
            ray_nav_cuda_df,
            label_cuda,
            out_format,
        )

        plot_navigation_track_data(
            args,
            det_name,
            plot_factory,
            "helix",
            helix_nav_df,
            label_cpu,
            helix_nav_cuda_df,
            label_cuda,
            out_format,
        )


# ------------------------------------------------------------------------------

if __name__ == "__main__":
    __main__()

# ------------------------------------------------------------------------------
