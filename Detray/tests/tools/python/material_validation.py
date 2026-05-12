# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# detray includes
from impl import plot_material_scan as mat_plotter
from impl import read_material_data
from plotting import pyplot_factory as plt_factory
from options import (
    common_options,
    detector_io_options,
    uniform_track_generator_options,
    propagation_options,
    plotting_options,
)
from options import (
    parse_common_options,
    parse_detector_io_options,
    parse_plotting_options,
)
from utils import read_detector_name
from utils import add_track_generator_args, add_propagation_args, add_detector_io_args

# python includes
import argparse
import json
import os
import subprocess
import sys


def __main__():

    # ---------------------------------------------------------------arg parsing

    descr = "Detray Material Validation"

    # Define options
    parent_parsers = [
        common_options(descr),
        detector_io_options(),
        uniform_track_generator_options(),
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
        default="./validation_data/material",
        type=str,
    )
    parser.add_argument(
        "--material_tol",
        "-mt",
        help=("Tolerance for material comparisons [%]"),
        default=1,
        type=float,
    )
    parser.add_argument(
        "--overlaps_tol",
        "-ot",
        help=("Tolerance for considering surfaces to be overlapping [mm]"),
        default=0.0001,
        type=float,
    )
    parser.add_argument(
        "--cuda",
        help=("Run the CUDA material validation."),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--sycl",
        help=("Run the SYCL material validation."),
        action="store_true",
        default=False,
    )

    args = parser.parse_args()

    logging = parse_common_options(args, descr)
    parse_detector_io_options(args, logging)
    _, out_dir, out_format = parse_plotting_options(args, logging)

    # IO path for data files
    datadir = args.datadir.strip("/")

    # Check bin path
    bindir = args.bindir.strip("/")
    cpu_validation = bindir + "/detray_material_validation"
    cuda_validation = bindir + "/detray_material_validation_cuda"

    if not os.path.isdir(bindir) or not os.path.isfile(cpu_validation):
        logging.error(f"Material validation binaries were not found! ({args.bindir})")
        sys.exit(1)

    # -----------------------------------------------------------------------run

    # Pass on the options for the validation tools
    args_list = [
        "--data_dir",
        datadir,
        "--material_tol",
        str(args.material_tol),
        "--overlaps_tol",
        str(args.overlaps_tol),
    ]

    # Add parsed options to argument list
    add_detector_io_args(args_list, args)
    add_track_generator_args(args_list, args)
    add_propagation_args(args_list, args)

    logging.debug(args_list)

    if "--material_file" not in args_list:
        logging.error(
            "Detector material is required! Please add it using the '--material_file' option"
        )
        sys.exit(1)

    # Run the host validation and produce the truth data
    logging.debug("Running CPU material validation")
    subprocess.run([cpu_validation] + args_list)

    # Run the device validation (if it has been built)
    if args.cuda and os.path.isfile(cuda_validation):
        logging.debug("Running CUDA material validation")
        subprocess.run([cuda_validation] + args_list)

    elif args.cuda:
        logging.error("Could not find CUDA material validation executable")

    if args.sycl:
        logging.error("SYCL material validation is not implemented")

    # ----------------------------------------------------------------------plot

    logging.info("Generating data plots...\n")

    det_name = read_detector_name(args.geometry_file, logging)
    logging.debug("Detector: " + det_name)

    # Check the data path (should have been created when running the validation)
    if not os.path.isdir(datadir):
        logging.error(f"Data directory was not found! ({args.datadir})")
        sys.exit(1)

    df_scan, df_cpu, df_cuda = read_material_data(datadir, logging, det_name, args.cuda)

    plot_factory = plt_factory(out_dir, logging)

    # The histograms are not re-weighted (if the rays are not evenly distributed
    # the material in some bins might be artificially high)!
    mat_plotter.X0_vs_eta_phi(
        df_scan, "material_scan", det_name, plot_factory, out_format
    )
    mat_plotter.L0_vs_eta_phi(
        df_scan, "material_scan", det_name, plot_factory, out_format
    )
    mat_plotter.X0_vs_eta(df_scan, "material_scan", det_name, plot_factory, out_format)
    mat_plotter.L0_vs_eta(df_scan, "material_scan", det_name, plot_factory, out_format)

    # Navigation material Traces
    # CPU
    mat_plotter.X0_vs_eta_phi(
        df_cpu, "cpu_material_trace", det_name, plot_factory, out_format
    )
    mat_plotter.L0_vs_eta_phi(
        df_cpu, "cpu_material_trace", det_name, plot_factory, out_format
    )
    mat_plotter.X0_vs_eta(
        df_cpu, "cpu_material_trace", det_name, plot_factory, out_format
    )
    mat_plotter.L0_vs_eta(
        df_cpu, "cpu_material_trace", det_name, plot_factory, out_format
    )

    # Comparison between scan and navigator trace in sX0
    mat_plotter.compare_mat(
        df_scan, df_cpu, "cpu_material", det_name, plot_factory, out_format
    )

    # CUDA
    if args.cuda:
        mat_plotter.X0_vs_eta_phi(
            df_cuda, "cuda_material_trace", det_name, plot_factory, out_format
        )
        mat_plotter.L0_vs_eta_phi(
            df_cuda, "cuda_material_trace", det_name, plot_factory, out_format
        )
        mat_plotter.X0_vs_eta(
            df_cuda, "cuda_material_trace", det_name, plot_factory, out_format
        )
        mat_plotter.L0_vs_eta(
            df_cuda, "cuda_material_trace", det_name, plot_factory, out_format
        )
        mat_plotter.compare_mat(
            df_scan, df_cuda, "cuda_material", det_name, plot_factory, out_format
        )


# ------------------------------------------------------------------------------

if __name__ == "__main__":
    __main__()

# ------------------------------------------------------------------------------
