# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import argparse
import os
import sys

# ------------------------------------------------------------------------------
# Options parsing
# ------------------------------------------------------------------------------

""" Parent detector reader options that contain common options """


def detector_io_options():

    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument(
        "--geometry_file",
        "-geo",
        help=("Detector geometry file"),
        default="",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--grid_file",
        "-grid",
        help=("Detector surface grids file"),
        default="",
        type=str,
    )
    parser.add_argument(
        "--material_file", "-mat", help=("Detector material file"), default="", type=str
    )

    return parser


""" Parse detector reader options from commandline """


def parse_detector_io_options(args, logging):

    # Check detector files
    if not os.path.isfile(args.geometry_file):
        logging.error(f"Detector geometry file does not exist! ({args.geometry_file})")
        sys.exit(1)

    if args.grid_file and not os.path.isfile(args.grid_file):
        logging.error(f"Detector grid file does not exist! ({args.material_file})")
        sys.exit(1)

    if args.material_file and not os.path.isfile(args.material_file):
        logging.error(f"Detector material file does not exist! ({args.material_file})")
        sys.exit(1)


""" Fill a detector reader config from the parsed commandline options """


def fill_reader_config(args, config):

    if not args.geometry_file:
        raise ValueError("Please specify a geometry input file!")

    config.addFile(args.geometry_file)

    if args.material_file:
        config.addFile(args.material_file)

    if args.grid_file:
        config.addFile(args.grid_file)

    return config


def detector_writer_options():
    """Parent parser that contains the detector writer options."""

    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument(
        "--outdir",
        "-o",
        help=("Output directory for the detector files"),
        default="./toy_detector/",
        type=str,
    )
    parser.add_argument(
        "--write_material", help=("Toggle material output"), action="store_true"
    )
    parser.add_argument(
        "--write_grids", help=("Toggle grid output"), action="store_true"
    )
    parser.add_argument(
        "--replace_files",
        help=("Whether to replace existing files"),
        action="store_true",
    )
    parser.add_argument(
        "--compactify_json", help=("Not implemented"), action="store_true"
    )

    return parser


def fill_writer_config(args, config):
    """Fill a detector writer config from the parsed commandline options."""

    config.path = args.outdir
    config.compactifyJson = args.compactify_json
    config.writeMaterial = args.write_material
    config.writeGrids = args.write_grids
    config.replaceFiles = args.replace_files

    return config
