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
