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

""" Parent parser that contains plotting options """


def plotting_options():

    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument(
        "--inputdir",
        "-i",
        help=("Directory containing input data files."),
        default="./",
        type=str,
    )
    parser.add_argument(
        "--outdir",
        "-o",
        help=("Output directory for plots."),
        default="./plots/",
        type=str,
    )
    parser.add_argument(
        "--output_format",
        "-of",
        help=("Format of the plot files (svg|png|pdf)."),
        default="pdf",
        type=str,
    )

    return parser


""" Parse plotting options from commandline """


def parse_plotting_options(args, logging):

    # Check input path
    if args.inputdir and not os.path.isdir(args.inputdir):
        logging.error(f"Plot data director does not exist! ({args.inputdir})")
        sys.exit(1)

    # Check output path
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir, 0o755)

    if args.output_format not in ["svg", "png", "pdf"]:
        logging.error(f"Unknown output file format: {out_format}")
        sys.exit(1)

    return args.inputdir, args.outdir, args.output_format
