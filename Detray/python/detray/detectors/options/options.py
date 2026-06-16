# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import argparse
import logging
import os
import sys
from datetime import datetime

# ------------------------------------------------------------------------------
# Options parsing
# ------------------------------------------------------------------------------

""" Parent parser that contains logging options """


def add_logging_options(parser):

    parser.add_argument(
        "-v",
        "--info",
        action="store_const",
        dest="loglevel",
        const=logging.INFO,
        help=("increase log level (INFO)"),
    )
    parser.add_argument(
        "-vv",
        "--debug",
        action="store_const",
        dest="loglevel",
        const=logging.DEBUG,
        help=("increase log level (DEBUG)"),
    )

    return parser


""" Parse logging options from commandline and set up logging service"""


def parse_logging_options(args, prog_name=sys.argv[0]):
    # Configure log level
    log_level = logging.WARNING
    if args.loglevel:
        log_level = args.loglevel

    # Write log to terminal
    logging.basicConfig(
        format=("%(levelname)s (%(module)s): %(message)s"), level=log_level
    )

    logging.info(
        "\ndetray - "
        + str(datetime.now().strftime("%d/%m/%Y %H:%M"))
        + ': Running detector type generator "'
        + os.path.basename(prog_name)
        + '"\n'
    )


""" Parent parser that contains IO options """


def add_io_options(parser):

    parser.add_argument(
        "-o",
        "--output",
        help=("Metadata output directory"),
        type=str,
    )
    parser.add_argument(
        "-f",
        "--format",
        help=("Format the header file"),
        default=True,
        action=argparse.BooleanOptionalAction,
    )

    return parser


""" Parse IO options from commandline"""


def parse_io_options(args):
    # Check if the output directory exists
    if args.output and not os.path.isdir(args.output):
        print(f"Output directory does not exist! ({args.output})")
        sys.exit(1)
