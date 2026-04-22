# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import argparse
import logging
import sys
from datetime import datetime

# ------------------------------------------------------------------------------
# Options parsing
# ------------------------------------------------------------------------------

""" Parent parser that contains common options """


def common_options(prog_name=sys.argv[0]):

    parser = argparse.ArgumentParser(add_help=False, prog=prog_name)

    parser.add_argument(
        "--debug", "-d", help=("Enables debug logging"), action="store_true"
    )
    parser.add_argument("--logfile", help=("Write log in file"), default="", type=str)

    return parser


""" Parse common options from commandline """


def parse_common_options(args, prog_name=sys.argv[0]):

    # Set log level
    log_level = logging.INFO
    if args.debug:
        log_level = logging.DEBUG

    # Check logfile path
    if args.logfile != "":
        log_dir_name = os.path.dirname(args.logfile)

        if log_dir_name != "" and not os.path.isdir(log_dir_name):
            os.mkdir(log_dir_name, 0o755)

        # Write log in logfile
        logging.basicConfig(
            filename=args.logfile,
            format=("%(levelname)s (%(module)s): %(message)s"),
            level=log_level,
        )
    else:
        # Write log to terminal
        logging.basicConfig(
            format=("%(levelname)s (%(module)s): %(message)s"), level=log_level
        )

    logging.info(
        "\n--------------------------------------------------------\n"
        "Running "
        + prog_name
        + " "
        + str(datetime.now().strftime("%d/%m/%Y %H:%M"))
        + "\n--------------------------------------------------------\n"
    )

    return logging
