# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import detray

from detray.detectors import metadata, metadata_generator
from detray.detectors import Shape
from detray.detectors import (
    add_silicon_tracker_defaults,
)

import argparse
import logging
import sys

# --------------------------------------------------------------------------
# Generate the ITk metadata type
# --------------------------------------------------------------------------

""" Add all types needed to describe the ATLAS ITk detector """


def add_itk_types(md: metadata):
    logger = logging.getLogger(__name__)
    logger.info("Define types required by the ATLAS Inner Tracker (ITk):")

    # Strip stereo annulus shape
    logger.info("-> adding ITk strip detecor custom shape")
    md.add_sensitive(Shape.ANNULUS)

    # Add default types for silicon trackers (cylindrical detector shape)
    add_silicon_tracker_defaults(metadata=md, use_mat_maps=True)

    logger.info("Done")


def __main__():
    # Commandline options
    parser = argparse.ArgumentParser(prog=sys.argv[0])
    detray.detectors.add_logging_options(parser)
    detray.detectors.add_io_options(parser)

    args = parser.parse_args()
    detray.detectors.parse_logging_options(args)
    detray.detectors.parse_io_options(args)

    md = metadata("itk")

    add_itk_types(md)

    # Dump the metadata to header file
    metadata_generator(md, output=args.output, format_header=args.format)


if __name__ == "__main__":
    __main__()
