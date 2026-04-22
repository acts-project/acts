# Detray library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

import detray

from detray.detectors import metadata, metadata_generator
from detray.detectors import Shape
from detray.detectors import add_silicon_tracker_defaults

import argparse
import logging
import sys

# --------------------------------------------------------------------------
# Generate the Open Data Detector metadata type
# --------------------------------------------------------------------------

""" Add all types needed to describe the ACTS Open Data Detector (ODD) """


def add_odd_types(md: metadata):
    logger = logging.getLogger(__name__)
    logger.info("Define types required by the ACTS Open Data Detector (ODD):")

    # Add default types for silicon trackers (cylindrical detector shape)
    add_silicon_tracker_defaults(
        metadata=md, use_homogeneous_mat=True, use_mat_maps=True, add_trapezoid=True
    )

    # Add passive material surfaces (intersectable for cosmics)
    # md.add_passive(Shape.CYLINDER2D)
    # md.add_passive(Shape.RING)

    logger.info("Done")


def __main__():
    # Commandline options
    parser = argparse.ArgumentParser(prog=sys.argv[0])
    detray.detectors.add_logging_options(parser)
    detray.detectors.add_io_options(parser)

    args = parser.parse_args()
    detray.detectors.parse_logging_options(args)
    detray.detectors.parse_io_options(args)

    md = metadata("odd")

    add_odd_types(md)

    # Dump the metadata to header file
    if args.output:
        metadata_generator(md, args.output)
    else:
        metadata_generator(md)


if __name__ == "__main__":
    __main__()
