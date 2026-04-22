# Detray library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

import detray

from detray.detectors import metadata, metadata_generator
from detray.detectors import Shape, Accelerator, GridBin, GridSerializer
from detray.detectors import add_silicon_tracker_defaults

import argparse
import logging
import sys

# --------------------------------------------------------------------------
# Generate the toy detector metadata type
# --------------------------------------------------------------------------

""" Add all types needed to describe the detray toy detector """


def add_toy_types(md: metadata):
    logger = logging.getLogger(__name__)
    logger.info("Define types required by the toy detector:")

    # Add default types for silicon trackers (cylindrical detector shape)
    add_silicon_tracker_defaults(
        metadata=md, use_homogeneous_mat=True, use_mat_maps=True, add_trapezoid=True
    )

    # Add surface grids with static bin capacity
    toy_grid_bin = GridBin.STATIC
    toy_grid_bin.param["capacity"] = 1

    md.add_accel_struct(
        Accelerator.CONCENTRIC_CYLINDER_GRID2D,
        "sensitive",
        grid_bin=toy_grid_bin,
        type_id=1,
    )

    md.add_accel_struct(
        Accelerator.DISC_GRID2D, "sensitive", grid_bin=toy_grid_bin, type_id=2
    )

    logger.info("Done")


def __main__():
    # Commandline options
    parser = argparse.ArgumentParser(prog=sys.argv[0])
    detray.detectors.add_logging_options(parser)
    detray.detectors.add_io_options(parser)

    args = parser.parse_args()
    detray.detectors.parse_logging_options(args)
    detray.detectors.parse_io_options(args)

    md = metadata("toy")

    add_toy_types(md)

    # Dump the metadata to header file
    if args.output:
        metadata_generator(md, args.output)
    else:
        metadata_generator(md)


if __name__ == "__main__":
    __main__()
