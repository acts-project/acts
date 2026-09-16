# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# detray includes
from options import (
    common_options,
    detector_writer_options,
    toy_detector_options,
)
from options import (
    parse_common_options,
    fill_toy_detector_config,
    fill_writer_config,
)

# detray python bindings
import detray.core
import detray.io
import detray.tests

# python includes
import argparse


def run_generate_toy_detector(args):

    # Build the geometry
    toy_cfg = fill_toy_detector_config(args, detray.tests.ToyDetectorConfig())
    det, names = detray.tests.buildToyDetector(
        detray.core.HostMemoryResource(), toy_cfg
    )

    # Write to file
    writer_cfg = fill_writer_config(args, detray.io.DetectorWriterConfig())

    # Make sure material is written to file, if it was requested
    writer_cfg.writeMaterial = (
        args.homogeneous_material or args.material_maps or args.write_material
    )

    detray.io.writeDetector(det, names, writer_cfg)


def __main__():

    # ---------------------------------------------------------------arg parsing

    descr = "Detray Toy Detector Generation"

    # Define options
    parent_parsers = [
        common_options(descr),
        toy_detector_options(),
        detector_writer_options(),
    ]

    parser = argparse.ArgumentParser(description=descr, parents=parent_parsers)

    parser.add_argument(
        "--write_volume_graph",
        help=("Write the volume graph to file"),
        action="store_true",
    )

    args = parser.parse_args()

    logging = parse_common_options(args, descr)

    # -----------------------------------------------------------------------run

    logging.debug("Generating the toy detector")
    run_generate_toy_detector(args)

    logging.info(f"Wrote the detector files to '{args.outdir}'")

    # General options
    if args.write_volume_graph:
        raise NotImplementedError("Writing of volume graph not implemented")


# ------------------------------------------------------------------------------

if __name__ == "__main__":
    __main__()

# ------------------------------------------------------------------------------
