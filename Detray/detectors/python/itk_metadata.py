# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import argparse

import detray
from detray.detectors import metadata, metadata_generator
from detray.detectors.configurations.itk import add_itk_types


def __main__():
    # Commandline options
    parser = argparse.ArgumentParser()
    detray.detectors.add_logging_options(parser)
    detray.detectors.add_io_options(parser)

    args = parser.parse_args()
    detray.detectors.parse_logging_options(args)
    detray.detectors.parse_io_options(args)

    md = metadata("itk")

    add_itk_types(md)

    # Dump the metadata to header file
    metadata_generator(
        md, output=args.output, format_header=args.format, output_file=args.output_file
    )


if __name__ == "__main__":
    __main__()
