# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Generate a detector metadata header from a packaged or external configuration."""

import argparse
from importlib import import_module
from pathlib import Path
import runpy

from detray.detectors import (
    add_io_options,
    add_logging_options,
    metadata,
    metadata_generator,
    parse_io_options,
    parse_logging_options,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--detector",
        required=True,
        metavar="MODULE.FUNCTION or FILE:FUNCTION",
        help="Packaged configuration callable or external Python file and callable; the module/file stem names the detector",
    )
    add_logging_options(parser)
    add_io_options(parser)
    args = parser.parse_args()
    parse_logging_options(args)
    parse_io_options(args)

    try:
        if ":" in args.detector:
            filename, function_name = args.detector.rsplit(":", 1)
            detector_name = Path(filename).stem
            configure = runpy.run_path(filename)[function_name]
        else:
            module_name, function_name = args.detector.rsplit(".", 1)
            module = import_module(f"detray.detectors.configurations.{module_name}")
            detector_name = module_name.rsplit(".", 1)[-1]
            configure = getattr(module, function_name)
        if not callable(configure):
            raise ValueError(f"{args.detector} is not callable")
    except (ValueError, ImportError, AttributeError, KeyError, OSError) as error:
        parser.error(f"Invalid detector configuration {args.detector!r}: {error}")

    md = metadata(detector_name)
    configure(md)
    metadata_generator(
        md, output=args.output, format_header=args.format, output_file=args.output_file
    )


if __name__ == "__main__":
    main()
