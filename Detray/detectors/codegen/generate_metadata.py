# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Generate a detector metadata header from a packaged or external configuration.

A configuration is a module named ``<name>_detector`` that provides a
``configure(md)`` function. It produces ``<name>_metadata.hpp`` declaring
``detray::<name>_metadata``.
"""

import argparse
from importlib import import_module
from pathlib import Path
import runpy

from detray_metadata import (
    add_io_options,
    add_logging_options,
    metadata,
    metadata_generator,
    parse_io_options,
    parse_logging_options,
)


def detector_name(stem: str) -> str:
    """Map a configuration module name to the detector name."""
    name = stem.removesuffix("_detector")
    if not name or not name.isidentifier():
        raise ValueError(f"cannot derive a detector name from {stem!r}")
    return name


def load_configuration(selector: str):
    """Resolve NAME or FILE.py:FUNCTION to (detector name, callable)."""
    if ":" in selector:
        filename, function_name = selector.rsplit(":", 1)
        name = detector_name(Path(filename).stem)
        configure = runpy.run_path(filename)[function_name]
    else:
        name = selector
        module = import_module(f"detray_metadata.configurations.{name}_detector")
        configure = module.configure
    if not callable(configure):
        raise ValueError(f"{selector} is not callable")
    return name, configure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--detector",
        required=True,
        metavar="NAME or FILE:FUNCTION",
        help="Shipped detector name (e.g. 'odd') or an external Python file and "
        "callable; <name>_detector.py names the detector <name>",
    )
    add_logging_options(parser)
    add_io_options(parser)
    args = parser.parse_args()
    parse_logging_options(args)
    parse_io_options(args)

    try:
        name, configure = load_configuration(args.detector)
    except (ValueError, ImportError, AttributeError, KeyError, OSError) as error:
        parser.error(f"Invalid detector configuration {args.detector!r}: {error}")

    md = metadata(name)
    configure(md)
    metadata_generator(
        md, output=args.output, format_header=args.format, output_file=args.output_file
    )


if __name__ == "__main__":
    main()
