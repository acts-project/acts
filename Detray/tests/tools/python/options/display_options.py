# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import argparse
import os

# ------------------------------------------------------------------------------
# Options parsing
# ------------------------------------------------------------------------------


def display_options():
    """Parent parser that contains the detector display options."""

    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument(
        "--outdir",
        "-o",
        help=("Output directory for the svg files"),
        default="./plots/",
        type=str,
    )
    parser.add_argument(
        "--context", help=("Number of the geometry context"), default=0, type=int
    )
    parser.add_argument(
        "--r_axis", help=("Length of the radial axis [mm]"), default=1100.0, type=float
    )
    parser.add_argument(
        "--x_axis", help=("Length of the x axis [mm]"), default=None, type=float
    )
    parser.add_argument(
        "--y_axis", help=("Length of the y axis [mm]"), default=None, type=float
    )
    parser.add_argument(
        "--z_axis", help=("Length of the z axis [mm]"), default=3100.0, type=float
    )
    parser.add_argument("--font_size", help=("Font size"), default=35, type=int)
    parser.add_argument(
        "--search_window",
        help=("Size of the grid surface search window"),
        nargs=2,
        default=[1, 1],
        type=int,
    )
    parser.add_argument(
        "--volume_indices",
        help=("List of volumes that should be displayed by index"),
        nargs="+",
        default=[],
        type=int,
    )
    parser.add_argument(
        "--volume_names",
        help=("List of volumes that should be displayed by name"),
        nargs="+",
        default=[],
        type=str,
    )
    parser.add_argument(
        "--surfaces",
        help=("List of surfaces that should be displayed"),
        nargs="+",
        default=[],
        type=int,
    )
    parser.add_argument(
        "--hide_portals", help=("Hide portal surfaces"), action="store_true"
    )
    parser.add_argument(
        "--hide_passives", help=("Hide passive surfaces"), action="store_true"
    )
    parser.add_argument(
        "--hide_material", help=("Don't draw surface material"), action="store_true"
    )
    parser.add_argument(
        "--hide_eta_lines", help=("Hide eta lines"), action="store_true"
    )
    parser.add_argument("--show_info", help=("Show info boxes"), action="store_true")
    parser.add_argument(
        "--write_volume_graph",
        help=("Writes the volume graph to file"),
        action="store_true",
    )

    return parser


def parse_display_options(args, logging):
    """Parse the detector display options from the commandline."""

    # The axes default to the radial axis, but the gradient box is aligned with
    # the axis that was given explicitly
    if args.x_axis is None:
        args.x_axis = args.r_axis
        args.material_gradient_axis = args.z_axis
    else:
        args.material_gradient_axis = args.x_axis

    if args.y_axis is None:
        args.y_axis = args.r_axis

    # Create the output path
    os.makedirs(args.outdir, mode=0o755, exist_ok=True)

    return args.outdir
