# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import argparse

# ------------------------------------------------------------------------------
# Options parsing
# ------------------------------------------------------------------------------


def toy_detector_options():
    """Parent parser that contains the toy detector options."""

    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument(
        "--barrel_layers",
        help=("Number of barrel layers [0-4]"),
        default=4,
        type=int,
    )
    parser.add_argument(
        "--endcap_layers",
        help=("Number of endcap layers on either side [0-7]"),
        default=3,
        type=int,
    )
    parser.add_argument(
        "--homogeneous_material",
        help=("Generate homogeneous material description (default)"),
        action="store_true",
    )
    parser.add_argument(
        "--material_maps", help=("Generate material maps"), action="store_true"
    )

    return parser


def fill_toy_detector_config(args, config):
    """Fill a toy detector config from the parsed commandline options."""

    if args.homogeneous_material and args.material_maps:
        raise ValueError("Please specify only one material description")

    config.nBarrelLayers = args.barrel_layers
    config.nEndcapLayers = args.endcap_layers

    if args.homogeneous_material:
        config.useMaterialMaps = False
    if args.material_maps:
        config.useMaterialMaps = True

    return config
