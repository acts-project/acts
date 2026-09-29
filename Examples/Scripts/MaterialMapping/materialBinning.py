# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

"""Access current and legacy material binning at the JSON boundary."""


def binning(material):
    if "axis_specs" in material:
        return material["axis_specs"]
    return material.get("binUtility", {}).get("binningdata") or []


def binningDirection(axis):
    direction = axis.get("direction", axis.get("value"))
    # The scripts use the old names for comparisons. Cylinder phi and rphi
    # describe the same angular segmentation.
    return {
        "AxisX": "binX",
        "AxisY": "binY",
        "AxisZ": "binZ",
        "AxisR": "binR",
        "AxisPhi": "binPhi",
        "AxisRPhi": "binPhi",
        "binRPhi": "binPhi",
    }.get(direction, direction)


def configureBinning(material, config):
    """Copy bin counts by direction, retaining each surface's own ranges."""
    source = binning(config)
    target = binning(material)
    for index, axis in enumerate(target):
        direction = binningDirection(axis)
        if direction is None:
            configured = source[index]
        else:
            configured = next(
                candidate
                for candidate in source
                if binningDirection(candidate) == direction
            )
        if "bins" not in axis or "bins" not in configured:
            raise ValueError("Bin-count configuration requires equidistant axes")
        axis["bins"] = configured["bins"]
