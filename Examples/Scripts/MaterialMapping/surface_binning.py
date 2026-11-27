# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


def surface_binning(material):
    """Return mutable axis descriptions from current or legacy surface JSON."""
    if "axis_specs" in material:
        return material["axis_specs"]
    return material.get("binUtility", {}).get("binningdata") or []


def axis_direction(axis, index, surface_type):
    """Normalize old direction names and resolve positional surface specs."""
    direction = axis.get("direction", axis.get("value"))
    if direction is not None:
        return direction.removeprefix("Axis").removeprefix("bin")
    directions = {
        "CylinderSurface": ("RPhi", "Z"),
        "DiscSurface": ("R", "Phi"),
        "PlaneSurface": ("X", "Y"),
    }
    return directions[surface_type][index]


def configure_surface(entry, config):
    """Apply bin counts by direction, including legacy cylinder phi names."""
    material = entry["value"]["material"]
    configured = config["value"]["material"]
    surface_type = entry["value"]["type"]
    material["mapMaterial"] = configured["mapMaterial"]
    material["mappingType"] = configured["mappingType"]

    def key(axis, index):
        direction = axis_direction(axis, index, surface_type)
        return "Phi" if direction == "RPhi" else direction

    counts = {
        key(axis, index): axis["bins"]
        for index, axis in enumerate(surface_binning(configured))
        if "bins" in axis
    }
    for index, axis in enumerate(surface_binning(material)):
        direction = key(axis, index)
        if "bins" in axis and direction in counts:
            axis["bins"] = counts[direction]
