# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# --------------------------------------------------------------------------
# Function that will add types commonly needed for calorimeters
# --------------------------------------------------------------------------

from ..impl import metadata
from ..impl import Shape, Material, Accelerator

import logging

""" Types that are typically needed for calorimeters """


def add_calorimeter_defaults(
    md: metadata, use_mat_maps=False, use_homogeneous_mat=False
):
    # Don't run more than once on given metadata (prevents spurious log entries)
    if md in add_calorimeter_defaults.clients:
        return

    add_calorimeter_defaults.clients.append(md)

    logger = logging.getLogger(__name__)
    logger.info("Define calorimeter types:")

    logger.info("-> adding sensitive types")
    md.add_sensitive(Shape.RECTANGLE, type_id=0)
    md.add_sensitive(Shape.TRAPEZOID, type_id=1)

    logger.info("-> adding portal types")
    md.add_portal(Shape.CONCENTRIC_CYLINDER, type_id=2)
    md.add_portal(Shape.RING, type_id=3)

    # Acceleration Struct for portals and passives
    md.add_accel_struct(Accelerator.BRUTE_FORCE, "portal", is_default=True)
    md.add_accel_struct(Accelerator.BRUTE_FORCE, "passive")

    if use_mat_maps:
        logger.info("-> requested material map types")
        md.add_material(Material.CYLINDER_MAP3D)
    if use_homogeneous_mat:
        logger.info("-> requested homogeneous material types")
        md.add_material(Material.RAW)

    # Add acceleration structures (e.g. Frustum navigation) in the future...

    logger.info("Done")


add_calorimeter_defaults.clients = []
