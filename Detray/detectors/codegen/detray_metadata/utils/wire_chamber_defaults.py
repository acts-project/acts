# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# --------------------------------------------------------------------------
# Function that will add types commonly needed for wire chamber detectors
# --------------------------------------------------------------------------

from ..impl import metadata
from ..impl import Shape, Material, Accelerator

import logging

""" Types that are typically needed for wirechamber-like detectors """


def add_wire_chamber_defaults(
    md: metadata, use_mat_maps=False, use_homogeneous_mat=False
):
    # Don't run more than once on given metadata (prevents spurious log entries)
    if md in add_wire_chamber_defaults.clients:
        return

    add_wire_chamber_defaults.clients.append(md)

    logger = logging.getLogger(__name__)
    logger.info("Define wire chamber types:")

    # Sensitive shapes
    logger.info("-> adding sensitive types")
    md.add_sensitive(Shape.DRIFT_CELL, type_id=0)
    md.add_sensitive(Shape.STRAW_TUBE, type_id=1)

    # Cylindrical volume portals (barrel and endcap)
    logger.info("-> adding portal types")
    md.add_portal(Shape.CONCENTRIC_CYLINDER, type_id=2)
    md.add_portal(Shape.RING, type_id=3)

    # Acceleration struct for portals and passives
    md.add_accel_struct(Accelerator.BRUTE_FORCE, "portal", is_default=True)
    md.add_accel_struct(Accelerator.BRUTE_FORCE, "passive")

    # Surface acceleration structure for the wires
    md.add_accel_struct(Accelerator.CONCENTRIC_CYLINDER_GRID2D, "sensitive")

    if use_mat_maps:
        logger.info("-> requested material map types")
        # Map the material for the support structures
        md.add_material(Material.CONCENTIRC_CYLINDER_MAP2D)
        md.add_material(Material.DISC_MAP2D)
        # Experimental: map the all of the material above into 3D bins
        # md.add_material(Material.CYLINDER3D)

    if use_homogeneous_mat:
        logger.info("-> requested homogeneous material types")
        # Model the gas content
        md.add_material(Material.RAW)

    # Sensitive material
    md.add_material(Material.ROD)

    # Volume accelerator for layered cylindrical detectors
    logger.info("-> adding detector volume acceleration structure")
    md.add_accel_struct(Accelerator.CYLINDER_GRID3D, "volume", is_default=True)

    logger.info("Done")


add_wire_chamber_defaults.clients = []
