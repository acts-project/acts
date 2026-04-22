# Detray library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# --------------------------------------------------------------------------
# Function that will add types commonly needed for telescope detectors
# --------------------------------------------------------------------------

from ..impl import metadata
from ..impl import Shape, Material, Accelerator

import logging

""" Types that are typically needed for telescope detectors """


def add_telescope_detector_defaults(
    md: metadata, use_mat_maps=False, use_homogeneous_mat=False
):
    # Don't run more than once on given metadata (prevents spurious log entries)
    if md in add_telescope_detector_defaults.clients:
        return

    add_telescope_detector_defaults.clients.append(md)

    logger = logging.getLogger(__name__)
    logger.info("Define telescope detector types:")

    # Sensitive and portal shapes
    md.add_sensitive(Shape.RECTANGLE, type_id=0)
    md.add_portal(Shape.RECTANGLE, type_id=1)

    # Acceleration struct for portals and passives
    md.add_accel_struct(Accelerator.BRUTE_FORCE, "portal", is_default=True)
    md.add_accel_struct(Accelerator.BRUTE_FORCE, "passive")
    md.add_accel_struct(Accelerator.BRUTE_FORCE, "sensitive")

    # Map the material for the support structures onto the portals
    if use_mat_maps:
        logger.info("-> requested material map types")
        md.add_material(Material.RECTANGLE_MAP2D)

    # Sensitive material
    if use_homogeneous_mat:
        logger.info("-> requested homogeneous material types")
        md.add_material(Material.SLAB)

    # Acceleration struct for telescope detector volumes
    logger.info("-> adding detector volume acceleration structure")
    md.add_accel_struct(Accelerator.BRUTE_FORCE, "volume", is_default=True)

    logger.info("Done")


add_telescope_detector_defaults.clients = []
