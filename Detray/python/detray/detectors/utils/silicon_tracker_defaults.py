# Detray library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# --------------------------------------------------------------------------
# Function that will add types commonly needed for silicon trackers
# --------------------------------------------------------------------------

from ..impl import metadata
from ..impl import Shape, Material, Accelerator, GridBin

import logging

""" Types that are typically needed for silicon tracker detectors """


def add_silicon_tracker_defaults(
    metadata: metadata,
    use_mat_maps=False,
    use_homogeneous_mat=False,
    add_trapezoid=False,
):
    # Don't run more than once on given metadata (prevents spurious log entries)
    if metadata in add_silicon_tracker_defaults.clients:
        return

    add_silicon_tracker_defaults.clients.append(metadata)

    logger = logging.getLogger(__name__)
    logger.info("Define silicon tracker types:")

    # Cylindrical volume portals (barrel and endcap)
    logger.info("-> adding portal types")
    metadata.add_portal(Shape.CONCENTRIC_CYLINDER, type_id=2)
    metadata.add_portal(Shape.RING, type_id=3)

    # Acceleration struct for portals and passives
    metadata.add_accel_struct(
        Accelerator.BRUTE_FORCE, "portal", is_default=True, type_id=0
    )
    metadata.add_accel_struct(Accelerator.BRUTE_FORCE, "passive", type_id=0)

    # Barrel Detector
    logger.info("-> adding barrel section types")
    metadata.add_sensitive(Shape.RECTANGLE, type_id=0)
    metadata.add_accel_struct(
        Accelerator.CONCENTRIC_CYLINDER_GRID2D, "sensitive", type_id=1
    )
    if use_mat_maps:
        metadata.add_material(Material.CONCENTIRC_CYLINDER_MAP2D, type_id=0)

    # Endcap Detector
    logger.info("-> adding endcap section types")
    if add_trapezoid:
        metadata.add_sensitive(Shape.TRAPEZOID, type_id=1)
    metadata.add_accel_struct(Accelerator.DISC_GRID2D, "sensitive", type_id=2)
    if use_mat_maps:
        metadata.add_material(Material.DISC_MAP2D, type_id=1)

    # Material slabs can be used for both barrel and endcap surface shapes
    if use_homogeneous_mat:
        logger.info("-> requested homogeneous material types")
        metadata.add_material(Material.SLAB, type_id=3)

    # Volume accelerator for layered cylindrical detectors
    logger.info("-> adding detector volume acceleration structure")
    metadata.add_accel_struct(
        Accelerator.CYLINDER_GRID3D, "volume", grid_bin=GridBin.SINGLE, is_default=True
    )

    logger.info("Done")


add_silicon_tracker_defaults.clients = []
