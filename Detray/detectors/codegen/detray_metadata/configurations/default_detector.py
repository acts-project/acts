# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from detray_metadata import metadata
from detray_metadata import Type, Algebra, Shape, Material, Accelerator
from detray_metadata import (
    add_silicon_tracker_defaults,
    add_calorimeter_defaults,
    add_telescope_detector_defaults,
    add_wire_chamber_defaults,
)

from .odd_detector import configure as configure_odd

import logging


def configure(md: metadata):
    logger = logging.getLogger(__name__)

    # Collect the types required for a detector that can hold all detector types

    # Specify a particular algebra plugin (otherwise left as template param.)
    # md.set_algebra_plugin(Algebra.ARRAY, Type.SINGLE)

    # Make sure all of the defaults are added
    add_telescope_detector_defaults(md, use_mat_maps=True, use_homogeneous_mat=True)
    add_wire_chamber_defaults(md, use_mat_maps=True, use_homogeneous_mat=True)
    add_calorimeter_defaults(md, use_mat_maps=True, use_homogeneous_mat=True)
    add_silicon_tracker_defaults(
        md, use_mat_maps=True, use_homogeneous_mat=True, add_trapezoid=True
    )

    # Make sure all detectors can be represented by this metadata
    configure_odd(md)
    # Stereo annulus strip modules (e.g. ATLAS ITk)
    md.add_sensitive(Shape.ANNULUS)

    # Extra passive surface shapes (cylinder with two intersection solutions)
    md.add_passive(Shape.CYLINDER2D)
    md.add_passive(Shape.RING)

    # Add an acceleration struct and material map for generic 2D cylinders
    md.add_accel_struct(Accelerator.CYLINDER_GRID2D)
    md.add_material(Material.CYLINDER_MAP2D)

    #
    # Add some special types (e.g. for detector R&D)

    # Special surface types
    # md.add_sensitive(Shape.UNMASKED)

    # Add more material types for testing
    md.add_material(Material.ANNULUS_MAP2D)
    md.add_material(Material.TRAPEZOID_MAP2D)
    md.add_material(Material.CUBOID_MAP3D)
    md.add_material(Material.CYLINDER_MAP2D)

    # Make sure a default acceleration struct is chosen that can be used in all
    # detector types
    logger.info(
        "Overwrite default acceleration structures to use most generic types as defaults"
    )
    md.set_default_accel_struct(Accelerator.BRUTE_FORCE, "portal")
    md.set_default_accel_struct(Accelerator.CYLINDER_GRID3D, "volume")
