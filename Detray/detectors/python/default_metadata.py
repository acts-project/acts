# Detray library, part of the ACTS project (R&D line)
#
# (c) 2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

import detray

from detray.detectors import metadata, metadata_generator
from detray.detectors import Type, Algebra, Shape, Material, Accelerator
from detray.detectors import (
    add_silicon_tracker_defaults,
    add_calorimeter_defaults,
    add_telescope_detector_defaults,
    add_wire_chamber_defaults,
)

from itk_metadata import add_itk_types
from odd_metadata import add_odd_types

import argparse
import logging
import sys

# --------------------------------------------------------------------------
# Generate the default metadata type (can represent all detectors)
# --------------------------------------------------------------------------


def __main__():
    # Commandline options
    parser = argparse.ArgumentParser(prog=sys.argv[0])
    detray.detectors.add_logging_options(parser)
    detray.detectors.add_io_options(parser)

    args = parser.parse_args()
    detray.detectors.parse_logging_options(args)
    detray.detectors.parse_io_options(args)

    logger = logging.getLogger(__name__)

    # Collect the types required for a detector that can hold all detector types
    md = metadata("default")

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
    add_odd_types(md)
    add_itk_types(md)

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

    # Dump the metadata to header file
    if args.output:
        metadata_generator(md, args.output)
    else:
        metadata_generator(md)


if __name__ == "__main__":
    __main__()
