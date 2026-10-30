# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from detray.detectors import metadata
from detray.detectors import Shape, Accelerator, GridBin, GridSerializer
from detray.detectors import add_silicon_tracker_defaults

import logging


def add_toy_types(md: metadata):
    logger = logging.getLogger(__name__)
    logger.info("Define types required by the toy detector:")

    # Add default types for silicon trackers (cylindrical detector shape)
    add_silicon_tracker_defaults(
        metadata=md, use_homogeneous_mat=True, use_mat_maps=True, add_trapezoid=True
    )

    # Add surface grids with static bin capacity
    toy_grid_bin = GridBin.STATIC
    toy_grid_bin.param["capacity"] = 1

    md.add_accel_struct(
        Accelerator.CONCENTRIC_CYLINDER_GRID2D,
        "sensitive",
        grid_bin=toy_grid_bin,
        type_id=1,
    )

    md.add_accel_struct(
        Accelerator.DISC_GRID2D, "sensitive", grid_bin=toy_grid_bin, type_id=2
    )

    logger.info("Done")
