# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from detray.detectors import metadata
from detray.detectors import Shape
from detray.detectors import add_silicon_tracker_defaults

import logging


def add_odd_types(md: metadata):
    logger = logging.getLogger(__name__)
    logger.info("Define types required by the ACTS Open Data Detector (ODD):")

    # Add default types for silicon trackers (cylindrical detector shape)
    add_silicon_tracker_defaults(
        metadata=md, use_homogeneous_mat=True, use_mat_maps=True, add_trapezoid=True
    )

    # Add passive material surfaces (intersectable for cosmics)
    # md.add_passive(Shape.CYLINDER2D)
    # md.add_passive(Shape.RING)

    logger.info("Done")
