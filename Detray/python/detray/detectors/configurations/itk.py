# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from detray.detectors import metadata
from detray.detectors import Shape
from detray.detectors import (
    add_silicon_tracker_defaults,
)

import logging


def add_itk_types(md: metadata):
    logger = logging.getLogger(__name__)
    logger.info("Define types required by the ATLAS Inner Tracker (ITk):")

    # Strip stereo annulus shape
    logger.info("-> adding ITk strip detecor custom shape")
    md.add_sensitive(Shape.ANNULUS)

    # Add default types for silicon trackers (cylindrical detector shape)
    add_silicon_tracker_defaults(metadata=md, use_mat_maps=True)

    logger.info("Done")
