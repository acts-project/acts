# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from detray_metadata import metadata
from detray_metadata import Shape
from detray_metadata import add_wire_chamber_defaults

import logging


def configure(md: metadata):
    logger = logging.getLogger(__name__)
    logger.info("Define types required by the wire chamber detector:")

    # Add default types for wire chambers
    add_wire_chamber_defaults(md, use_mat_maps=True, use_homogeneous_mat=True)

    logger.info("Done")
