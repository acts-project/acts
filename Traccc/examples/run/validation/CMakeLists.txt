# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

traccc_add_executable( kalman_filter_validation "apps/kalman_filter_validation.cpp"
   LINK_LIBRARIES vecmem::core detray::io detray::detectors traccc::core
   traccc::io traccc::simulation traccc::performance traccc::options traccc_examples_common
)
