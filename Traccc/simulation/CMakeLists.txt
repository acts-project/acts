# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Set up the "build" of the traccc::io library.
traccc_add_library( traccc_simulation simulation TYPE INTERFACE
  # Public headers
  "include/traccc/simulation/event_generators.hpp"
  "include/traccc/simulation/measurement_smearer.hpp"
  "include/traccc/simulation/simulator.hpp"
  "include/traccc/simulation/smearing_writer.hpp"
)
target_link_libraries(
    traccc_simulation
    INTERFACE
        traccc::core
        traccc::io
        detray::core
        detray::io
        detray::test_common
        detray::test_utils
)
