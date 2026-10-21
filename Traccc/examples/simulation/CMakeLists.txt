# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Simulation with the default metadata and json inputs
traccc_add_executable(simulate "simulate.cpp"
    LINK_LIBRARIES vecmem::core traccc::io traccc::core
    traccc::options traccc::simulation detray::io covfie::core Boost::filesystem
)

# Simulation with the detectors pre-built in detray
traccc_add_executable(simulate_telescope "simulate_telescope.cpp"
    LINK_LIBRARIES vecmem::core traccc::io traccc::core
    traccc::options traccc::simulation detray::core detray::io detray::test_common covfie::core Boost::filesystem
)

traccc_add_executable(simulate_toy_detector "simulate_toy_detector.cpp"
    LINK_LIBRARIES vecmem::core traccc::io traccc::core
    traccc::options traccc::simulation detray::io detray::test_common covfie::core Boost::filesystem
)

traccc_add_executable(simulate_wire_chamber "simulate_wire_chamber.cpp"
    LINK_LIBRARIES vecmem::core traccc::io traccc::core
    traccc::options traccc::simulation detray::io detray::test_common covfie::core Boost::filesystem
)
