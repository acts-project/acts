# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# External(s).
find_package(HIPToolkit REQUIRED)

#
# Set up the "throughput applications".
#
traccc_add_library( traccc_examples_hip examples_hip SHARED
   "include/traccc/examples/hip/full_chain_algorithm.hpp"
   "src/full_chain_algorithm.cpp"
   "include/traccc/examples/hip/tbb_await.hpp"
   "src/tbb_await.cpp"
)
target_link_libraries(
    traccc_examples_hip
    PUBLIC
        HIP::hiprt
        vecmem::core
        vecmem::hip
        detray::core
        detray::detectors
        traccc::core
        traccc::device_common
        traccc::hip_utils
        traccc::hip
        traccc::examples_common
    PRIVATE TBB::tbb
)

traccc_add_executable( throughput_st_hip "apps/throughput_st.cpp"
   LINK_LIBRARIES vecmem::core vecmem::hip
                  detray::detectors detray::io traccc::io traccc::performance
                  traccc::core traccc::device_common traccc::hip_utils traccc::hip
                  traccc::options traccc::examples_hip
)

traccc_add_executable( throughput_mt_hip "apps/throughput_mt.cpp"
   LINK_LIBRARIES TBB::tbb vecmem::core vecmem::hip
                  detray::detectors detray::io traccc::io traccc::performance
                  traccc::core traccc::device_common traccc::hip_utils traccc::hip
                  traccc::options traccc::examples_hip
)
