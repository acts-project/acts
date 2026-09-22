# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# External(s).
find_package(HIPToolkit REQUIRED)

#
# Set up the "throughput applications".
#
traccc_add_library( traccc_examples_hip examples_hip SHARED
   "include/traccc/examples/hip/full_chain_algorithm.hpp"
   "src/full_chain_algorithm.cpp"
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
