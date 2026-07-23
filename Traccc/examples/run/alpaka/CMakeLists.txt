# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2023-2025 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Project include(s).
include(traccc-alpaka-functions)

set(EXTRA_LIBS)

if(alpaka_ACC_GPU_CUDA_ENABLE)
    list(APPEND EXTRA_LIBS vecmem::cuda)
elseif(alpaka_ACC_GPU_HIP_ENABLE)
    find_package(HIPToolkit REQUIRED)
    list(APPEND EXTRA_LIBS HIP::hiprt vecmem::hip)
elseif(alpaka_ACC_SYCL_ENABLE)
    list(APPEND EXTRA_LIBS vecmem::sycl)
endif()

set(LIBRARIES
    vecmem::core
    traccc::io
    traccc::performance
    traccc::core
    traccc::device_common
    traccc::alpaka
    traccc::options
    traccc_examples_common
    ${EXTRA_LIBS}
)
set(DETRAY detray::io detray::detectors)

traccc_add_executable( seq_example_alpaka "apps/seq_example_alpaka.cpp"
    LINK_LIBRARIES ${LIBRARIES} ${DETRAY}
)
traccc_add_executable( seeding_example_alpaka "apps/seeding_example_alpaka.cpp"
    LINK_LIBRARIES ${LIBRARIES}
)

#
# Set up the "throughput applications".
#
traccc_add_library( traccc_examples_alpaka examples_alpaka SHARED
   "include/traccc/examples/alpaka/full_chain_algorithm.hpp"
   "src/full_chain_algorithm.cpp"
)
target_link_libraries(
    traccc_examples_alpaka
    PUBLIC
        vecmem::core
        detray::core
        detray::detectors
        traccc::core
        traccc::device_common
        traccc::alpaka
        traccc_examples_common
        ${EXTRA_LIBS}
)

traccc_add_executable( throughput_st_alpaka "apps/throughput_st.cpp"
   LINK_LIBRARIES ${LIBRARIES} ${DETRAY} traccc_examples_alpaka
)

traccc_add_executable( throughput_mt_alpaka "apps/throughput_mt.cpp"
   LINK_LIBRARIES TBB::tbb ${LIBRARIES} ${DETRAY} traccc_examples_alpaka
)

add_test(
    NAME throughput_st_alpaka
    COMMAND
        traccc_throughput_st_alpaka --finding-run-smoother=none
        --read-bfield-from-file --cold-run-events=1 --processed-events=1
)

add_test(
    NAME throughput_mt_alpaka
    COMMAND
        traccc_throughput_mt_alpaka --finding-run-smoother=none
        --read-bfield-from-file --cold-run-events=1 --processed-events=1
)
