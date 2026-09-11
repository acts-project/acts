# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2021-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# We need the SYCL language in this directory.
include(vecmem-check-language)
vecmem_check_language( SYCL )
enable_language(SYCL)

# Project include(s).
include(traccc-compiler-options-cpp)
include(traccc-compiler-options-sycl)

# Set up the build of the traccc::sycl library.
traccc_add_library( traccc_sycl sycl TYPE SHARED
  # Clusterization algorithm(s).
  "include/traccc/sycl/clusterization/clusterization_algorithm.hpp"
  "src/clusterization/clusterization_algorithm.sycl"
  "include/traccc/sycl/clusterization/measurement_sorting_algorithm.hpp"
  "src/clusterization/measurement_sorting_algorithm.sycl"
  # Seeding algorithm(s).
  "include/traccc/sycl/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"
  "src/seeding/silicon_pixel_spacepoint_formation_algorithm.sycl"
  "include/traccc/sycl/seeding/triplet_seeding_algorithm.hpp"
  "src/seeding/triplet_seeding_algorithm.sycl"
  "include/traccc/sycl/seeding/seed_parameter_estimation_algorithm.hpp"
  "src/seeding/seed_parameter_estimation_algorithm.sycl"
  # GBTS seeding algorithm(s).
  "include/traccc/sycl/gbts_seeding/gbts_seeding_algorithm.hpp"
  "src/gbts_seeding/gbts_seeding_algorithm.sycl"
  # Track finding algorithm(s).
  "include/traccc/sycl/finding/combinatorial_kalman_filter_algorithm.hpp"
  "src/finding/combinatorial_kalman_filter_algorithm.cpp"
  "src/finding/combinatorial_kalman_filter_algorithm.sycl"
  # Track fitting algorithm(s).
  "include/traccc/sycl/fitting/kalman_fitting_algorithm.hpp"
  "src/fitting/kalman_fitting_algorithm.cpp"
  "src/fitting/kalman_fitting_algorithm.sycl"
  # Utilities.
  "include/traccc/sycl/utils/algorithm_base.hpp"
  "src/utils/algorithm_base.sycl"
  "src/utils/magnetic_field_types.hpp"
  "src/utils/barrier.hpp"
  "src/utils/calculate1DimNdRange.hpp"
  "src/utils/calculate1DimNdRange.sycl"
  "src/utils/get_queue.hpp"
  "src/utils/get_queue.sycl"
  "src/utils/thread_id.hpp"
  "src/utils/oneDPL.hpp"
  "src/sanity/contiguous_on.hpp"
  "src/sanity/ordered_on.hpp"
)
target_link_libraries(
    traccc_sycl
    PUBLIC
        detray::core
        vecmem::core
        traccc::core
        traccc::device_common
        traccc::sycl_utils
    PRIVATE vecmem::sycl covfie::sycl oneDPL
)

# Set up Thrust specifically for the traccc::sycl library.
thrust_create_target( traccc::sycl_thrust
  HOST CPP
  DEVICE CPP
)
target_link_libraries(traccc_sycl PRIVATE traccc::sycl_thrust)
