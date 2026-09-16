# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2025-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Enable HIP as a language.
enable_language(HIP)

# Project include(s).
include(traccc-compiler-options-cpp)
include(traccc-compiler-options-hip)

# Find HIP/ROCm toolchain.
find_package(HIPToolkit REQUIRED)

# Set up the build of the traccc::hip library. In a bit of a complicated way.
# All the device code is built into a separate object library, since rocThrust
# is not quite prepared for being linked against libraries that would hold both
# C++ and HIP source files.
add_library(
    traccc_hip_device
    OBJECT
    # Utility code.
    "src/utils/magnetic_field_types.hpp"
    # Clusterization code.
    "include/traccc/hip/clusterization/clusterization_algorithm.hpp"
    "src/clusterization/clusterization_algorithm.hip"
    "src/clusterization/kernels/ccl_kernel.hpp"
    "src/clusterization/kernels/ccl_kernel.hip"
    "src/clusterization/kernels/reify_cluster_data.hpp"
    "src/clusterization/kernels/reify_cluster_data.hip"
    "include/traccc/hip/clusterization/measurement_sorting_algorithm.hpp"
    "src/clusterization/measurement_sorting_algorithm.hip"
    # Seed finding code.
    "include/traccc/hip/seeding/seed_parameter_estimation_algorithm.hpp"
    "src/seeding/seed_parameter_estimation_algorithm.hip"
    "include/traccc/hip/seeding/triplet_seeding_algorithm.hpp"
    "src/seeding/triplet_seeding_algorithm.hip"
    "include/traccc/hip/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"
    "src/seeding/silicon_pixel_spacepoint_formation_algorithm.hip"
    # Track finding code.
    "include/traccc/hip/finding/combinatorial_kalman_filter_algorithm.hpp"
    "src/finding/combinatorial_kalman_filter_algorithm_build_measurement_ranges_buffer.hip"
    "src/finding/combinatorial_kalman_filter_algorithm_build_tracks_kernel.hip"
    "src/finding/combinatorial_kalman_filter_algorithm_condense_tracks_kernel.hip"
    "src/finding/combinatorial_kalman_filter_algorithm_create_device_detector.hip"
    "src/finding/combinatorial_kalman_filter_algorithm_fill_finding_duplicate_removal_sort_keys_kernel.hip"
    "src/finding/combinatorial_kalman_filter_algorithm_fill_finding_propagation_sort_keys_kernel.hip"
    "src/finding/combinatorial_kalman_filter_algorithm_find_tracks_kernel.hip"
    "src/finding/combinatorial_kalman_filter_algorithm_gather_best_tips_per_measurement_kernel.hip"
    "src/finding/combinatorial_kalman_filter_algorithm_gather_measurement_votes_kernel.hip"
    "src/finding/combinatorial_kalman_filter_algorithm_input_is_valid.hip"
    "src/finding/combinatorial_kalman_filter_algorithm_progressive_kalman_filter_kernel.hip"
    "src/finding/combinatorial_kalman_filter_algorithm_propagate_to_next_surface_kernel.hip"
    "src/finding/combinatorial_kalman_filter_algorithm_remove_duplicates_kernel.hip"
    "src/finding/combinatorial_kalman_filter_algorithm_sort_param_ids_by_keys.hip"
    "src/finding/combinatorial_kalman_filter_algorithm_sort_param_ids_by_last_measurement.hip"
    "src/finding/combinatorial_kalman_filter_algorithm_update_tip_length_buffer_kernel.hip"
    # Track fitting code.
    "include/traccc/hip/fitting/kalman_fitting_algorithm.hpp"
    "src/fitting/kalman_fitting_algorithm_fit_backward_kernel.hip"
    "src/fitting/kalman_fitting_algorithm_fit_forward_kernel.hip"
    "src/fitting/kalman_fitting_algorithm_fit_prelude_kernel.hip"
    "src/fitting/kalman_fitting_algorithm_prepare_track_fit_order.hip"
)
target_include_directories(
    traccc_hip_device
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_link_libraries(
    traccc_hip_device
    PUBLIC traccc::core traccc::hip_utils vecmem::core covfie::core covfie::hip
    PRIVATE HIP::hiprt traccc::device_common
)
set_target_properties(
    traccc_hip_device
    PROPERTIES POSITION_INDEPENDENT_CODE TRUE
)

# Figure out how to use Thrust with the library. Unfortunately modern rocThrust
# versions seem to have lost the ability to use a CUDA backend reliably.
if("${CMAKE_HIP_PLATFORM}" STREQUAL "amd")
    target_link_libraries(traccc_hip_device PRIVATE roc::rocthrust)
elseif("${CMAKE_HIP_PLATFORM}" STREQUAL "nvidia")
    thrust_create_target(traccc::hip_device_thrust
        HOST CPP
        DEVICE CUDA
    )
    target_link_libraries(traccc_hip_device PRIVATE traccc::hip_device_thrust)
endif()

traccc_add_library( traccc_hip hip TYPE SHARED
   # Track finding code.
   "src/finding/combinatorial_kalman_filter_algorithm.cpp"
   # Track fitting code.
   "src/fitting/kalman_fitting_algorithm.cpp"
   # Utility code.
   "include/traccc/hip/utils/algorithm_base.hpp"
   "src/utils/algorithm_base.cpp"
   "src/utils/barrier.hpp"
   "src/utils/global_index.hpp"
   "src/utils/hip_error_handling.hpp"
   "src/utils/hip_error_handling.cpp"
   "src/utils/utils.hpp"
   "src/utils/utils.cpp"
   # Sanity checking code.
   "src/sanity/contiguous_on.hpp"
   "src/sanity/ordered_on.hpp"
)
target_link_libraries(
    traccc_hip
    PUBLIC traccc::core traccc::hip_utils vecmem::core
    PRIVATE HIP::hiprt traccc::device_common traccc_hip_device
)
