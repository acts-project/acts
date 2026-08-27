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
    # Clusterization code.
    "include/traccc/hip/clusterization/clusterization_algorithm.hpp"
    "src/clusterization/clusterization_algorithm.hip"
    "src/clusterization/kernels/ccl_kernel.hpp"
    "src/clusterization/kernels/ccl_kernel.hip"
    "src/clusterization/kernels/reify_cluster_data.hpp"
    "src/clusterization/kernels/reify_cluster_data.hip"
    "include/traccc/hip/clusterization/measurement_sorting_algorithm.hpp"
    "src/clusterization/measurement_sorting_algorithm.hip"
)
target_include_directories(
    traccc_hip_device
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include"
)
target_link_libraries(
    traccc_hip_device
    PUBLIC traccc::core traccc::hip_utils vecmem::core
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
