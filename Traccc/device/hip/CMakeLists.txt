# TRACCC library, part of the ACTS project (R&D line)
#
# (c) 2025-2026 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Enable HIP as a language.
enable_language( HIP )

# Project include(s).
include( traccc-compiler-options-cpp )
include( traccc-compiler-options-hip )

# Find HIP/ROCm toolchain.
find_package( HIPToolkit REQUIRED )

# Set up the build of the traccc::hip library.
traccc_add_library( traccc_hip hip TYPE SHARED
   # Utility code.
   "include/traccc/hip/utils/algorithm_base.hpp"
   "src/utils/algorithm_base.cpp"
   "include/traccc/hip/utils/stream.hpp"
   "src/utils/stream.cpp"
   "src/utils/opaque_stream.hpp"
   "src/utils/opaque_stream.cpp"
   "src/utils/hip_error_handling.hpp"
   "src/utils/hip_error_handling.cpp"
   "src/utils/utils.hpp"
   "src/utils/utils.cpp"
   # Clusterization code.
   "include/traccc/hip/clusterization/clusterization_algorithm.hpp"
   "src/clusterization/clusterization_algorithm.hip"
   "src/clusterization/kernels/ccl_kernel.hpp"
   "src/clusterization/kernels/ccl_kernel.hip"
   "src/clusterization/kernels/reify_cluster_data.hpp"
   "src/clusterization/kernels/reify_cluster_data.hip"
)
target_link_libraries( traccc_hip
   PUBLIC traccc::core vecmem::core
   PRIVATE HIP::hiprt traccc::device_common )
